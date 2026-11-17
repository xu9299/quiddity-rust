//! Positive rational Bezier hulls of a finite NURBS interval. Subdivision only tightens
//! enclosing control hulls; sampled points select work, never replace an outer bound.

use crate::geom::V3;
use crate::nurbs::NurbsCurve;

type H = [f64; 3];
type Point = (f64, f64);

fn mix<const N: usize>(a: [f64; N], b: [f64; N], t: f64) -> [f64; N] {
    std::array::from_fn(|i| (1.0 - t) * a[i] + t * b[i])
}

fn point(p: H) -> Option<Point> {
    (p[2] > 0.0 && p.into_iter().all(f64::is_finite))
        .then_some((p[0] / p[2], p[1] / p[2]))
        .filter(|p| p.0.is_finite() && p.1.is_finite())
}

pub(super) fn split<const N: usize>(net: &[[f64; N]], t: f64) -> (Vec<[f64; N]>, Vec<[f64; N]>) {
    let mut work = net.to_vec();
    let mut left = vec![work[0]];
    let mut right = vec![work[work.len() - 1]];
    for level in 1..net.len() {
        for i in 0..net.len() - level {
            work[i] = mix(work[i], work[i + 1], t);
        }
        left.push(work[0]);
        right.push(work[net.len() - level - 1]);
    }
    right.reverse();
    (left, right)
}

fn insert<const N: usize>(
    points: &mut Vec<[f64; N]>,
    knots: &mut Vec<f64>,
    degree: usize,
    t: f64,
) -> Option<()> {
    let multiplicity = knots.iter().filter(|&&u| u == t).count();
    if multiplicity > degree {
        return Some(());
    }
    let k = knots.iter().rposition(|&u| u <= t)?;
    if k < degree || k < multiplicity || k - multiplicity >= points.len() {
        return None;
    }
    let mut out = vec![[0.0; N]; points.len() + 1];
    out[..=k - degree].copy_from_slice(&points[..=k - degree]);
    out[k - multiplicity + 1..].copy_from_slice(&points[k - multiplicity..]);
    for i in k - degree + 1..=k - multiplicity {
        let denominator = knots[i + degree] - knots[i];
        let alpha = (t - knots[i]) / denominator;
        if !alpha.is_finite() || !(0.0..=1.0).contains(&alpha) {
            return None;
        }
        out[i] = mix(points[i - 1], points[i], alpha);
    }
    if out
        .iter()
        .any(|p| p[N - 1] <= 0.0 || p.iter().any(|x| !x.is_finite()))
    {
        return None;
    }
    knots.insert(k + 1, t);
    *points = out;
    Some(())
}

fn hull(net: &[H]) -> Option<(Point, Point)> {
    let mut low = (f64::INFINITY, f64::INFINITY);
    let mut high = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &p in net {
        let p = point(p)?;
        low = (low.0.min(p.0), low.1.min(p.1));
        high = (high.0.max(p.0), high.1.max(p.1));
    }
    Some((low, high))
}

/// Return homogeneous Bezier nets in the finite interval's directed path order. The final
/// coordinate of `project` is one; it becomes the normalized positive rational weight.
/// A wrapped closed interval retains its starting phase and every traversed domain piece.
pub(super) fn beziers<const N: usize>(
    curve: &NurbsCurve,
    interval: (f64, f64),
    project: impl Fn(V3) -> [f64; N],
) -> Option<Vec<Vec<[f64; N]>>> {
    let degree = curve.degree();
    if degree == 0
        || curve.control_points().len() > 4096
        || !curve.weights().iter().all(|w| w.is_finite() && *w > 0.0)
    {
        return None;
    }
    let maximum_weight = curve.weights().iter().copied().fold(0.0, f64::max);
    let mut points: Vec<_> = curve
        .control_points()
        .iter()
        .zip(curve.weights())
        .map(|(&p, &w)| {
            let q = project(p);
            q.map(|x| x * (w / maximum_weight))
        })
        .collect();
    if points
        .iter()
        .any(|p| p[N - 1] <= 0.0 || p.iter().any(|x| !x.is_finite()))
    {
        return None;
    }
    let (lo, hi) = curve.domain();
    let (a, b) = (interval.0.min(interval.1), interval.0.max(interval.1));
    if !a.is_finite() || !b.is_finite() || hi <= lo {
        return None;
    }
    let ranges = if a >= lo && b <= hi {
        vec![(a, b)]
    } else if curve.is_closed() {
        let period = hi - lo;
        if b - a > 2.0 * period || (a - lo).abs() > 4.0 * period {
            return None;
        }
        let mut start = lo + (a - lo).rem_euclid(period);
        let mut remaining = b - a;
        let mut ranges = Vec::new();
        for _ in 0..3 {
            let length = remaining.min(hi - start);
            ranges.push((start, start + length));
            remaining -= length;
            if remaining <= 0.0 {
                break;
            }
            start = lo;
        }
        if remaining > 0.0 {
            return None;
        }
        ranges
    } else {
        return None;
    };
    let mut knots = curve.knots().to_vec();
    let mut cuts: Vec<_> = knots
        .iter()
        .copied()
        .filter(|u| *u >= lo && *u <= hi)
        .collect();
    cuts.dedup();
    for t in cuts {
        let required = if t == lo || t == hi {
            degree + 1
        } else {
            degree
        };
        while knots.iter().filter(|&&u| u == t).count() < required {
            insert(&mut points, &mut knots, degree, t)?;
            if points.len() > 16384 {
                return None;
            }
        }
    }
    let mut nets = Vec::new();
    for &(a, b) in &ranges {
        for k in degree..points.len() {
            let (u, v) = (knots[k], knots[k + 1]);
            if u >= v || u < lo || v > hi {
                continue;
            }
            let (a, b) = (a.max(u), b.min(v));
            if a > b || (a == b && interval.0 != interval.1) {
                continue;
            }
            let mut net = points[k - degree..=k].to_vec();
            let (alpha, beta) = ((a - u) / (v - u), (b - u) / (v - u));
            if beta < 1.0 {
                net = split(&net, beta).0;
            }
            if alpha > 0.0 && beta > 0.0 {
                net = split(&net, alpha / beta).1;
            }
            nets.push(net);
        }
    }
    if interval.1 < interval.0 {
        nets.reverse();
        for net in &mut nets {
            net.reverse();
        }
    }
    (!nets.is_empty()).then_some(nets)
}

pub(super) fn bounds(
    curve: &NurbsCurve,
    interval: (f64, f64),
    project: impl Fn(V3) -> Point,
    accuracy: f64,
) -> Option<(Point, Point)> {
    let degree = curve.degree();
    let mut queue: Vec<_> = beziers(curve, interval, |p| {
        let q = project(p);
        [q.0, q.1, 1.0]
    })?
    .into_iter()
    .map(|net| (net, 0))
    .collect();
    let size = curve
        .control_points()
        .iter()
        .copied()
        .map(project)
        .flat_map(|p| [p.0.abs(), p.1.abs()])
        .fold(1.0, f64::max);
    if !size.is_finite() {
        return None;
    }
    let mut known_low = (f64::INFINITY, f64::INFINITY);
    let mut known_high = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut observe = |p: H| -> Option<()> {
        let p = point(p)?;
        known_low = (known_low.0.min(p.0), known_low.1.min(p.1));
        known_high = (known_high.0.max(p.0), known_high.1.max(p.1));
        Some(())
    };
    for (net, _) in &queue {
        observe(net[0])?;
        observe(*net.last()?)?;
    }
    let mut low = (f64::INFINITY, f64::INFINITY);
    let mut high = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let accuracy = accuracy.max(64.0 * f64::EPSILON * size);
    let mut subdivisions = 0;
    while let Some((net, depth)) = queue.pop() {
        let (a, b) = hull(&net)?;
        let unresolved = a.0 < known_low.0 - accuracy
            || a.1 < known_low.1 - accuracy
            || b.0 > known_high.0 + accuracy
            || b.1 > known_high.1 + accuracy;
        if unresolved && depth < 32 && subdivisions < 4096 {
            let (left, right) = split(&net, 0.5);
            let p = point(*left.last()?)?;
            known_low = (known_low.0.min(p.0), known_low.1.min(p.1));
            known_high = (known_high.0.max(p.0), known_high.1.max(p.1));
            queue.push((left, depth + 1));
            queue.push((right, depth + 1));
            subdivisions += 1;
        } else {
            // Exhaustion retains this entire hull. It never accepts sampled extrema.
            low = (low.0.min(a.0), low.1.min(a.1));
            high = (high.0.max(b.0), high.1.max(b.1));
        }
    }
    let pad = 256.0 * f64::EPSILON * size * (degree + 1).pow(2) as f64 * 34.0;
    Some(((low.0 - pad, low.1 - pad), (high.0 + pad, high.1 + pad)))
}
