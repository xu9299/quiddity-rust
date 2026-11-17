//! Positive rational Bezier hulls of a finite NURBS interval. Subdivision only tightens
//! enclosing control hulls; sampled points select work, never replace an outer bound.

use crate::geom::V3;
use crate::nurbs::NurbsCurve;

type H = [f64; 3];
type Point = (f64, f64);

fn mix(a: H, b: H, t: f64) -> H {
    std::array::from_fn(|i| (1.0 - t) * a[i] + t * b[i])
}

fn point(p: H) -> Option<Point> {
    (p[2] > 0.0 && p.into_iter().all(f64::is_finite))
        .then_some((p[0] / p[2], p[1] / p[2]))
        .filter(|p| p.0.is_finite() && p.1.is_finite())
}

fn split(net: &[H], t: f64) -> (Vec<H>, Vec<H>) {
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

fn insert(points: &mut Vec<H>, knots: &mut Vec<f64>, degree: usize, t: f64) -> Option<()> {
    let multiplicity = knots.iter().filter(|&&u| u == t).count();
    if multiplicity > degree {
        return Some(());
    }
    let k = knots.iter().rposition(|&u| u <= t)?;
    if k < degree || k < multiplicity || k - multiplicity >= points.len() {
        return None;
    }
    let mut out = vec![[0.0; 3]; points.len() + 1];
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
    if out.iter().copied().any(|p| point(p).is_none()) {
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

pub(super) fn bounds(
    curve: &NurbsCurve,
    interval: (f64, f64),
    project: impl Fn(V3) -> Point,
    accuracy: f64,
) -> Option<(Point, Point)> {
    let degree = curve.degree();
    if degree == 0
        || curve.control_points().len() > 4096
        || !curve.weights().iter().all(|w| w.is_finite() && *w > 0.0)
    {
        return None;
    }
    let mut points: Vec<_> = curve
        .control_points()
        .iter()
        .zip(curve.weights())
        .map(|(&p, &w)| {
            let q = project(p);
            [q.0 * w, q.1 * w, w]
        })
        .collect();
    let (original_low, original_high) = hull(&points)?;
    let size = [
        original_low.0,
        original_low.1,
        original_high.0,
        original_high.1,
    ]
    .into_iter()
    .map(f64::abs)
    .fold(1.0, f64::max);
    let (lo, hi) = curve.domain();
    let (a, b) = (interval.0.min(interval.1), interval.0.max(interval.1));
    if !a.is_finite() || !b.is_finite() || hi <= lo {
        return None;
    }
    let ranges = if a >= lo && b <= hi {
        vec![(a, b)]
    } else if curve.is_closed() {
        let period = hi - lo;
        if b - a >= period * (1.0 - 64.0 * f64::EPSILON) {
            vec![(lo, hi)]
        } else {
            let start = lo + (a - lo).rem_euclid(period);
            let end = start + (b - a);
            if end <= hi {
                vec![(start, end)]
            } else {
                vec![(start, hi), (lo, (lo + end - hi).min(hi))]
            }
        }
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
    let mut queue = Vec::new();
    for k in degree..points.len() {
        let (u, v) = (knots[k], knots[k + 1]);
        if u >= v || u < lo || v > hi {
            continue;
        }
        for &(a, b) in &ranges {
            let (a, b) = (a.max(u), b.min(v));
            if a > b {
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
            queue.push((net, 0));
        }
    }
    if queue.is_empty() {
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
