//! Finite rational/conic correspondence without assuming equal parameter speeds. Positive
//! Bezier nets bound the conic's implicit residual. Bernstein derivative signs and an open
//! radial hemisphere prove each piece's angular order and winding. Reparameterizing by
//! angle then bounds the complete finite paths, including their endpoint phase error.

use super::spline_bounds;
use crate::geom::{self, Curve, V3};

type H = [f64; 4];

fn axes(curve: &Curve) -> Option<(V3, V3, V3)> {
    match curve {
        Curve::Circle { frame, radius } if *radius > 0.0 => Some((
            frame.origin,
            geom::scale(frame.x, *radius),
            geom::scale(frame.y, *radius),
        )),
        Curve::Ellipse {
            frame,
            major,
            minor,
        } if *major > 0.0 && *minor > 0.0 => Some((
            frame.origin,
            geom::scale(frame.x, *major),
            geom::scale(frame.y, *minor),
        )),
        _ => None,
    }
}

fn choose(n: usize, k: usize) -> f64 {
    (0..k.min(n - k)).fold(1.0, |value, i| value * (n - i) as f64 / (i + 1) as f64)
}

fn angle(a: H, b: H) -> f64 {
    (a[0] * b[1] - a[1] * b[0]).atan2(a[0] * b[0] + a[1] * b[1])
}

struct Piece {
    radial: f64,
    normal: f64,
    span: f64,
}

fn piece(net: &[H], sense: f64) -> Option<Piece> {
    let p = net.len().checked_sub(1)?;
    if p == 0 {
        return None;
    }
    let maximum_weight = net.iter().map(|h| h[3]).fold(0.0, f64::max);
    let net: Vec<_> = net.iter().map(|h| h.map(|v| v / maximum_weight)).collect();
    let minimum_weight = net.iter().map(|h| h[3]).fold(1.0, f64::min);
    if minimum_weight <= 1e-8 || net.iter().flatten().any(|x| !x.is_finite()) {
        return None;
    }
    let magnitude = net
        .iter()
        .flat_map(|h| [h[0].abs(), h[1].abs()])
        .fold(1.0, f64::max);
    let round = 256.0 * f64::EPSILON * (p + 1).pow(2) as f64 * magnitude * magnitude;
    let (a, b) = (net[0], *net.last()?);
    let na = a[0].hypot(a[1]);
    let nb = b[0].hypot(b[1]);
    if na <= round || nb <= round {
        return None;
    }
    let direction = [a[0] / na + b[0] / nb, a[1] / na + b[1] / nb];
    if net
        .iter()
        .any(|h| h[0] * direction[0] + h[1] * direction[1] <= round)
    {
        return None;
    }
    // The product of Bernstein polynomials has nonnegative combination coefficients.
    // Its control coefficients bound the entire implicit residual and derivative sign.
    let mut residual = 0.0_f64;
    for k in 0..=2 * p {
        let mut coefficient = 0.0;
        for i in k.saturating_sub(p)..=k.min(p) {
            let j = k - i;
            let factor = choose(p, i) * choose(p, j) / choose(2 * p, k);
            coefficient +=
                factor * (net[i][0] * net[j][0] + net[i][1] * net[j][1] - net[i][3] * net[j][3]);
        }
        residual = residual.max(coefficient.abs());
    }
    for k in 0..2 * p {
        let mut coefficient = 0.0;
        for i in k.saturating_sub(p - 1)..=k.min(p) {
            let j = k - i;
            let dx = p as f64 * (net[j + 1][0] - net[j][0]);
            let dy = p as f64 * (net[j + 1][1] - net[j][1]);
            let factor = choose(p, i) * choose(p - 1, j) / choose(2 * p - 1, k);
            coefficient += factor * (net[i][0] * dy - net[i][1] * dx);
        }
        if sense * coefficient <= round * (p + 1) as f64 {
            return None;
        }
    }
    let residual = (residual + round) / minimum_weight.powi(2);
    if residual >= 0.5 {
        return None;
    }
    let span = angle(a, b);
    if sense * span <= 0.0 || !span.is_finite() {
        return None;
    }
    Some(Piece {
        radial: residual / (1.0 + (1.0 - residual).sqrt()),
        normal: net.iter().map(|h| (h[2] / h[3]).abs()).fold(0.0, f64::max),
        span,
    })
}

/// A complete-path distance bound under angular reparameterization. Unsupported or
/// unproved angular order returns `None`; a large bound never implies acceptance.
pub(super) fn deviation(
    reference: &Curve,
    interval: (f64, f64),
    other: &Curve,
    other_interval: (f64, f64),
) -> Option<f64> {
    let (origin, x, y) = axes(reference)?;
    let Curve::Nurbs(curve) = other else {
        return None;
    };
    let (xx, yy, xy) = (geom::dot(x, x), geom::dot(y, y), geom::dot(x, y));
    let determinant = xx * yy - xy * xy;
    if determinant <= 0.25 * xx * yy || !determinant.is_finite() {
        return None;
    }
    let maximum = (0.5 * (xx + yy + ((xx - yy).powi(2) + 4.0 * xy * xy).sqrt())).sqrt();
    let minimum = determinant.sqrt() / maximum;
    let normal = geom::unit(geom::cross(x, y))?;
    let project = |p| {
        let q = geom::sub(p, origin);
        let (a, b) = (geom::dot(q, x), geom::dot(q, y));
        [
            (a * yy - b * xy) / determinant,
            (b * xx - a * xy) / determinant,
            geom::dot(q, normal),
            1.0,
        ]
    };
    let nets = spline_bounds::beziers(curve, other_interval, project)?;
    let first = *nets.first()?.first()?;
    let span = interval.1 - interval.0;
    if !span.is_finite() || span == 0.0 || span.abs() > 2.0 * std::f64::consts::TAU {
        return None;
    }
    let sense = span.signum();
    let mut queue: Vec<_> = nets.into_iter().rev().map(|n| (n, 0)).collect();
    let (mut radial, mut off_plane, mut actual_span) = (0.0_f64, 0.0_f64, 0.0);
    let (mut joins, mut last) = (0.0_f64, None::<H>);
    let mut count = 0;
    while let Some((net, depth)) = queue.pop() {
        let Some(bound) = piece(&net, sense) else {
            if depth >= 12 || count + queue.len() >= 1024 {
                return None;
            }
            let (left, right) = spline_bounds::split(&net, 0.5);
            queue.push((right, depth + 1));
            queue.push((left, depth + 1));
            continue;
        };
        if let Some(previous) = last {
            let jump = angle(previous, net[0]).abs();
            if !jump.is_finite() || jump > 0.01 {
                return None;
            }
            // At a periodic join, finite endpoint roundoff can leave a tiny connector.
            // Its angular drift is charged instead of assuming exact closure.
            joins += jump;
        }
        last = net.last().copied();
        radial = radial.max(bound.radial);
        off_plane = off_plane.max(bound.normal);
        actual_span += bound.span;
        count += 1;
    }
    let phase = (first[1].atan2(first[0]) - interval.0)
        .sin()
        .atan2((first[1].atan2(first[0]) - interval.0).cos());
    let angular = phase.abs().max((phase + actual_span - span).abs()) + joins;
    let size = curve
        .control_points()
        .iter()
        .flatten()
        .chain(origin.iter())
        .map(|x| x.abs())
        .fold(1.0, f64::max);
    // Include affine projection, knot insertion, Bernstein arithmetic and angle summation.
    // The matrix guard above keeps the inverse coordinate map uniformly conditioned.
    let round = 1024.0
        * f64::EPSILON
        * (curve.degree() + 1).pow(2) as f64
        * (count + 1) as f64
        * size
        * (1.0 + maximum / minimum);
    let bound = off_plane + maximum * (radial + angular) + round;
    (bound.is_finite() && bound >= 0.0).then_some(bound)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Frame;
    use crate::nurbs::NurbsCurve;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    fn circle() -> Curve {
        Curve::Circle {
            frame: Frame {
                origin: [0.0; 3],
                x: [1.0, 0.0, 0.0],
                y: [0.0, 1.0, 0.0],
                z: [0.0, 0.0, 1.0],
            },
            radius: 1.0,
        }
    }

    fn periodic(scale: f64) -> Curve {
        let points = (0..7)
            .map(|i| {
                let a = i as f64 * PI / 3.0;
                let r = if i % 2 == 0 { scale } else { 2.0 * scale };
                [r * a.cos(), r * a.sin(), 0.0]
            })
            .collect();
        Curve::Nurbs(
            NurbsCurve::new(
                2,
                points,
                vec![1.0, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0],
                vec![
                    -TAU / 3.0,
                    0.0,
                    0.0,
                    TAU / 3.0,
                    TAU / 3.0,
                    2.0 * TAU / 3.0,
                    2.0 * TAU / 3.0,
                    TAU,
                    TAU,
                    4.0 * TAU / 3.0,
                ],
            )
            .unwrap(),
        )
    }

    #[test]
    fn periodic_rational_conics_retain_finite_phase_direction_and_winding() {
        let c = circle();
        let n = periodic(1.0);
        assert!(deviation(&c, (0.0, TAU), &n, (0.0, TAU)).unwrap() < 1e-9);
        assert!(deviation(&c, (TAU, 0.0), &n, (TAU, 0.0)).unwrap() < 1e-9);
        // A wrapped whole interval starts at the corresponding geometric phase.
        let u = 0.3;
        let a = c.parameter(n.value(u));
        assert!(deviation(&c, (a, a + TAU), &n, (u, u + TAU)).unwrap() < 1e-9);
        assert!(deviation(&c, (0.0, TAU), &n, (TAU, 0.0)).is_none());
        assert!(deviation(&c, (0.0, PI), &n, (0.0, TAU)).unwrap() > 3.0);
        assert!(deviation(&c, (0.0, TAU), &periodic(1.01), (0.0, TAU)).unwrap() >= 0.01);
    }

    #[test]
    fn finite_quarter_arc_does_not_assume_equal_rational_parameter_speed() {
        let n = Curve::Nurbs(
            NurbsCurve::new(
                2,
                vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                vec![1.0, 0.5_f64.sqrt(), 1.0],
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            )
            .unwrap(),
        );
        let c = circle();
        for (a, b) in [(0.0, 1.0), (0.1, 0.7), (0.7, 0.1)] {
            let interval = (c.parameter(n.value(a)), c.parameter(n.value(b)));
            let d = deviation(&c, interval, &n, (a, b)).unwrap();
            assert!(d < 1e-9, "{a}..{b}: {d}");
        }
        assert!(deviation(&c, (0.0, -3.0 * FRAC_PI_2), &n, (0.0, 1.0)).is_none());
    }
}
