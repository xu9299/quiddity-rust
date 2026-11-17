//! Finite polynomial cylinder pcurves projected onto a held plane.
//!
//! The plane/cylinder intersection supplies an analytic reference. At a fixed angular
//! parameter, the projected pcurve differs from it in the projected axial direction.
//! Its homotopy Jacobian is therefore the axial residual times the intersection's
//! transverse derivative; the quadratic homotopy term vanishes. Positive Bernstein
//! hulls and Taylor remainders bound that complete finite integral. This is not a fit,
//! a sampled acceptance test, or an allowance based on the file uncertainty.

use crate::geom::{self, Curve, Frame, Surface, V3};

const DEGREE: usize = 8;

fn split(net: &[f64], t: f64) -> (Vec<f64>, Vec<f64>) {
    let mut net = net.to_vec();
    let (mut left, mut right) = (vec![net[0]], vec![*net.last().unwrap()]);
    while net.len() > 1 {
        net = net
            .windows(2)
            .map(|w| (1.0 - t) * w[0] + t * w[1])
            .collect();
        left.push(net[0]);
        right.push(*net.last().unwrap());
    }
    right.reverse();
    (left, right)
}

fn restrict(net: &[f64], a: f64, b: f64) -> Vec<f64> {
    let left = split(net, b).0;
    split(&left, a / b).1
}

fn elevate(mut net: Vec<f64>) -> Vec<f64> {
    while net.len() <= DEGREE {
        let n = net.len();
        let mut out = vec![net[0]];
        out.extend((1..n).map(|i| {
            let t = i as f64 / n as f64;
            t * net[i - 1] + (1.0 - t) * net[i]
        }));
        out.push(*net.last().unwrap());
        net = out;
    }
    net
}

fn choose(n: usize, k: usize) -> f64 {
    (0..k.min(n - k)).fold(1.0, |v, i| v * (n - i) as f64 / (i + 1) as f64)
}

fn trig_net(a: f64, b: f64, cosine: bool) -> (Vec<f64>, f64) {
    let middle = 0.5 * (a + b);
    let (sin, cos) = middle.sin_cos();
    let derivatives = if cosine {
        [cos, -sin, -cos, sin]
    } else {
        [sin, cos, -sin, -cos]
    };
    let mut powers = [0.0; DEGREE + 1];
    let mut coefficient = 1.0;
    for k in 0..=DEGREE {
        if k != 0 {
            coefficient *= (b - a) / k as f64;
        }
        for (j, value) in powers.iter_mut().enumerate().take(k + 1) {
            *value +=
                derivatives[k % 4] * coefficient * choose(k, j) * (-0.5_f64).powi((k - j) as i32);
        }
    }
    let net = (0..=DEGREE)
        .map(|i| {
            (0..=i)
                .map(|j| powers[j] * choose(i, j) / choose(DEGREE, j))
                .sum()
        })
        .collect();
    let remainder = (0.5 * (b - a).abs()).powi(9) / 362880.0;
    let round = 512.0 * f64::EPSILON * (DEGREE + 1) as f64 * (1.0 + middle.abs());
    (net, remainder + round)
}

/// Return a projected area charge and a distance from the nominal finite conic to the
/// analytic intersection. The caller separately charges that latter complete-path tube.
pub(super) fn planar_error(
    reference: &Curve,
    interval: (f64, f64),
    endpoints: [V3; 2],
    plane: &Frame,
    cylinder: (&Frame, f64),
    points: &[[f64; 2]],
    eligibility: f64,
) -> Option<(f64, f64)> {
    let (frame, radius) = cylinder;
    let degree = points.len().checked_sub(1)?;
    if !(2..=DEGREE).contains(&degree)
        || radius <= 0.0
        || !radius.is_finite()
        || points.iter().flatten().any(|x| !x.is_finite())
        || points.iter().any(|p| p[0].abs() > 64.0)
        || ![interval.0, interval.1].into_iter().all(f64::is_finite)
        || (interval.1 - interval.0).abs() >= std::f64::consts::PI
    {
        return None;
    }
    let (origin, x, y) = super::corresponding::conic(reference)?;
    if endpoints
        .iter()
        .flatten()
        .chain(&origin)
        .chain(&x)
        .chain(&y)
        .chain(&plane.origin)
        .chain(&plane.z)
        .chain(&frame.origin)
        .chain(&frame.x)
        .chain(&frame.y)
        .chain(&frame.z)
        .any(|x| !x.is_finite())
        || [plane.z, frame.x, frame.y, frame.z]
            .into_iter()
            .any(|v| !(0.5..=2.0).contains(&geom::norm(v)))
        || eligibility.is_nan()
        || eligibility < 0.0
    {
        return None;
    }
    let normal = geom::unit(plane.z)?;
    let denominator = geom::dot(normal, frame.z);
    if denominator.abs() < 0.25 {
        return None;
    }
    let u: Vec<_> = points.iter().map(|p| p[0]).collect();
    let v: Vec<_> = points.iter().map(|p| p[1]).collect();
    let width = u[degree] - u[0];
    if width == 0.0
        || width.abs() >= std::f64::consts::PI
        || u.windows(2).any(|w| (w[1] - w[0]) * width <= 0.0)
    {
        return None;
    }
    let alpha = geom::dot(x, frame.y).atan2(geom::dot(x, frame.x));
    let determinant = geom::dot(x, frame.x) * geom::dot(y, frame.y)
        - geom::dot(x, frame.y) * geom::dot(y, frame.x);
    if !determinant.is_finite() || determinant == 0.0 {
        return None;
    }
    let sense = determinant.signum();
    let start_u = geom::nearest_turn(alpha + sense * interval.0, 0.5 * (u[0] + u[degree]));
    let end_u = start_u + sense * (interval.1 - interval.0);
    let ta = ((start_u - u[0]) / width).clamp(0.0, 1.0);
    let tb = ((end_u - u[0]) / width).clamp(0.0, 1.0);
    if (tb - ta) * (interval.1 - interval.0) <= 0.0 {
        return None;
    }
    // Retain a finite spline end when the nominal conic extends beyond it. The
    // complete angular phase difference is charged, rather than rounded to zero.
    let phase_error = (start_u - (u[0] + width * ta))
        .abs()
        .max((end_u - (u[0] + width * tb)).abs());
    let (low, high) = (ta.min(tb), ta.max(tb));
    let magnitude = endpoints
        .iter()
        .flatten()
        .chain(origin.iter())
        .chain(x.iter())
        .chain(y.iter())
        .chain(plane.origin.iter())
        .chain(frame.origin.iter())
        .chain(points.iter().flatten())
        .map(|x| x.abs())
        .fold(radius.max(1.0), f64::max);
    // Two de Casteljau restrictions and degree elevation have linear affine depth.
    // Frame transforms and the guarded plane division add fixed depth. With the
    // bounded frame axes above, 64*(degree+1) error units enclose these operations;
    // the factor two includes their coordinate sums and conversion to a distance.
    let operations = 64.0 * (degree + 1) as f64;
    let gamma = operations * f64::EPSILON / (1.0 - operations * f64::EPSILON);
    let round = 2.0 * gamma * magnitude / denominator.abs();
    let surface = Surface::Cylinder {
        frame: *frame,
        radius,
    };
    let evaluate = |t| {
        let a = split(&u, t).0[degree];
        let b = split(&v, t).0[degree];
        surface.value(a, b)
    };
    let endpoint_error = [ta, tb]
        .into_iter()
        .zip(endpoints)
        .map(|(t, p)| geom::dist(evaluate(t), p))
        .fold(0.0_f64, f64::max)
        + round;
    if endpoint_error > eligibility {
        return None;
    }
    let constant = geom::dot(normal, geom::sub(frame.origin, plane.origin)) / denominator;
    let cosine = radius * geom::dot(normal, frame.x) / denominator;
    let sine = radius * geom::dot(normal, frame.y) / denominator;
    let ideal_origin = geom::sub(frame.origin, geom::scale(frame.z, constant));
    let ideal_x = geom::sub(geom::scale(frame.x, radius), geom::scale(frame.z, cosine));
    let ideal_y = geom::sub(geom::scale(frame.y, radius), geom::scale(frame.z, sine));
    let (sa, ca) = alpha.sin_cos();
    let nominal_x = geom::sub(geom::scale(x, ca), geom::scale(y, sense * sa));
    let nominal_y = geom::add(geom::scale(x, sa), geom::scale(y, sense * ca));
    let projected_length = |p| geom::norm(geom::sub(p, geom::scale(normal, geom::dot(p, normal))));
    let conic_distance = projected_length(geom::sub(origin, ideal_origin))
        + projected_length(geom::sub(nominal_x, ideal_x))
        + projected_length(geom::sub(nominal_y, ideal_y))
        + (geom::norm(x) + geom::norm(y)) * phase_error
        + round;
    let factor_sine = -radius * geom::dot(normal, geom::cross(frame.z, frame.x));
    let factor_cosine = radius * geom::dot(normal, geom::cross(frame.z, frame.y));
    let mut area = 0.0;
    for i in 0..64 {
        let a = low + (high - low) * i as f64 / 64.0;
        let b = low + (high - low) * (i + 1) as f64 / 64.0;
        let (ua, ub) = (u[0] + width * a, u[0] + width * b);
        let (cos, cos_error) = trig_net(ua, ub, true);
        let (sin, sin_error) = trig_net(ua, ub, false);
        let z = elevate(restrict(&v, a, b));
        let residual = z
            .iter()
            .zip(&cos)
            .zip(&sin)
            .map(|((&z, &c), &s)| (z + constant + cosine * c + sine * s).abs())
            .fold(0.0_f64, f64::max)
            + cosine.abs() * cos_error
            + sine.abs() * sin_error
            + round;
        let jacobian = cos
            .iter()
            .zip(&sin)
            .map(|(&c, &s)| (factor_cosine * c + factor_sine * s).abs())
            .fold(0.0_f64, f64::max)
            + factor_cosine.abs() * cos_error
            + factor_sine.abs() * sin_error
            + round;
        area += residual * jacobian * (ub - ua).abs();
    }
    let angular_error = u
        .iter()
        .enumerate()
        .map(|(i, &p)| (p - (u[0] + width * i as f64 / degree as f64)).abs())
        .fold(0.0_f64, f64::max)
        + 64.0 * f64::EPSILON * u.iter().map(|x| x.abs()).fold(1.0_f64, f64::max);
    let radial_norm = geom::norm(frame.x) + geom::norm(frame.y);
    let displacement = radius * radial_norm * angular_error + round;
    let derivative = |net: &[f64]| {
        net.windows(2)
            .map(|w| degree as f64 * (w[1] - w[0]).abs())
            .fold(0.0_f64, f64::max)
    };
    let length = (high - low)
        * (radius * radial_norm * derivative(&u) + geom::norm(frame.z) * derivative(&v))
        + round;
    area += 2.0 * displacement * length
        + std::f64::consts::PI * displacement * displacement
        + 16.0 * endpoint_error * endpoint_error
        + round * (1.0 + length);
    (area.is_finite() && conic_distance.is_finite()).then_some((area, conic_distance))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_restriction_and_trig_enclosures_keep_the_complete_interval() {
        assert_eq!(
            restrict(&[1.0, 2.0, 4.0], 0.25, 0.75),
            [25.0 / 16.0, 35.0 / 16.0, 49.0 / 16.0]
        );
        for cosine in [false, true] {
            let (net, error) = trig_net(0.0, 1.0, cosine);
            let exact = if cosine { 1.0 } else { 0.0 };
            assert!((net[0] - exact).abs() <= error);
            let exact = if cosine { 1.0_f64.cos() } else { 1.0_f64.sin() };
            assert!((net[DEGREE] - exact).abs() <= error);
        }
    }

    #[test]
    fn an_interior_bulge_keeps_its_area_charge_despite_equal_endpoints() {
        let axis = Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        };
        let h = 0.5_f64.sqrt();
        let plane = Frame {
            origin: [0.0; 3],
            x: [h, 0.0, h],
            y: [0.0, 1.0, 0.0],
            z: [-h, 0.0, h],
        };
        let reference = Curve::Ellipse {
            frame: plane,
            major: 2.0_f64.sqrt(),
            minor: 1.0,
        };
        let a = 0.01_f64;
        let endpoints = [reference.value(-a), reference.value(a)];
        let mut points = [[-a, a.cos()], [0.0, 2.0 - a.cos()], [a, a.cos()]];
        let bound = |points: &[[f64; 2]]| {
            planar_error(
                &reference,
                (-a, a),
                endpoints,
                &plane,
                (&axis, 1.0),
                points,
                1e-7,
            )
            .unwrap()
            .0
        };
        assert!(bound(&points) < 1e-8);
        points[1][1] += 0.001;
        // The added height is 0.002*t*(1-t), with integral 0.001/3. The
        // projected transverse derivative is at least cos(a)/sqrt(2).
        let exact_lower = 2.0 * a * a.cos() * h * 0.001 / 3.0;
        assert!(bound(&points) >= exact_lower);
        assert!(bound(&points) > 1e-6);
        points.reverse();
        assert!(
            planar_error(
                &reference,
                (-a, a),
                endpoints,
                &plane,
                (&axis, 1.0),
                &points,
                1e-7
            )
            .is_none()
        );
    }

    #[test]
    fn winding_unresolved_angular_order_and_wrong_finite_ends_refuse_the_bound() {
        let frame = Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        };
        let reference = Curve::Circle { frame, radius: 1.0 };
        let endpoints = [reference.value(0.0), reference.value(1.0)];
        for points in [
            [[0.0, 0.0], [2.0, 0.0], [1.0, 0.0]],
            [[0.0, 0.0], [4.0, 0.0], [8.0, 0.0]],
            [[0.0, 0.0], [0.5, 0.1], [1.0, 0.1]],
        ] {
            assert!(
                planar_error(
                    &reference,
                    (0.0, 1.0),
                    endpoints,
                    &frame,
                    (&frame, 1.0),
                    &points,
                    1e-7
                )
                .is_none()
            );
        }
    }
}
