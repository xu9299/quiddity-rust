//! A sufficient floating-point finite-value certificate for native analytic trims.
//! Failure always retains the sampled check; foreign residuals and finite NURBS
//! carrier domains are outside its scope. Each positive source knot span is
//! covered by two closed intervals. The last active span must be positive because
//! the evaluator selects it at the upper endpoint even after repeated knots.
use super::{ParamCurve, parameters::Axis};
use crate::geom::{Curve, Surface};
use crate::nurbs::NurbsCurve;

/// A native analytic carrier has no finite chart boundary. A successful
/// certificate bounds every source evaluation and every lift, over the whole
/// source domain, rather than assuming that positive weights imply a finite
/// floating-point quotient. Failure leaves the original sampled checks intact.
pub(super) fn native(curve: &ParamCurve, surface: &Surface, interval: (f64, f64)) -> bool {
    if curve.transfer.is_some()
        || curve.axes != [Axis::ID; 2]
        || ![interval.0, interval.1]
            .iter()
            .all(|x| x.is_finite() && x.abs() < 1e50)
        || !analytic_lift_finite(surface, 0.)
    {
        return false;
    }
    let Curve::Nurbs(source) = &*curve.source else {
        return false;
    };
    finite_value_bound(source).is_some_and(|bound| analytic_lift_finite(surface, bound))
}

fn finite_input_guard(curve: &NurbsCurve) -> bool {
    const LIMIT: f64 = 1e50;
    const FLOOR: f64 = 1e-50;
    let bounded = |x: f64| x.is_finite() && x.abs() < LIMIT;
    let knots = curve.knots();
    let n = curve.control_points().len();
    curve.degree() <= 3
        && curve.control_points().iter().flatten().all(|&x| bounded(x))
        && curve.weights().iter().all(|&w| w > FLOOR && w < LIMIT)
        && knots.iter().all(|&x| bounded(x))
        && knots.windows(2).all(|p| p[0] == p[1] || p[1] - p[0] >= FLOOR)
        // find_span uses n - 1 at the upper endpoint. A zero last span can
        // make all basis values zero even in otherwise accepted knot arrays.
        && knots[n - 1] < knots[n]
}

#[derive(Clone, Copy, Debug)]
struct Nonnegative {
    lo: f64,
    hi: f64,
}

impl Nonnegative {
    const ZERO: Self = Self { lo: 0., hi: 0. };
    fn exact(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn rounded(lo: f64, hi: f64) -> Option<Self> {
        let value = Self {
            lo: lo.next_down().max(0.),
            hi: hi.next_up(),
        };
        (lo >= 0. && lo <= hi && value.hi.is_finite()).then_some(value)
    }
    fn add(self, other: Self) -> Option<Self> {
        Self::rounded(self.lo + other.lo, self.hi + other.hi)
    }
    fn mul(self, other: Self) -> Option<Self> {
        Self::rounded(self.lo * other.lo, self.hi * other.hi)
    }
    fn basis_division(self, denominator: Self) -> Option<Self> {
        // Match the actual basis evaluator's small-denominator branch.
        if denominator.hi < 1e-300 {
            return Some(Self::ZERO);
        }
        let lo = if denominator.lo < 1e-300 {
            0.
        } else {
            self.lo / denominator.hi
        };
        let hi = self.hi / denominator.lo.max(1e-300);
        Self::rounded(lo, hi)
    }
}

// Propagate enclosures through the actual basis recurrence and homogeneous
// sums. No partition-of-unity or floating-point normalization is assumed.
fn finite_value_bound(curve: &NurbsCurve) -> Option<f64> {
    if !finite_input_guard(curve) {
        return None;
    }
    let p = curve.degree();
    let n = curve.control_points().len();
    let knots = curve.knots();
    let mut result = 0_f64;
    for span in p..n {
        let (a, b) = (knots[span], knots[span + 1]);
        if a == b {
            continue;
        }
        let middle = a + (b - a) * 0.5;
        for (a, b) in [(a, middle), (middle, b)] {
            let mut basis = [Nonnegative::ZERO; 4];
            let mut left = [Nonnegative::ZERO; 4];
            let mut right = [Nonnegative::ZERO; 4];
            basis[0] = Nonnegative::exact(1.);
            for j in 1..=p {
                left[j] = Nonnegative::rounded(a - knots[span + 1 - j], b - knots[span + 1 - j])?;
                right[j] = Nonnegative::rounded(knots[span + j] - b, knots[span + j] - a)?;
                let mut saved = Nonnegative::ZERO;
                for r in 0..j {
                    let denominator = right[r + 1].add(left[j - r])?;
                    let temp = basis[r].basis_division(denominator)?;
                    basis[r] = saved.add(right[r + 1].mul(temp)?)?;
                    saved = left[j - r].mul(temp)?;
                }
                basis[j] = saved;
            }
            let mut denominator = Nonnegative::ZERO;
            let mut numerator = [0_f64; 3];
            for (j, basis) in basis.iter().take(p + 1).enumerate() {
                let i = span - p + j;
                let weight = Nonnegative::exact(curve.weights()[i]).mul(*basis)?;
                denominator = denominator.add(weight)?;
                for (axis, bound) in numerator.iter_mut().enumerate() {
                    let term = (curve.control_points()[i][axis].abs() * weight.hi).next_up();
                    *bound = (*bound + term).next_up();
                    if !bound.is_finite() {
                        return None;
                    }
                }
            }
            if denominator.lo <= 0. {
                return None;
            }
            for numerator in numerator {
                let bound = (numerator / denominator.lo).next_up();
                if !bound.is_finite() {
                    return None;
                }
                result = result.max(bound);
            }
        }
    }
    Some(result)
}

// The five analytic carriers have no finite parameter domain. For source
// coordinates below 1e50, each local value is below 3e50: sin/cos are bounded,
// and only the cone adds a radius to a coordinate. Frame coefficients below
// 1e50 keep each product below 3e100 and the four-term world sum below 1e102.
// These loose bounds leave far more room than floating-point rounding needs.
// This certifies finite evaluation only; it says nothing about another face.
fn analytic_lift_finite(surface: &crate::geom::Surface, source_bound: f64) -> bool {
    use crate::geom::Surface;
    let bounded = |x: f64| x.is_finite() && x.abs() < 1e50;
    if !bounded(source_bound) {
        return false;
    }
    let (frame, dimensions) = match surface {
        Surface::Plane { frame } => (frame, true),
        Surface::Cylinder { frame, radius } | Surface::Sphere { frame, radius } => {
            (frame, bounded(*radius))
        }
        Surface::Cone {
            frame,
            radius,
            semi_angle,
        } => (frame, bounded(*radius) && bounded(*semi_angle)),
        Surface::Torus {
            frame,
            major,
            minor,
        } => (frame, bounded(*major) && bounded(*minor)),
        _ => return false,
    };
    dimensions
        && [frame.origin, frame.x, frame.y, frame.z]
            .iter()
            .flatten()
            .all(|&x| bounded(x))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    #[test]
    fn bounded_low_degree_curves_remain_finite_at_knots_and_between_them() {
        let mut accepted = 0;
        let mut deferred = 0;
        let mut evaluations = 0;
        for degree in 0..=3 {
            for n in degree + 1..=degree + 7 {
                for shape in 0..32_u64 {
                    for knot_scale in [1e-49, 1., 1e47] {
                        for pattern in 0..3 {
                            let mut seed = shape + 1;
                            let mut value = -10.;
                            let knots: Vec<_> = (0..n + degree + 1)
                                .map(|_| {
                                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                                    if seed >> 62 != 0 {
                                        value += 1.;
                                    }
                                    value * knot_scale
                                })
                                .collect();
                            let points: Vec<_> = (0..n)
                                .map(|i| {
                                    let x = (i as f64 * 0.5).sin();
                                    [x * 1e49, -x * 1e-49, i as f64]
                                })
                                .collect();
                            let weights: Vec<_> = (0..n)
                                .map(|i| match pattern {
                                    0 => 1.,
                                    1 if i % 2 == 0 => 1e-49,
                                    2 if i % 2 == 0 => 1e49,
                                    _ => 0.5,
                                })
                                .collect();
                            let Some(curve) = NurbsCurve::new(degree, points, weights, knots)
                            else {
                                continue;
                            };
                            if !finite_value_bound(&curve).is_some() {
                                deferred += 1;
                                continue;
                            }
                            accepted += 1;
                            let bound = finite_value_bound(&curve).unwrap();
                            let (lo, hi) = curve.domain();
                            let mut probes: Vec<_> =
                                (0..=64).map(|i| lo + (hi - lo) * i as f64 / 64.).collect();
                            for &knot in curve.knots() {
                                probes.extend([knot.next_down(), knot, knot.next_up()]);
                            }
                            for t in probes {
                                assert!(
                                    curve.value(t).iter().all(|x| x.is_finite()),
                                    "degree {degree}, n {n}, shape {shape}, scale {knot_scale}, pattern {pattern}, t {t}"
                                );
                                assert!(curve.value(t).iter().all(|x| x.abs() <= bound));
                                evaluations += 1;
                            }
                        }
                    }
                }
            }
        }
        assert!(accepted > 3000 && deferred > 1000);
        println!("Accepted {accepted} curves, deferred {deferred}, checked {evaluations} values");
    }

    #[test]
    fn malformed_end_spans_and_unbounded_operands_are_deferred() {
        let end_zero = NurbsCurve::new(
            1,
            vec![[1., 2., 3.]; 3],
            vec![1.; 3],
            vec![0., 0., 1., 1., 2.],
        )
        .unwrap();
        assert!(end_zero.value(1.).iter().any(|x| !x.is_finite()));
        assert!(!finite_value_bound(&end_zero).is_some());
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e50, -1e50] {
            let bad = NurbsCurve::new(
                2,
                vec![[value, 0., 0.]; 3],
                vec![1.; 3],
                vec![0., 0., 0., 1., 1., 1.],
            )
            .unwrap();
            assert!(!finite_value_bound(&bad).is_some());
        }
        for weight in [f64::NAN, f64::INFINITY, 1e50, 1e-50, f64::MIN_POSITIVE] {
            let bad = NurbsCurve::new(
                2,
                vec![[1., 0., 0.]; 3],
                vec![weight; 3],
                vec![0., 0., 0., 1., 1., 1.],
            )
            .unwrap();
            assert!(!finite_value_bound(&bad).is_some());
        }
        let tiny = NurbsCurve::new(
            2,
            vec![[1., 0., 0.]; 3],
            vec![1.; 3],
            vec![0., 0., 0., 1e-51, 1e-51, 1e-51],
        )
        .unwrap();
        assert!(!finite_value_bound(&tiny).is_some());
        let high = NurbsCurve::new(
            4,
            vec![[1., 0., 0.]; 5],
            vec![1.; 5],
            vec![0., 0., 0., 0., 0., 1., 1., 1., 1., 1.],
        )
        .unwrap();
        assert!(!finite_value_bound(&high).is_some());
    }

    #[test]
    fn native_certificate_is_limited_to_bounded_analytic_lifts() {
        use crate::geom::Frame;
        let source = Arc::new(Curve::Nurbs(
            NurbsCurve::new(
                2,
                vec![[-1e40, 1e40, 0.], [1e40, -1e40, 0.], [1e40, 1e40, 0.]],
                vec![1., 0.75, 1.],
                vec![0., 0., 0., 1., 1., 1.],
            )
            .unwrap(),
        ));
        let curve = ParamCurve {
            source: source.clone(),
            axes: [Axis::ID; 2],
            transfer: None,
            inversion_pool: None,
        };
        let mut evaluations = 0;
        for size in [1e-49, 1., 1e20, 1e48] {
            let frame = Frame {
                origin: [2. * size, -size, 3. * size],
                x: [size, 0., 0.],
                y: [0., -size, 0.],
                z: [0., 0., size],
            };
            let surfaces = [
                Surface::Plane { frame },
                Surface::Cylinder {
                    frame,
                    radius: size,
                },
                Surface::Cone {
                    frame,
                    radius: size,
                    semi_angle: 0.45,
                },
                Surface::Sphere {
                    frame,
                    radius: size,
                },
                Surface::Torus {
                    frame,
                    major: 2. * size,
                    minor: size,
                },
            ];
            for surface in &surfaces {
                assert!(native(&curve, surface, (0., 1.)));
                for i in 0..=1024 {
                    let uv = source.value(i as f64 / 1024.);
                    assert!(
                        surface
                            .value(uv[0] * 1. + 0., uv[1] * 1. + 0.)
                            .iter()
                            .all(|x| x.is_finite())
                    );
                    evaluations += 1;
                }
                for interval in [(f64::NAN, 1.), (0., f64::INFINITY), (0., 1e50)] {
                    assert!(!native(&curve, surface, interval));
                }
            }
        }
        assert_eq!(evaluations, 20500);
    }

    #[test]
    fn foreign_curves_axes_and_nonfinite_carriers_keep_the_original_checks() {
        use crate::geom::Frame;
        let frame = Frame {
            origin: [0.; 3],
            x: [1., 0., 0.],
            y: [0., 1., 0.],
            z: [0., 0., 1.],
        };
        let plane = Surface::Plane { frame };
        let source = Arc::new(Curve::Nurbs(
            NurbsCurve::new(1, vec![[0.; 3], [1.; 3]], vec![1.; 2], vec![0., 0., 1., 1.]).unwrap(),
        ));
        let mut curve = ParamCurve {
            source,
            axes: [Axis::ID; 2],
            transfer: None,
            inversion_pool: None,
        };
        assert!(native(&curve, &plane, (0., 1.)));
        for axis in [
            Axis::Linear {
                scale: 2.,
                offset: 0.,
            },
            Axis::Angle { radians: 1. },
        ] {
            curve.axes[0] = axis;
            assert!(!native(&curve, &plane, (0., 1.)));
        }
        curve.axes = [Axis::ID; 2];
        for value in [1e50, 1e200, f64::NAN, f64::INFINITY] {
            let bad = Surface::Plane {
                frame: Frame {
                    origin: [value; 3],
                    ..frame
                },
            };
            assert!(!native(&curve, &bad, (0., 1.)));
            assert!(!native(
                &curve,
                &Surface::Cylinder {
                    frame,
                    radius: value
                },
                (0., 1.)
            ));
            assert!(!native(
                &curve,
                &Surface::Cone {
                    frame,
                    radius: 1.,
                    semi_angle: value
                },
                (0., 1.)
            ));
            assert!(!native(
                &curve,
                &Surface::Torus {
                    frame,
                    major: 1.,
                    minor: value
                },
                (0., 1.)
            ));
        }
        assert!(!native(
            &curve,
            &Surface::Other { kind: "unresolved" },
            (0., 1.)
        ));
        curve.transfer = Some(Box::new(super::super::Transfer {
            carrier: Arc::new(plane.clone()),
            face: Arc::new(plane.clone()),
            carrier_spans: (f64::INFINITY, f64::INFINITY),
            face_spans: (f64::INFINITY, f64::INFINITY),
            face_collapsed: (None, None),
        }));
        assert!(!native(&curve, &plane, (0., 1.)));
    }
}
