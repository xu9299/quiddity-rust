use super::*;
use crate::geom;

fn bits((p, u, v): (V3, V3, V3)) -> [[u64; 3]; 3] {
    [p, u, v].map(|value| value.map(f64::to_bits))
}

fn surface(p: usize, q: usize, inner: &[f64]) -> NurbsSurface {
    let knots = |degree| {
        std::iter::repeat_n(0.0, degree + 1)
            .chain(inner.iter().copied())
            .chain(std::iter::repeat_n(1.0, degree + 1))
            .collect()
    };
    let count = |degree| degree + 1 + inner.len();
    let points = (0..count(p))
        .map(|i| {
            (0..count(q))
                .map(|j| {
                    [
                        i as f64 - 2.0,
                        j as f64 - 3.0,
                        ((i * j) as f64 * 0.03).sin(),
                    ]
                })
                .collect()
        })
        .collect();
    let weights = (0..count(p))
        .map(|i| {
            (0..count(q))
                .map(|j| 0.5 + ((i + 2 * j) % 7) as f64 * 0.125)
                .collect()
        })
        .collect();
    NurbsSurface::new(p, q, points, weights, knots(p), knots(q)).unwrap()
}

#[test]
fn span_polynomials_cover_degrees_repeated_knots_and_endpoints() {
    for degree in 0..SIZE {
        let knots: Vec<_> = std::iter::repeat_n(0.0, degree + 1)
            .chain([0.2, 0.4, 0.4, 0.7])
            .chain(std::iter::repeat_n(1.0, degree + 1))
            .collect();
        for span in degree..degree + 5 {
            if knots[span] == knots[span + 1] {
                continue;
            }
            let poly = Polynomial::new(&knots, degree, span).unwrap();
            for step in 0..=32 {
                let t = knots[span] + (knots[span + 1] - knots[span]) * step as f64 / 32.0;
                let (expected, derivative) =
                    super::super::basis_and_derivative(&knots, degree, span, t);
                let (actual, actual_derivative) = poly.evaluate(degree, t);
                for i in 0..=degree {
                    assert!((actual[i] - expected[i]).abs() < 2e-12);
                    assert!(
                        (actual_derivative[i] - derivative[i]).abs()
                            < 2e-11 * (1.0 + derivative[i].abs())
                    );
                }
                assert!((actual.iter().sum::<f64>() - 1.0).abs() < 2e-12);
                assert!(actual_derivative.iter().sum::<f64>().abs() < 2e-10);
            }
        }
    }
}

#[test]
fn contraction_preserves_both_partials_across_spans_and_fixed_parameters() {
    let params: Vec<_> = [
        -1.0,
        -0.0,
        0.2,
        0.4,
        0.40000000001,
        0.7,
        1.0,
        2.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ]
    .into_iter()
    .chain((0..21).map(|i| i as f64 / 20.0))
    .collect();
    let mut maximum = 0.0_f64;
    for (p, q) in [(0, 1), (1, 3), (3, 2), (5, 1), (5, 5), (7, 7)] {
        let surface = surface(p, q, &[0.2, 0.4, 0.4, 0.7]);
        for along_v in [false, true] {
            let mut evaluator = IntegralEvaluation::new(&surface, along_v);
            for &fixed in &params {
                for &varying in &params {
                    let (u, v) = if along_v {
                        (fixed, varying)
                    } else {
                        (varying, fixed)
                    };
                    let actual = evaluator.value_and_partials(u, v);
                    let expected = surface.value_and_partials(u, v);
                    for (a, b) in [actual.0, actual.1, actual.2]
                        .iter()
                        .flatten()
                        .zip([expected.0, expected.1, expected.2].iter().flatten())
                    {
                        let error = (a - b).abs() / (1.0 + b.abs());
                        maximum = maximum.max(error);
                        assert!(
                            error < 2e-11,
                            "degree {p}/{q}, direction {along_v}, ({u}, {v}): {a} vs {b}"
                        );
                    }
                }
            }
        }
    }
    println!("maximum normalized position/partial difference: {maximum:e}");
}

#[test]
fn unsupported_degrees_keep_original_bits() {
    for (p, q) in [(8, 5), (5, 8), (25, 25)] {
        let surface = surface(p, q, &[0.2, 0.4, 0.7]);
        for along_v in [false, true] {
            let mut evaluator = IntegralEvaluation::new(&surface, along_v);
            let mut original = SurfaceEvaluation::new(&surface);
            assert!(evaluator.prepared.is_none());
            for u in [-1.0, 0.0, 0.2, 0.61, 1.0, 2.0] {
                for v in [0.0, 0.27, 0.4, 1.0] {
                    assert_eq!(
                        bits(evaluator.value_and_partials(u, v)),
                        bits(original.value_and_partials(u, v))
                    );
                }
            }
        }
    }
}

#[test]
fn polynomial_budget_and_ill_resolved_intervals_use_basis_recurrence() {
    assert!(Polynomial::new(&[1.0, 1.0 + 2.0 * f64::EPSILON], 0, 0).is_none());
    assert!(Polynomial::new(&[0.0, 0.0], 0, 0).is_none());
    assert!(Polynomial::new(&[0.0, 1.0], SIZE, 0).is_none());
    let n = AXIS_BYTES / std::mem::size_of::<Option<Polynomial>>() + 1;
    let knots: Vec<_> = (0..=n).map(|i| i as f64 / n as f64).collect();
    let mut axis = Axis::new(&knots, 0, n);
    assert!(axis.polynomials.is_none());
    let actual = axis.prepare(0.31, &knots, 0, n).unwrap();
    assert_eq!(actual.basis[0], 1.0);
    assert_eq!(actual.derivative[0], 0.0);
    let knots = [1.0, 1.0, 1.0 + 2.0 * f64::EPSILON, 1.0 + 2.0 * f64::EPSILON];
    let mut axis = Axis::new(&knots, 1, 2);
    assert!(axis.polynomials.as_ref().unwrap()[1].is_none());
    let actual = axis.prepare(1.0, &knots, 1, 2).unwrap();
    let (b, d) = super::super::basis_and_derivative(&knots, 1, 1, 1.0);
    assert_eq!(&*actual.basis, &*b);
    assert_eq!(&*actual.derivative, &*d);
}

#[test]
fn nonfinite_values_and_weighted_overflow_keep_original_behavior() {
    for (coordinate, weight) in [(f64::INFINITY, 1.0), (f64::NAN, 1.0), (1e300, 1e300)] {
        let surface = NurbsSurface::new(
            1,
            1,
            vec![vec![[coordinate, 1.0, 0.0]; 2]; 2],
            vec![vec![weight; 2]; 2],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
        )
        .unwrap();
        for along_v in [false, true] {
            let mut evaluator = IntegralEvaluation::new(&surface, along_v);
            let mut original = SurfaceEvaluation::new(&surface);
            for (u, v) in [
                (0.0, 0.0),
                (0.3, 0.8),
                (1.0, 1.0),
                (f64::NAN, 0.5),
                (f64::INFINITY, 0.5),
            ] {
                assert_eq!(
                    bits(evaluator.value_and_partials(u, v)),
                    bits(original.value_and_partials(u, v))
                );
            }
        }
    }
    let surface = surface(5, 5, &[0.4]);
    let mut evaluator = IntegralEvaluation::new(&surface, false);
    let mut original = SurfaceEvaluation::new(&surface);
    assert_eq!(
        bits(evaluator.value_and_partials(f64::NAN, 0.5)),
        bits(original.value_and_partials(f64::NAN, 0.5))
    );
}

#[test]
fn analytic_cylinder_area_flux_and_area_vector() {
    let points = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        .map(|[x, y]| vec![[x, y, 0.0], [x, y, 2.0]])
        .to_vec();
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let surface = NurbsSurface::new(
        2,
        1,
        points,
        vec![vec![1.0; 2], vec![w; 2], vec![1.0; 2]],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
    )
    .unwrap();
    // Independent composite Simpson integration, with analytic five-term expectations.
    let panels = 4096;
    for along_v in [false, true] {
        let mut evaluator = IntegralEvaluation::new(&surface, along_v);
        let mut total = [0.0; 5];
        for i in 0..=panels {
            let (p, du, dv) = evaluator.value_and_partials(i as f64 / panels as f64, 0.37);
            let n = geom::cross(du, dv);
            let y = [geom::norm(n), geom::dot(p, n), n[0], n[1], n[2]];
            let weight = if i == 0 || i == panels {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for k in 0..5 {
                total[k] += weight * y[k] / (3.0 * panels as f64);
            }
        }
        for (a, b) in
            total
                .into_iter()
                .zip([std::f64::consts::PI, std::f64::consts::PI, 2.0, 2.0, 0.0])
        {
            assert!((a - b).abs() < 2e-12, "direction {along_v}: {a} vs {b}");
        }
    }
}

#[test]
fn finite_weight_contrasts_near_span_ends_keep_the_analytic_sweep() {
    // A flat sweep with C(u) = (u + (w - 1) u^p) / (1 + (w - 1) u^p).
    // Its last control row alone has weight w. A near-zero basis rounding error
    // can leave every intermediate finite while moving C outside its [0, 1] hull.
    for degree in [3, 5, 7] {
        for weight in [1.0, 16.0, 256.0, 512.0, 1e8, 1e16] {
            let surface = NurbsSurface::new(
                degree,
                1,
                (0..=degree)
                    .map(|i| {
                        vec![
                            [i as f64 / degree as f64, 0.0, 0.0],
                            [i as f64 / degree as f64, 1.0, 0.0],
                        ]
                    })
                    .collect(),
                (0..=degree)
                    .map(|i| vec![if i == degree { weight } else { 1.0 }; 2])
                    .collect(),
                std::iter::repeat_n(0.0, degree + 1)
                    .chain(std::iter::repeat_n(1.0, degree + 1))
                    .collect(),
                vec![0.0, 0.0, 1.0, 1.0],
            )
            .unwrap();
            for along_v in [false, true] {
                let mut evaluator = IntegralEvaluation::new(&surface, along_v);
                for u in [0.0_f64, 1e-7, 0.001, 0.005, 0.02, 0.5, 1.0] {
                    let z = (weight - 1.0) * u.powi(degree as i32);
                    let dz = (weight - 1.0) * degree as f64 * u.powi(degree as i32 - 1);
                    let denominator = 1.0 + z;
                    let expected = [
                        [(u + z) / denominator, 0.37, 0.0],
                        [
                            (denominator + dz * (1.0 - u)) / denominator.powi(2),
                            0.0,
                            0.0,
                        ],
                        [0.0, 1.0, 0.0],
                    ];
                    let (p, du, dv) = evaluator.value_and_partials(u, 0.37);
                    for (actual, expected) in
                        [p, du, dv].iter().flatten().zip(expected.iter().flatten())
                    {
                        assert!(
                            (actual - expected).abs() < 2e-11 * (1.0 + expected.abs()),
                            "degree {degree}, weight {weight}, u {u}, direction {along_v}: {actual} vs {expected}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
#[ignore = "serial timing experiment; no performance assertion"]
fn benchmark_synthetic_inner_sweeps() {
    use std::hint::black_box;
    use std::time::Instant;
    // Analytic control grids generated here; no external geometry is required.
    let spans = 128;
    let repeats = 32;
    let evaluations = repeats * 128 * spans * 16;
    let inner: Vec<_> = (1..spans)
        .flat_map(|i| [i as f64 / spans as f64; 3])
        .collect();
    for q in [1, 5] {
        let mut surface = surface(5, q, &inner);
        // Use one span along the fixed axis, as in swept or narrow spline patches.
        surface.knots_v = std::iter::repeat_n(0.0, q + 1)
            .chain(std::iter::repeat_n(1.0, q + 1))
            .collect();
        for row in &mut surface.control_points {
            row.truncate(q + 1);
        }
        for row in &mut surface.weights {
            row.truncate(q + 1);
        }
        for optimized in [false, true, true, false] {
            let start = Instant::now();
            let mut original = (!optimized).then(|| SurfaceEvaluation::new(&surface));
            let mut candidate = optimized.then(|| IntegralEvaluation::new(&surface, false));
            let mut checksum = [0.0; 9];
            for _ in 0..repeats {
                for fixed in 0..128 {
                    let v = (fixed as f64 + 0.5) / 128.0;
                    for span in 0..spans {
                        for node in 0..16 {
                            let u = (span as f64 + (node as f64 + 0.5) / 16.0) / spans as f64;
                            let (p, du, dv) = if optimized {
                                candidate.as_mut().unwrap().value_and_partials(u, v)
                            } else {
                                original.as_mut().unwrap().value_and_partials(u, v)
                            };
                            for (sum, value) in
                                checksum.iter_mut().zip([p, du, dv].into_iter().flatten())
                            {
                                *sum += value;
                            }
                        }
                    }
                }
            }
            black_box(checksum);
            println!(
                "{{\"degree_v\":{q},\"optimized\":{optimized},\"evaluations\":{evaluations},\"seconds\":{},\"checksum\":{checksum:?}}}",
                start.elapsed().as_secs_f64()
            );
        }
    }
}
