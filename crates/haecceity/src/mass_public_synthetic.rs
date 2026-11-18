//! Public analytic triangular parabolic patch with deliberately dense spline knots.
//! These fixtures are constructed from polynomial identities, with no external geometry.
use super::*;
use std::cell::Cell;

fn patch(spans: usize) -> Part {
    let c = 0.5;
    let offset = [1000., 2000., 3000.];
    let mut knots = vec![0., 0., 0.];
    knots.extend((1..spans).map(|i| i as f64 / spans as f64));
    knots.extend([1., 1., 1.]);
    let controls: Vec<_> = (0..spans + 2)
        .map(|i| {
            let a = knots[i + 1];
            let b = knots[i + 2];
            // Quadratic B-spline coefficients of u and u^2 are (a+b)/2 and a*b.
            [(a + b) / 2., c * a * b]
        })
        .collect();
    let surface = crate::nurbs::NurbsSurface::new(
        2,
        1,
        controls
            .iter()
            .map(|p| {
                [0., 1.]
                    .map(|v| [offset[0] + p[0], offset[1] + v, offset[2] + p[1]])
                    .to_vec()
            })
            .collect(),
        vec![vec![1.; 2]; spans + 2],
        knots.clone(),
        vec![0., 0., 1., 1.],
    )
    .unwrap();
    let points = [
        surface.value(0., 0.),
        surface.value(1., 1.),
        surface.value(0., 1.),
    ];
    let curves = [
        Curve::Nurbs(
            crate::nurbs::NurbsCurve::new(
                2,
                controls
                    .iter()
                    .map(|p| [offset[0] + p[0], offset[1] + p[0], offset[2] + p[1]])
                    .collect(),
                vec![1.; spans + 2],
                knots.clone(),
            )
            .unwrap(),
        ),
        Curve::Nurbs(
            crate::nurbs::NurbsCurve::new(
                2,
                controls
                    .iter()
                    .rev()
                    .map(|p| [offset[0] + p[0], offset[1] + 1., offset[2] + p[1]])
                    .collect(),
                vec![1.; spans + 2],
                knots,
            )
            .unwrap(),
        ),
        Curve::Line {
            origin: points[2],
            dir: geom::sub(points[0], points[2]),
        },
    ];
    let edges = curves
        .into_iter()
        .enumerate()
        .map(|(i, curve)| crate::brep::Edge {
            curve,
            start: points[i],
            end: points[(i + 1) % 3],
            vertices: (i, (i + 1) % 3),
            same_sense: true,
            samples: vec![points[i], points[(i + 1) % 3]],
        })
        .collect();
    Part::new(
        vec![crate::brep::Face {
            surface: Surface::Freeform {
                kind: "analytic synthetic",
                surface: Box::new(surface),
            },
            reversed: false,
            loops: vec![crate::brep::Loop {
                edges: vec![(0, true), (1, true), (2, true)],
                vertex: None,
            }],
            solid: None,
            pcurves: vec![],
        }],
        edges,
        vec![],
    )
}

#[test]
#[ignore = "manual analytic scaling study"]
fn analytic_dense_patch_pairs() {
    let output = std::env::var_os("QUIDDITY_QUADRATURE_REPORT").map(std::path::PathBuf::from);
    let c = 0.5_f64;
    // S(u,v)=(u,v,c*u^2), 0<=u<=v<=1, translated by offset.
    let full = 0.5 * (1. + 4. * c * c).sqrt() + (2. * c).asinh() / (4. * c);
    let weighted = ((1. + 4. * c * c).powf(1.5) - 1.) / (12. * c * c);
    let exact = [
        full - weighted,
        -c / 12. - 1000. * c / 3. + 1500.,
        -c / 3.,
        0.,
        0.5,
    ];
    let mut integrals = [0.; 6];
    integrals[0] = full;
    integrals[1] = weighted;
    for n in 2..6 {
        integrals[n] = (2. * 2_f64.sqrt() - (n - 1) as f64 * integrals[n - 2]) / (n + 2) as f64;
    }
    let [i0, i1, i2, i3, i4, i5] = integrals;
    let local_mean = [
        (i1 - i2) / exact[0],
        (i0 - i2) / (2. * exact[0]),
        c * (i2 - i3) / exact[0],
    ];
    let raw = [
        [i2 - i3, (i1 - i3) / 2., c * (i3 - i4)],
        [(i1 - i3) / 2., (i0 - i3) / 3., c * (i2 - i4) / 2.],
        [c * (i3 - i4), c * (i2 - i4) / 2., c * c * (i4 - i5)],
    ];
    let central: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| raw[i][j] - exact[0] * local_mean[i] * local_mean[j])
    });
    let mut rows = Vec::new();
    for spans in [32, 128, 512] {
        let mut part = patch(spans);
        for u in [0., 0.125, 0.37, 0.625, 1.] {
            for v in [0., 0.2, 1.] {
                let (p, su, sv) = point_and_partials(&part.faces[0].surface, u, v);
                assert!(geom::dist(p, [1000. + u, 2000. + v, 3000. + c * u * u]) < 1e-10);
                assert!(geom::dist(su, [1., 0., 2. * c * u]) < 1e-10);
                assert!(geom::dist(sv, [0., 1., 0.]) < 1e-10);
            }
        }
        assert!(part.uv_loops(0).is_some());
        assert!(part.domain(0).is_some());
        assert!(part.loop_placement(0).is_some());
        assert!(part.short_axis_walk_is_closed(0));
        for canonical in if spans == 128 {
            [false, true]
        } else {
            [true, false]
        } {
            let calls = Cell::new(0_u64);
            let density = |p, su, sv| {
                calls.set(calls.get() + 1);
                mass_density(p, su, sv)
            };
            let start = std::time::Instant::now();
            let values = part
                .integrate_face_axis(
                    0,
                    &density,
                    IntegralRule {
                        components: Components::Prefix(5),
                        observe_mass: None,
                    },
                    &mut None,
                    !canonical,
                )
                .unwrap();
            let seconds = start.elapsed().as_secs_f64();
            let errors: [f64; 5] =
                std::array::from_fn(|k| (values[k] - exact[k]).abs() / exact[k].abs().max(1.));
            rows.push(serde_json::json!({"kind":"mass","spans":spans,"canonical_axis":canonical,"seconds":seconds,
                "density_calls":calls.get(),"values":values,"exact":exact,"normalized_errors":errors}));
            if let Some(path) = &output {
                std::fs::write(path, serde_json::to_vec_pretty(&rows).unwrap()).unwrap();
            }
            assert!(
                errors.into_iter().all(|e| e < 1e-9),
                "spans={spans}, canonical={canonical}, values={values:?}, exact={exact:?}"
            );
            assert!(part.cache[0].mass.get().is_none());
            if spans == 32 && !canonical {
                // Exercise the actual public moment path with automatic axis selection.
                let start = std::time::Instant::now();
                let moments = part.compute_face_moments(0).unwrap();
                let seconds = start.elapsed().as_secs_f64();
                let centroid: [f64; 3] =
                    std::array::from_fn(|i| [1000., 2000., 3000.][i] + local_mean[i]);
                let centroid_error = geom::dist(moments.centroid, centroid);
                let second_error = (0..3)
                    .flat_map(|i| {
                        (0..3).map(move |j| (moments.second_moments[i][j] - central[i][j]).abs())
                    })
                    .fold(0., f64::max);
                let vector_error = geom::dist(moments.area_vector, [-c / 3., 0., 0.5]);
                rows.push(serde_json::json!({"kind":"moments","spans":spans,"canonical_axis":canonical,"seconds":seconds,
                    "area":moments.area,"centroid":moments.centroid,
                    "second_moments":moments.second_moments,"expected_centroid":centroid,"expected_second_moments":central,
                    "centroid_error":centroid_error,"second_error":second_error,"vector_error":vector_error}));
                if let Some(path) = &output {
                    std::fs::write(path, serde_json::to_vec_pretty(&rows).unwrap()).unwrap();
                }
                assert!(
                    (moments.area - exact[0]).abs() < 1e-10
                        && centroid_error < 1e-9
                        && second_error < 1e-10
                        && vector_error < 1e-10
                );
                part.cache[0].mass.take();
                part.cache[0].moments.take();
            }
        }
    }
    println!("{}", serde_json::to_string(&rows).unwrap());
}
