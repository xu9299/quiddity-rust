//! Closed-form planar centroids under translations and both loop directions.
use super::*;

fn patch(offset: [f64; 3], rational: bool, reverse: bool) -> Part {
    use crate::brep::{Edge, Face, Loop};
    use crate::geom::{Curve, Frame};
    let local = if rational {
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 0.]]
    } else {
        [[0., 0., 0.], [1., 1., 0.], [0., 1., 0.]]
    };
    let points = local.map(|p| geom::add(offset, p));
    let poles = if rational {
        [[1., 0., 0.], [1., 1., 0.], [0., 1., 0.]]
    } else {
        [[0., 0., 0.], [0.5, 0., 0.], [1., 1., 0.]]
    };
    let curved = Curve::Nurbs(
        crate::nurbs::NurbsCurve::new(
            2,
            poles.map(|p| geom::add(offset, p)).to_vec(),
            if rational {
                vec![1., 0.5_f64.sqrt(), 1.]
            } else {
                vec![1.; 3]
            },
            vec![0., 0., 0., 1., 1., 1.],
        )
        .unwrap(),
    );
    let curves = [
        curved,
        Curve::Line {
            origin: points[1],
            dir: geom::sub(points[2], points[1]),
        },
        Curve::Line {
            origin: points[2],
            dir: geom::sub(points[0], points[2]),
        },
    ];
    let edges = curves
        .into_iter()
        .enumerate()
        .map(|(i, curve)| Edge {
            curve,
            start: points[i],
            end: points[(i + 1) % 3],
            vertices: (i, (i + 1) % 3),
            same_sense: true,
            samples: vec![points[i], points[(i + 1) % 3]],
        })
        .collect();
    Part::new(
        vec![Face {
            surface: Surface::Plane {
                frame: Frame {
                    origin: offset,
                    x: [1., 0., 0.],
                    y: [0., 1., 0.],
                    z: [0., 0., 1.],
                },
            },
            reversed: false,
            loops: vec![Loop {
                edges: if reverse {
                    (0..3).rev().map(|i| (i, false)).collect()
                } else {
                    (0..3).map(|i| (i, true)).collect()
                },
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
fn analytic_planar_centroids_do_not_amplify_translation_roundoff() {
    let mut failures = Vec::new();
    for offset in [[0.; 3], [123.456, -78.9, 41.3], [1e6, -2e6, 3e5]] {
        for rational in [false, true] {
            for reverse in [false, true] {
                let part = patch(offset, rational, reverse);
                let local = if rational {
                    [4. / (3. * std::f64::consts::PI); 2]
                } else {
                    [3. / 8., 3. / 5.]
                };
                let expected = geom::add(offset, [local[0], local[1], 0.]);
                let centre = part.face_centre(0).unwrap();
                let error = geom::dist(centre, expected);
                let scale = expected.into_iter().map(f64::abs).fold(1., f64::max);
                let tolerance = 32. * f64::EPSILON * scale;
                eprintln!(
                    "CENTROID_CONTROL rational={rational} reverse={reverse} offset={offset:?} centre={centre:?} error={error:e} tolerance={tolerance:e}"
                );
                if error > tolerance {
                    failures.push((offset, rational, reverse, error, tolerance));
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "analytic centroid deviations: {failures:?}"
    );
}
