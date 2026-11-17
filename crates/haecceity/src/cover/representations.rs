//! Independent STEP trim representations mapped to the held analytic surface.
//!
//! A file uncertainty is only an eligibility check for coherent endpoints. It never changes
//! the coverage budget. An unevaluable pcurve cannot be silently omitted from that check.

use super::FaceRef;
use crate::brep::{Pcurve, RawCurve2d};
use crate::geom::{self, Curve, Frame, Surface};

/// The queried face's pcurves, and an explicitly preferred pcurve on an adjacent face.
/// The 3D curve is held separately by the edge. A preferred representation can have another
/// parameter speed, but STEP requires its sense to agree with the 3D representation.
pub(super) fn alternatives((part, face): FaceRef<'_>, edge: usize) -> Option<Vec<Curve>> {
    let mut out = Vec::new();
    let units = part.file_units().unwrap_or((1.0, 1.0));
    for &(index, ref curve) in part.file_pcurves(face) {
        if index == edge {
            out.push(mapped(curve, &part.faces[face].surface, units)?);
        }
    }
    if part.file_pcurves(face).is_empty() {
        for (index, curve) in &part.faces[face].pcurves {
            if *index == edge {
                let Pcurve::Line { point, dir } = curve else {
                    return None;
                };
                out.push(mapped(
                    &RawCurve2d::Line {
                        point: [point.0, point.1],
                        direction: [dir.0, dir.1],
                        magnitude: 1.0,
                    },
                    &part.faces[face].surface,
                    (1.0, 1.0),
                )?);
            }
        }
    }
    for &adjacent in &part.edge_faces()[edge] {
        if adjacent == face {
            continue;
        }
        for (index, curve) in part.primary_pcurves(adjacent) {
            if *index == edge {
                out.push(mapped(curve, &part.faces[adjacent].surface, units)?);
            }
        }
    }
    Some(out)
}

fn mapped(raw: &RawCurve2d, surface: &Surface, (mm, rad): (f64, f64)) -> Option<Curve> {
    match surface {
        Surface::Plane { frame } => {
            let point = |p: [f64; 2]| frame.to_world([p[0] * mm, p[1] * mm, 0.0]);
            let conic = |p: [f64; 2], reference: Option<[f64; 2]>, major, minor| {
                let d = reference.unwrap_or([1.0, 0.0]);
                let x = geom::unit([d[0], d[1], 0.0])?;
                let x = frame.dir_to_world(x);
                Some(Curve::Ellipse {
                    frame: Frame {
                        origin: point(p),
                        x,
                        y: geom::cross(frame.z, x),
                        z: frame.z,
                    },
                    major: major * mm,
                    minor: minor * mm,
                })
            };
            match raw {
                RawCurve2d::Unresolved => None,
                RawCurve2d::Line {
                    point: p,
                    direction: d,
                    magnitude,
                } => Some(Curve::Line {
                    origin: point(*p),
                    dir: geom::unit(frame.dir_to_world([
                        d[0] * magnitude * mm,
                        d[1] * magnitude * mm,
                        0.0,
                    ]))?,
                }),
                RawCurve2d::Circle {
                    point,
                    reference,
                    radius,
                } => conic(*point, *reference, *radius, *radius),
                RawCurve2d::Ellipse {
                    point,
                    reference,
                    semi_axis_1,
                    semi_axis_2,
                } => conic(*point, *reference, *semi_axis_1, *semi_axis_2),
                RawCurve2d::BSpline {
                    degree,
                    points,
                    multiplicities,
                    knots,
                    weights,
                } => {
                    let degree = usize::try_from(*degree).ok()?;
                    let multiplicities = multiplicities
                        .iter()
                        .map(|&n| usize::try_from(n).ok())
                        .collect::<Option<Vec<_>>>()?;
                    let count = multiplicities
                        .iter()
                        .try_fold(0usize, |a, &b| a.checked_add(b))?;
                    if degree > crate::nurbs::MAX_DEGREE
                        || points.len() > 4096
                        || knots.len() != multiplicities.len()
                        || count != points.len().checked_add(degree + 1)?
                        || weights.as_ref().is_some_and(|w| {
                            w.len() != points.len() || w.iter().any(|x| !x.is_finite() || *x <= 0.0)
                        })
                    {
                        return None;
                    }
                    crate::nurbs::NurbsCurve::new(
                        degree,
                        points.iter().copied().map(point).collect(),
                        weights.clone().unwrap_or_else(|| vec![1.0; points.len()]),
                        knots
                            .iter()
                            .zip(multiplicities)
                            .flat_map(|(&k, n)| std::iter::repeat_n(k, n))
                            .collect(),
                    )
                    .map(Curve::Nurbs)
                }
            }
        }
        Surface::Cylinder { frame, radius } => {
            let RawCurve2d::Line {
                point,
                direction,
                magnitude,
            } = raw
            else {
                return None;
            };
            let (u, v) = (point[0] * rad, point[1] * mm);
            let (du, dv) = (
                direction[0] * magnitude * rad,
                direction[1] * magnitude * mm,
            );
            if ![u, v, du, dv, *radius].into_iter().all(f64::is_finite) {
                return None;
            }
            if du == 0.0 && dv != 0.0 {
                Some(Curve::Line {
                    origin: surface.value(u, v),
                    dir: geom::scale(frame.z, dv.signum()),
                })
            } else if dv == 0.0 && du != 0.0 {
                let x = frame.dir_to_world([u.cos(), u.sin(), 0.0]);
                let z = geom::scale(frame.z, du.signum());
                Some(Curve::Circle {
                    frame: Frame {
                        origin: frame.to_world([0.0, 0.0, v]),
                        x,
                        y: geom::cross(z, x),
                        z,
                    },
                    radius: *radius,
                })
            } else {
                None
            }
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame() -> Frame {
        Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        }
    }

    #[test]
    fn planar_conics_keep_their_file_units_and_direction() {
        let curve = mapped(
            &RawCurve2d::Ellipse {
                point: [1.0, 2.0],
                reference: Some([0.0, 3.0]),
                semi_axis_1: 4.0,
                semi_axis_2: 2.0,
            },
            &Surface::Plane { frame: frame() },
            (1000.0, 1.0),
        )
        .unwrap();
        assert!(geom::dist(curve.value(0.0), [1000.0, 6000.0, 0.0]) < 1e-10);
        assert!(geom::dist(curve.derivative(0.0), [-2000.0, 0.0, 0.0]) < 1e-10);
    }

    #[test]
    fn cylindrical_parameter_lines_are_mapped_without_a_chord() {
        let surface = Surface::Cylinder {
            frame: frame(),
            radius: 2.0,
        };
        let arc = mapped(
            &RawCurve2d::Line {
                point: [0.0, 3.0],
                direction: [-1.0, 0.0],
                magnitude: 1.0,
            },
            &surface,
            (1.0, 1.0),
        )
        .unwrap();
        assert!(geom::dist(arc.value(0.0), [2.0, 0.0, 3.0]) < 1e-12);
        assert!(geom::dist(arc.derivative(0.0), [0.0, -2.0, 0.0]) < 1e-12);
        assert!(mapped(&RawCurve2d::Unresolved, &surface, (1.0, 1.0)).is_none());
        assert!(
            mapped(
                &RawCurve2d::Line {
                    point: [0.0; 2],
                    direction: [1.0, 1.0],
                    magnitude: 1.0,
                },
                &surface,
                (1.0, 1.0)
            )
            .is_none()
        );
    }

    #[test]
    fn file_conic_pcurves_and_their_master_flag_survive_a_brep_round_trip() {
        use std::io::Read;
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures/patches.json.gz");
        let mut text = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        let fixtures: serde_json::Value = serde_json::from_str(&text).unwrap();
        let input = fixtures["patches"][4]["patch"].as_str().unwrap();
        let part = crate::read_step(input.as_bytes()).unwrap();
        assert!(
            part.file_pcurves(0)
                .iter()
                .any(|(_, p)| matches!(p, RawCurve2d::Ellipse { .. }))
        );
        assert!(
            part.primary_pcurves(0)
                .any(|(_, p)| matches!(p, RawCurve2d::Ellipse { .. }))
        );
        let (bytes, _) =
            crate::brep_write::write_step(&part, &crate::brep_write::WriteOptions::default())
                .unwrap();
        let reopened = crate::read_step(&bytes).unwrap();
        assert_eq!(part.file_pcurves(0), reopened.file_pcurves(0));
        assert_eq!(
            part.primary_pcurves(0).collect::<Vec<_>>(),
            reopened.primary_pcurves(0).collect::<Vec<_>>()
        );
    }
}
