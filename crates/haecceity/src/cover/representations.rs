//! Independent STEP trim representations mapped to the held analytic surface.
//!
//! A file uncertainty is only an eligibility check for coherent endpoints. It never changes
//! the coverage budget. An unevaluable pcurve cannot be silently omitted from that check.

use super::FaceRef;
use crate::brep::{Pcurve, RawCurve2d};
use crate::geom::{self, Curve, Frame, Surface};

/// A degree-one, two-pole positive rational spline is a finite straight segment.
/// Clamp its foot to that domain instead of extending a file pcurve as an infinite line.
pub(super) fn straight_foot(curve: &Curve, p: geom::V3) -> Option<geom::V3> {
    let (origin, direction, extent) = match curve {
        Curve::Line { origin, dir } => (*origin, *dir, None),
        Curve::Nurbs(c)
            if c.degree() == 1
                && c.control_points().len() == 2
                && c.weights().iter().all(|w| w.is_finite() && *w > 0.0) =>
        {
            let p = c.control_points();
            let dir = geom::sub(p[1], p[0]);
            (p[0], dir, Some(geom::norm(dir)))
        }
        _ => return None,
    };
    let dir = geom::unit(direction)?;
    let mut t = geom::dot(geom::sub(p, origin), dir);
    if let Some(extent) = extent {
        t = t.clamp(0.0, extent);
    }
    Some(geom::add(origin, geom::scale(dir, t)))
}

pub(super) fn straight_origin(curve: &Curve) -> Option<geom::V3> {
    match curve {
        Curve::Line { origin, .. } => straight_foot(curve, *origin),
        Curve::Nurbs(c) => straight_foot(curve, c.control_points()[0]),
        _ => None,
    }
}

/// The queried face's pcurves, and an explicitly preferred pcurve on an adjacent face.
/// The 3D curve is held separately by the edge. A preferred representation can have another
/// parameter speed, but STEP requires its sense to agree with the 3D representation.
pub(super) fn alternatives((part, face): FaceRef<'_>, edge: usize) -> Option<Vec<Curve>> {
    bounded_alternatives((part, face), edge)?
        .into_iter()
        .map(|r| match r {
            Representation::Exact(curve) => Some(curve),
            Representation::Hull(_) => None,
        })
        .collect()
}

pub(super) enum Representation {
    Exact(Curve),
    Hull(Vec<geom::V3>),
}

/// Bound the image of a finite positive pcurve in the native spline parameterization.
/// Converted swept parameterizations and unbounded parameter lines cannot use this fallback.
fn bounded(raw: &RawCurve2d, surface: &Surface, units: (f64, f64)) -> Option<Representation> {
    if let Some(curve) = mapped(raw, surface, units) {
        return Some(Representation::Exact(curve));
    }
    let Surface::Freeform {
        surface,
        kind: "BSPLINE" | "BEZIER",
    } = surface
    else {
        return None;
    };
    let identity = Frame {
        origin: [0.0; 3],
        x: [1.0, 0.0, 0.0],
        y: [0.0, 1.0, 0.0],
        z: [0.0, 0.0, 1.0],
    };
    let Curve::Nurbs(pc) = mapped(raw, &Surface::Plane { frame: identity }, (1.0, 1.0))? else {
        return None;
    };
    Some(Representation::Hull(surface_hull(surface, &pc)?))
}

fn surface_hull(
    surface: &crate::nurbs::NurbsSurface,
    pc: &crate::nurbs::NurbsCurve,
) -> Option<Vec<geom::V3>> {
    let mut uv_low = [f64::INFINITY; 2];
    let mut uv_high = [f64::NEG_INFINITY; 2];
    for p in pc.control_points() {
        if p.iter().any(|x| !x.is_finite()) {
            return None;
        }
        for k in 0..2 {
            uv_low[k] = uv_low[k].min(p[k]);
            uv_high[k] = uv_high[k].max(p[k]);
        }
    }
    let mut variation = 1.0;
    for (k, degree, knots) in [
        (0, surface.degree_u, &surface.knots_u),
        (1, surface.degree_v, &surface.knots_v),
    ] {
        let (a, b) = (knots[degree], knots[knots.len() - degree - 1]);
        if ![a, b].into_iter().all(f64::is_finite) || a >= b {
            return None;
        }
        let distance = (a - uv_low[k]).max(uv_high[k] - b).max(0.0);
        if distance == 0.0 {
            continue;
        }
        if !knots[..=degree].iter().all(|t| *t == a)
            || !knots[knots.len() - degree - 1..].iter().all(|t| *t == b)
        {
            return None;
        }
        let span = knots
            .windows(2)
            .filter_map(|w| (w[1] > w[0]).then_some(w[1] - w[0]))
            .fold(f64::INFINITY, f64::min);
        if !span.is_finite() || distance > 0.01 * span {
            return None;
        }
        // On an end Bezier span, sum |B_i(s)| = (|s|+|1-s|)^degree. Knot insertion
        // keeps every homogeneous control in the original positive control hull.
        variation *= (1.0 + 2.0 * distance / span).powi(degree as i32)
            * (1.0 + 64.0 * f64::EPSILON * (degree + 1) as f64);
    }
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    let mut w_low = f64::INFINITY;
    let mut w_high = 0.0_f64;
    let count: usize = surface.control_points.iter().map(Vec::len).sum();
    if count == 0 || count > 16384 {
        return None;
    }
    for (points, weights) in surface.control_points.iter().zip(&surface.weights) {
        for (&p, &w) in points.iter().zip(weights) {
            if p.into_iter().any(|x| !x.is_finite()) || !w.is_finite() || w <= 0.0 {
                return None;
            }
            w_low = w_low.min(w);
            w_high = w_high.max(w);
            for k in 0..3 {
                low[k] = low[k].min(p[k]);
                high[k] = high[k].max(p[k]);
            }
        }
    }
    let round = 128.0 * f64::EPSILON * (surface.degree_u + surface.degree_v + 2).pow(2) as f64;
    let negative = 0.5 * (variation - 1.0).max(0.0) + round;
    let denominator = w_low - negative * (w_high - w_low) - round * w_high;
    if !variation.is_finite() || !denominator.is_finite() || denominator <= 0.0 {
        return None;
    }
    // Signed normalized coefficients have total one. Their negative mass is at most
    // negative*w_high/denominator, so each coordinate extends by at most that times its range.
    let expansion = negative * w_high / denominator;
    for k in 0..3 {
        let pad = expansion * (high[k] - low[k]) + round * low[k].abs().max(high[k].abs()).max(1.0);
        low[k] -= pad;
        high[k] += pad;
        if !low[k].is_finite() || !high[k].is_finite() {
            return None;
        }
    }
    Some(
        (0..8)
            .map(|i| std::array::from_fn(|k| if i & (1 << k) == 0 { low[k] } else { high[k] }))
            .collect(),
    )
}

pub(super) fn bounded_alternatives(
    (part, face): FaceRef<'_>,
    edge: usize,
) -> Option<Vec<Representation>> {
    let mut out = Vec::new();
    let units = part.file_units().unwrap_or((1.0, 1.0));
    for &(index, ref curve) in part.file_pcurves(face) {
        if index == edge {
            out.push(bounded(curve, &part.faces[face].surface, units)?);
        }
    }
    if part.file_pcurves(face).is_empty() {
        for (index, curve) in &part.faces[face].pcurves {
            if *index == edge {
                let Pcurve::Line { point, dir } = curve else {
                    return None;
                };
                out.push(bounded(
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
                out.push(bounded(curve, &part.faces[adjacent].surface, units)?);
            }
        }
    }
    Some(out)
}

/// A monotone, two-pole horizontal cylindrical pcurve describes a finite circular arc.
/// Quadratic rational Bezier spans retain its domain and winding; an untrimmed circle
/// would lose that evidence. The parameter speed need not match the source spline.
fn cylindrical_arc(surface: &Surface, first: [f64; 2], last: [f64; 2]) -> Option<Curve> {
    let Surface::Cylinder { frame, radius } = surface else {
        return None;
    };
    let angle = last[0] - first[0];
    if !first.into_iter().chain(last).all(f64::is_finite)
        || !radius.is_finite()
        || *radius <= 0.0
        || first[1] != last[1]
        || angle == 0.0
        || angle.abs() > std::f64::consts::TAU
    {
        return None;
    }
    let count = (angle.abs() / std::f64::consts::FRAC_PI_2).ceil() as usize;
    let mut points = Vec::new();
    let mut weights = Vec::new();
    let mut knots = vec![0.0; 3];
    for i in 0..count {
        let a = first[0] + angle * i as f64 / count as f64;
        let b = first[0] + angle * (i + 1) as f64 / count as f64;
        let middle = 0.5 * (a + b);
        let weight = (0.5 * (b - a)).cos();
        if i == 0 {
            points.push(surface.value(a, first[1]));
            weights.push(1.0);
        }
        points.extend([
            frame.to_world([
                radius * middle.cos() / weight,
                radius * middle.sin() / weight,
                first[1],
            ]),
            surface.value(b, first[1]),
        ]);
        weights.extend([weight, 1.0]);
        knots.extend(std::iter::repeat_n((i + 1) as f64 / count as f64, 2));
    }
    knots.push(1.0);
    crate::nurbs::NurbsCurve::new(2, points, weights, knots).map(Curve::Nurbs)
}

pub(super) fn mapped(raw: &RawCurve2d, surface: &Surface, (mm, rad): (f64, f64)) -> Option<Curve> {
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
            if let RawCurve2d::BSpline {
                degree: 1,
                points,
                multiplicities,
                knots,
                weights,
            } = raw
                && points.len() == 2
                && multiplicities == &[2, 2]
                && knots.len() == 2
                && knots.iter().all(|x| x.is_finite())
                && knots[0] < knots[1]
                && weights
                    .as_ref()
                    .is_none_or(|w| w.len() == 2 && w.iter().all(|x| x.is_finite() && *x > 0.0))
            {
                if points[0][1] == points[1][1] && points[0][0] != points[1][0] {
                    return cylindrical_arc(
                        surface,
                        [points[0][0] * rad, points[0][1] * mm],
                        [points[1][0] * rad, points[1][1] * mm],
                    );
                }
                if points[0][0] != points[1][0] {
                    return None;
                }
                return crate::nurbs::NurbsCurve::new(
                    1,
                    points
                        .iter()
                        .map(|p| surface.value(p[0] * rad, p[1] * mm))
                        .collect(),
                    weights.clone().unwrap_or_else(|| vec![1.0; 2]),
                    vec![knots[0], knots[0], knots[1], knots[1]],
                )
                .map(Curve::Nurbs);
            }
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
    fn horizontal_cylinder_pcurves_keep_their_finite_arc_and_sense() {
        let surface = Surface::Cylinder {
            frame: frame(),
            radius: 71.5,
        };
        let reference = Curve::Circle {
            frame: Frame {
                origin: [0.0, 0.0, 60.0],
                ..frame()
            },
            radius: 71.5,
        };
        for (a, b) in [(180.0, 360.0), (360.0, 180.0)] {
            let pc = RawCurve2d::BSpline {
                degree: 1,
                points: vec![[a, 30.0], [b, 30.0]],
                multiplicities: vec![2, 2],
                knots: vec![3.0, 7.0],
                weights: Some(vec![1.0, 3.0]),
            };
            let rad = std::f64::consts::PI / 180.0;
            let image = mapped(&pc, &surface, (2.0, rad)).unwrap();
            let bound = super::super::conic_bounds::deviation(
                &reference,
                (a * rad, b * rad),
                &image,
                (0.0, 1.0),
            )
            .unwrap();
            assert!(bound < 1e-7);
            assert!(
                super::super::conic_bounds::deviation(
                    &reference,
                    (b * rad, a * rad),
                    &image,
                    (0.0, 1.0),
                )
                .is_none()
            );
        }
        for (a, b) in [([0.0, 0.0], [0.1, 1.0]), ([0.0, 0.0], [7.0, 0.0])] {
            let pc = RawCurve2d::BSpline {
                degree: 1,
                points: vec![a, b],
                multiplicities: vec![2, 2],
                knots: vec![0.0, 1.0],
                weights: None,
            };
            assert!(mapped(&pc, &surface, (1.0, 1.0)).is_none());
        }
    }

    #[test]
    fn clamped_native_surface_extensions_remain_enclosed() {
        use crate::nurbs::{NurbsCurve, NurbsSurface};
        let surface = NurbsSurface::new(
            1,
            1,
            vec![
                vec![[0.0, 0.0, 0.0], [0.0, 1.0, 2.0]],
                vec![[1.0, 0.0, 1.0], [1.0, 1.0, 3.0]],
            ],
            vec![vec![1.0; 2]; 2],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
        )
        .unwrap();
        let pc = NurbsCurve::new(
            1,
            vec![[-1e-6, -2e-6, 0.0], [1.0 + 1e-6, 1.0 + 2e-6, 0.0]],
            vec![1.0; 2],
            vec![0.0, 0.0, 1.0, 1.0],
        )
        .unwrap();
        let hull = surface_hull(&surface, &pc).unwrap();
        // This independent bilinear surface is exactly (u,v,u+2v), including its extension.
        for p in [[-1e-6, -2e-6, -5e-6], [1.0 + 1e-6, 1.0 + 2e-6, 3.0 + 5e-6]] {
            for k in 0..3 {
                assert!(hull.iter().map(|q| q[k]).fold(f64::INFINITY, f64::min) <= p[k]);
                assert!(hull.iter().map(|q| q[k]).fold(f64::NEG_INFINITY, f64::max) >= p[k]);
            }
        }
        let far = NurbsCurve::new(
            1,
            vec![[-1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            vec![1.0; 2],
            vec![0.0, 0.0, 1.0, 1.0],
        )
        .unwrap();
        assert!(surface_hull(&surface, &far).is_none());
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
            crate::brep_write::write_step(&part, &crate::brep_write::WriteOptions::new("coverage"))
                .unwrap();
        let reopened = crate::read_step(&bytes).unwrap();
        assert_eq!(part.file_pcurves(0), reopened.file_pcurves(0));
        assert_eq!(
            part.primary_pcurves(0).collect::<Vec<_>>(),
            reopened.primary_pcurves(0).collect::<Vec<_>>()
        );
    }

    #[test]
    fn straight_rational_pcurves_keep_their_finite_domain() {
        let curve = mapped(
            &RawCurve2d::BSpline {
                degree: 1,
                points: vec![[0.0, 0.0], [1.0, 0.0]],
                multiplicities: vec![2, 2],
                knots: vec![0.0, 1.0],
                weights: Some(vec![1.0, 3.0]),
            },
            &Surface::Plane { frame: frame() },
            (1.0, 1.0),
        )
        .unwrap();
        assert_eq!(
            straight_foot(&curve, [2.0, 0.0, 0.0]),
            Some([1.0, 0.0, 0.0])
        );
        assert_eq!(
            straight_foot(&curve, [-1.0, 0.0, 0.0]),
            Some([0.0, 0.0, 0.0])
        );
        assert_eq!(
            straight_foot(&curve, [0.5, 1.0, 0.0]),
            Some([0.5, 0.0, 0.0])
        );
        let cylinder = mapped(
            &RawCurve2d::BSpline {
                degree: 1,
                points: vec![[0.0, 0.0], [0.0, 2.0]],
                multiplicities: vec![2, 2],
                knots: vec![0.0, 1.0],
                weights: None,
            },
            &Surface::Cylinder {
                frame: frame(),
                radius: 1.0,
            },
            (1.0, 1.0),
        )
        .unwrap();
        assert_eq!(
            straight_foot(&cylinder, [1.0, 0.0, 3.0]),
            Some([1.0, 0.0, 2.0])
        );
    }
}
