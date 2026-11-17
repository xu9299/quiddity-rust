//! Sufficient finite containment proofs for a curved target in a convex planar support,
//! and for local cutouts in a full cylindrical band. Neither path proves a gap.

use super::{FaceRef, conic_bounds, linear, representations, same_surface};
use crate::brep::{Edge, Loop};
use crate::geom::{self, Curve, Frame, Surface, V3};
use crate::sampling::{edge_interval, extremes_along};

fn points(edge: &Edge, curve: &Curve, direction: V3) -> Option<Vec<V3>> {
    let mut out = vec![edge.start, edge.end];
    if let Some(a) = representations::straight_foot(curve, edge.start) {
        out.extend([a, representations::straight_foot(curve, edge.end)?]);
    } else if let Curve::Nurbs(c) = curve {
        if !c.weights().iter().all(|w| w.is_finite() && *w > 0.0) {
            return None;
        }
        out.extend(c.control_points());
    } else {
        let interval = edge_interval(
            curve,
            edge.start,
            edge.end,
            edge.same_sense,
            edge.vertices.0 == edge.vertices.1,
        );
        out.extend([curve.value(interval.0), curve.value(interval.1)]);
        out.extend(extremes_along(curve, interval, &[direction]));
    }
    out.iter().flatten().all(|x| x.is_finite()).then_some(out)
}

pub(super) fn covered(patch: FaceRef<'_>, supports: &[FaceRef<'_>]) -> Option<bool> {
    supports
        .iter()
        .copied()
        .filter(|s| same_surface(&patch.0.faces[patch.1].surface, &s.0.faces[s.1].surface))
        .find_map(|support| planar(patch, support).or_else(|| cylindrical(patch, support)))
}

fn planar(patch: FaceRef<'_>, (part, face): FaceRef<'_>) -> Option<bool> {
    let frame = linear::canonical_frame((part, face))?;
    let fc = &part.faces[face];
    if fc.loops.len() != 1 {
        return None;
    }
    let lp = &fc.loops[0];
    if lp.edges.len() < 3 || lp.edges.len() > 512 {
        return None;
    }
    let mut vertices = Vec::new();
    let mut radius = 0.0_f64;
    let mut size = frame.origin.into_iter().map(f64::abs).fold(1.0, f64::max);
    let mut first = None;
    let mut last = None;
    for &(index, forward) in &lp.edges {
        if part.unresolved_edges().contains(&index) {
            return None;
        }
        let e = &part.edges[index];
        let (a, b) = if forward {
            (e.start, e.end)
        } else {
            (e.end, e.start)
        };
        if last.is_some_and(|p| p != a) {
            return None;
        }
        first.get_or_insert(a);
        last = Some(b);
        vertices.push(a);
        let mut curves = vec![e.curve.clone()];
        curves.extend(representations::alternatives((part, face), index)?);
        for curve in curves {
            for p in [a, b] {
                let q = representations::straight_foot(&curve, p)?;
                radius = radius.max(geom::dist(p, q));
                size = p.into_iter().chain(q).map(f64::abs).fold(size, f64::max);
            }
        }
    }
    if first != last {
        return None;
    }
    let project = |p| {
        let d = geom::sub(p, frame.origin);
        (geom::dot(d, frame.x), geom::dot(d, frame.y))
    };
    let polygon: Vec<_> = vertices.iter().copied().map(project).collect();
    let cross = |a: (f64, f64), b: (f64, f64)| a.0 * b.1 - a.1 * b.0;
    let difference = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0, a.1 - b.1);
    let signed: f64 = (0..polygon.len())
        .map(|i| cross(polygon[i], polygon[(i + 1) % polygon.len()]))
        .sum();
    if !signed.is_finite() || signed == 0.0 {
        return None;
    }
    let sign = signed.signum();
    for i in 0..polygon.len() {
        let (a, b, c) = (
            polygon[i],
            polygon[(i + 1) % polygon.len()],
            polygon[(i + 2) % polygon.len()],
        );
        if sign * cross(difference(b, a), difference(c, b)) <= 0.0 {
            return None;
        }
        for &p in &polygon {
            if sign * cross(difference(b, a), difference(p, a)) < 0.0 {
                return None;
            }
        }
    }
    let target = &patch.0.faces[patch.1];
    if target.loops.is_empty() {
        return None;
    }
    for i in 0..polygon.len() {
        let d = difference(polygon[(i + 1) % polygon.len()], polygon[i]);
        let length = d.0.hypot(d.1);
        if !length.is_finite() || length == 0.0 {
            return None;
        }
        let inward = geom::add(
            geom::scale(frame.x, -sign * d.1 / length),
            geom::scale(frame.y, sign * d.0 / length),
        );
        let origin = vertices[i];
        for lp in &target.loops {
            if lp.edges.is_empty() {
                return None;
            }
            let mut first = None;
            let mut last = None;
            for &(index, forward) in &lp.edges {
                if patch.0.unresolved_edges().contains(&index) {
                    return None;
                }
                let e = &patch.0.edges[index];
                let (a, b) = if forward {
                    (e.start, e.end)
                } else {
                    (e.end, e.start)
                };
                if last.is_some_and(|p| p != a) {
                    return None;
                }
                first.get_or_insert(a);
                last = Some(b);
                let mut curves = vec![e.curve.clone()];
                curves.extend(representations::alternatives(patch, index)?);
                for curve in curves {
                    for p in points(e, &curve, inward)? {
                        let magnitude = p.into_iter().map(f64::abs).fold(size, f64::max);
                        let padding = 256.0 * f64::EPSILON * magnitude;
                        let clearance = geom::dot(geom::sub(p, origin), inward);
                        if !clearance.is_finite() || clearance <= radius + padding {
                            return None;
                        }
                    }
                }
            }
            if first != last {
                return None;
            }
        }
    }
    Some(true)
}

struct Band {
    loop_index: usize,
    outer_low: f64,
    inner_low: f64,
    inner_high: f64,
    outer_high: f64,
    seam_error: f64,
}

fn range(edge: &Edge, curve: &Curve, frame: &Frame) -> Option<(f64, f64)> {
    let samples = points(edge, curve, frame.z)?;
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for p in samples {
        let v = geom::dot(geom::sub(p, frame.origin), frame.z);
        low = low.min(v);
        high = high.max(v);
    }
    let a = curve.value(curve.parameter(edge.start));
    let b = curve.value(curve.parameter(edge.end));
    let magnitude = edge
        .start
        .into_iter()
        .chain(edge.end)
        .chain(a)
        .chain(b)
        .map(f64::abs)
        .fold(1.0, f64::max);
    let error =
        geom::dist(a, edge.start).max(geom::dist(b, edge.end)) + 256.0 * f64::EPSILON * magnitude;
    (low.is_finite() && high.is_finite() && error.is_finite())
        .then_some((low - error, high + error))
}

fn closed_contour(
    (part, face): FaceRef<'_>,
    index: usize,
    forward: bool,
    frame: &Frame,
) -> Option<((f64, f64), f64)> {
    let edge = &part.edges[index];
    if edge.start != edge.end || edge.vertices.0 != edge.vertices.1 {
        return None;
    }
    let mut curves = vec![edge.curve.clone()];
    curves.extend(representations::alternatives((part, face), index)?);
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    let mut winding = None;
    for curve in curves {
        if matches!(curve, Curve::Nurbs(_)) {
            let interval = |c| edge_interval(c, edge.start, edge.end, edge.same_sense, true);
            let d = conic_bounds::deviation(
                &edge.curve,
                interval(&edge.curve),
                &curve,
                interval(&curve),
            )?;
            let Curve::Circle { radius, .. } = edge.curve else {
                return None;
            };
            if d >= 0.25 * radius {
                return None;
            }
            let r = range(edge, &edge.curve, frame)?;
            let axial = d * geom::norm(frame.z);
            low = low.min(r.0 - axial);
            high = high.max(r.1 + axial);
            continue;
        }
        let (circle, radius) = match curve {
            Curve::Circle { frame, radius } => (frame, radius),
            // A planar STEP circle maps to this conic representation. Exact equal
            // axes preserve the entire circle, its sense and every finite endpoint.
            Curve::Ellipse {
                frame,
                major,
                minor,
            } if major == minor => (frame, major),
            _ => return None,
        };
        let axis = geom::dot(circle.z, frame.z);
        let off = geom::sub(circle.origin, frame.origin);
        let across = geom::sub(off, geom::scale(frame.z, geom::dot(off, frame.z)));
        if radius <= 0.0 || axis.abs() < 1.0 - 1e-12 || geom::norm(across) > 0.25 * radius {
            return None;
        }
        let sense = axis.signum()
            * (if edge.same_sense { 1.0 } else { -1.0 })
            * (if forward { 1.0 } else { -1.0 });
        if winding.is_some_and(|s| s != sense) {
            return None;
        }
        winding = Some(sense);
        let r = range(edge, &curve, frame)?;
        low = low.min(r.0);
        high = high.max(r.1);
    }
    Some(((low, high), winding?))
}

fn band((part, face): FaceRef<'_>, frame: &Frame, radius: f64) -> Option<Band> {
    for (loop_index, lp) in part.faces[face].loops.iter().enumerate() {
        if lp.edges.len() != 4 {
            continue;
        }
        let mut counts = std::collections::BTreeMap::<usize, Vec<bool>>::new();
        for &(i, f) in &lp.edges {
            counts.entry(i).or_default().push(f);
        }
        let Some((&seam, walk)) = counts.iter().find(|(_, v)| v.len() == 2 && v[0] != v[1]) else {
            continue;
        };
        if counts.len() != 3 || walk.len() != 2 {
            continue;
        }
        let mut first = None;
        let mut last = None;
        for &(i, f) in &lp.edges {
            if part.unresolved_edges().contains(&i) {
                return None;
            }
            let e = &part.edges[i];
            let (a, b) = if f {
                (e.start, e.end)
            } else {
                (e.end, e.start)
            };
            if last.is_some_and(|p| p != a) {
                return None;
            }
            first.get_or_insert(a);
            last = Some(b);
        }
        if first != last {
            return None;
        }
        let mut contours = Vec::new();
        let mut winding = 0.0;
        for (&i, f) in &counts {
            if i == seam {
                continue;
            }
            if f.len() != 1 {
                return None;
            }
            let (r, w) = closed_contour((part, face), i, f[0], frame)?;
            contours.push(r);
            winding += w;
        }
        if contours.len() != 2 || winding != 0.0 {
            return None;
        }
        contours.sort_by(|a, b| (a.0 + a.1).total_cmp(&(b.0 + b.1)));
        let e = &part.edges[seam];
        let a = representations::straight_foot(&e.curve, e.start)?;
        let b = representations::straight_foot(&e.curve, e.end)?;
        let across = |p| {
            let d = geom::sub(p, frame.origin);
            geom::sub(d, geom::scale(frame.z, geom::dot(d, frame.z)))
        };
        let (u, v) = (across(a), across(b));
        let d = geom::sub(v, u);
        let t = if geom::dot(d, d) > 0.0 {
            (-geom::dot(u, d) / geom::dot(d, d)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        if geom::norm(geom::add(u, geom::scale(d, t))) < 0.75 * radius {
            return None;
        }
        let mut displacement = geom::dist(a, e.start).max(geom::dist(b, e.end));
        for c in representations::alternatives((part, face), seam)? {
            displacement = displacement
                .max(geom::dist(a, representations::straight_foot(&c, e.start)?))
                .max(geom::dist(b, representations::straight_foot(&c, e.end)?));
        }
        let magnitude = a.into_iter().chain(b).map(f64::abs).fold(1.0, f64::max);
        displacement += 256.0 * f64::EPSILON * magnitude;
        if displacement >= 0.25 * radius {
            return None;
        }
        // At rho >= R/2 the physical parameter map has operator norm at most two.
        // The seam uncertainty remains an area charge, even though its topology is shared.
        let seam_error = 8.0 * displacement * geom::dist(a, b)
            + 4.0 * std::f64::consts::PI * displacement * displacement;
        return Some(Band {
            loop_index,
            outer_low: contours[0].0,
            inner_low: contours[0].1,
            inner_high: contours[1].0,
            outer_high: contours[1].1,
            seam_error,
        });
    }
    None
}

fn cosine_low(a: f64, b: f64, low: f64, high: f64) -> f64 {
    let amplitude = a.hypot(b);
    if high - low >= std::f64::consts::TAU {
        return -amplitude;
    }
    let phase = b.atan2(a);
    let start = geom::nearest_turn(low, phase);
    let stop = start + (high - low);
    let mut minimum = (a * start.cos() + b * start.sin()).min(a * stop.cos() + b * stop.sin());
    for k in -1..=1 {
        let at = phase + std::f64::consts::PI + k as f64 * std::f64::consts::TAU;
        if start <= at && at <= stop {
            minimum = -amplitude;
        }
    }
    minimum
}

fn hole((part, face): FaceRef<'_>, lp: &Loop, frame: &Frame, radius: f64) -> Option<(f64, f64)> {
    let first = part.edges[lp.edges.first()?.0].start;
    let off = geom::sub(first, frame.origin);
    let radial = geom::unit(geom::sub(
        off,
        geom::scale(frame.z, geom::dot(off, frame.z)),
    ))?;
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    let mut last = None;
    let mut start = None;
    for &(index, forward) in &lp.edges {
        if part.unresolved_edges().contains(&index) {
            return None;
        }
        let e = &part.edges[index];
        let (a, b) = if forward {
            (e.start, e.end)
        } else {
            (e.end, e.start)
        };
        if last.is_some_and(|p| p != a) {
            return None;
        }
        start.get_or_insert(a);
        last = Some(b);
        let r = range(e, &e.curve, frame)?;
        low = low.min(r.0);
        high = high.max(r.1);
        let safely_inside = |p: V3| {
            let magnitude = p
                .into_iter()
                .chain(frame.origin)
                .map(f64::abs)
                .fold(1.0, f64::max);
            let projection = geom::dot(geom::sub(p, frame.origin), radial);
            projection.is_finite() && projection > 0.25 * radius + 512.0 * f64::EPSILON * magnitude
        };
        if !points(e, &e.curve, radial)?.into_iter().all(safely_inside) {
            return None;
        }
        if part.file_pcurves(face).is_empty() {
            for c in representations::alternatives((part, face), index)? {
                let r = range(e, &c, frame)?;
                low = low.min(r.0);
                high = high.max(r.1);
                if !points(e, &c, radial)?.into_iter().all(safely_inside) {
                    return None;
                }
            }
        }
        for &f in &part.edge_faces()[index] {
            let raw: Vec<_> = if f == face {
                part.file_pcurves(f)
                    .iter()
                    .filter(|(i, _)| *i == index)
                    .collect()
            } else {
                part.primary_pcurves(f)
                    .filter(|(i, _)| *i == index)
                    .collect()
            };
            for (_, pc) in raw {
                let Surface::Cylinder {
                    frame: g,
                    radius: r,
                } = &part.faces[f].surface
                else {
                    return None;
                };
                if let Some(c) = representations::mapped(
                    pc,
                    &part.faces[f].surface,
                    part.file_units().unwrap_or((1.0, 1.0)),
                ) {
                    let bound = range(e, &c, frame)?;
                    low = low.min(bound.0);
                    high = high.max(bound.1);
                    if !points(e, &c, radial)?.into_iter().all(safely_inside) {
                        return None;
                    }
                    continue;
                }
                let identity = Frame {
                    origin: [0.0; 3],
                    x: [1.0, 0.0, 0.0],
                    y: [0.0, 1.0, 0.0],
                    z: [0.0, 0.0, 1.0],
                };
                let c =
                    representations::mapped(pc, &Surface::Plane { frame: identity }, (1.0, 1.0))?;
                let Curve::Nurbs(c) = c else {
                    return None;
                };
                if !c.weights().iter().all(|w| w.is_finite() && *w > 0.0) {
                    return None;
                }
                let (mm, rad) = part.file_units().unwrap_or((1.0, 1.0));
                let mut uv_low = (f64::INFINITY, f64::INFINITY);
                let mut uv_high = (f64::NEG_INFINITY, f64::NEG_INFINITY);
                for p in c.control_points() {
                    uv_low = (uv_low.0.min(p[0] * rad), uv_low.1.min(p[1] * mm));
                    uv_high = (uv_high.0.max(p[0] * rad), uv_high.1.max(p[1] * mm));
                }
                let origin = geom::sub(g.origin, frame.origin);
                let slope = geom::dot(g.z, frame.z);
                let noise = r * (geom::dot(g.x, frame.z).hypot(geom::dot(g.y, frame.z)));
                let offset = geom::dot(origin, frame.z);
                let a = offset + slope * uv_low.1;
                let b = offset + slope * uv_high.1;
                let magnitude = g
                    .origin
                    .into_iter()
                    .chain(frame.origin)
                    .chain([r * uv_low.0, r * uv_high.0, uv_low.1, uv_high.1])
                    .map(f64::abs)
                    .fold(1.0, f64::max);
                let padding = 1024.0 * f64::EPSILON * magnitude;
                if ![a, b, noise, padding].into_iter().all(f64::is_finite) {
                    return None;
                }
                low = low.min(a.min(b) - noise - padding);
                high = high.max(a.max(b) + noise + padding);
                let h0 = geom::dot(g.z, radial) * uv_low.1;
                let h1 = geom::dot(g.z, radial) * uv_high.1;
                let radial_low = geom::dot(origin, radial)
                    + h0.min(h1)
                    + cosine_low(
                        r * geom::dot(g.x, radial),
                        r * geom::dot(g.y, radial),
                        uv_low.0,
                        uv_high.0,
                    );
                if !radial_low.is_finite() || radial_low <= 0.25 * radius + padding {
                    return None;
                }
            }
        }
    }
    if start != last {
        return None;
    }
    let padding = 512.0 * f64::EPSILON * low.abs().max(high.abs()).max(1.0);
    Some((low - padding, high + padding))
}

fn cylindrical(patch: FaceRef<'_>, support: FaceRef<'_>) -> Option<bool> {
    let Surface::Cylinder { frame, radius } = &patch.0.faces[patch.1].surface else {
        return None;
    };
    if support.0.faces[support.1].loops.len() != 1 || *radius <= 0.0 {
        return None;
    }
    let source = band(support, frame, *radius)?;
    let target = band(patch, frame, *radius)?;
    let circumference = std::f64::consts::TAU * radius;
    let mut lower = circumference * (target.inner_high - target.inner_low) - target.seam_error;
    for (i, lp) in patch.0.faces[patch.1].loops.iter().enumerate() {
        if i == target.loop_index {
            continue;
        }
        let (a, b) = hole(patch, lp, frame, *radius)?;
        if a < target.outer_low || b > target.outer_high {
            return None;
        }
        lower -= circumference * (b - a);
    }
    if !lower.is_finite() || lower <= 0.0 || source.inner_low >= source.inner_high {
        return None;
    }
    let missing = circumference
        * ((source.inner_low - target.outer_low).max(0.0)
            + (target.outer_high - source.inner_high).max(0.0))
        + source.seam_error;
    (missing.is_finite() && missing <= 1e-9 * lower).then_some(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::{Face, Part, RawCurve2d};
    use crate::cover::{Coverage, assess};
    use crate::nurbs::NurbsCurve;

    fn frame() -> Frame {
        Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        }
    }

    fn edge(curve: Curve, start: V3, end: V3, a: usize, b: usize) -> Edge {
        Edge {
            curve,
            start,
            end,
            vertices: (a, b),
            same_sense: true,
            samples: vec![],
        }
    }

    fn part(surface: Surface, edges: Vec<Edge>, loops: Vec<Loop>) -> Part {
        Part::new(
            vec![Face {
                surface,
                reversed: false,
                solid: None,
                pcurves: vec![],
                loops,
            }],
            edges,
            vec![],
        )
    }

    fn rectangle() -> Part {
        let points = [
            [-2.0, -2.0, 0.0],
            [2.0, -2.0, 0.0],
            [2.0, 2.0, 0.0],
            [-2.0, 2.0, 0.0],
        ];
        let edges = (0..4)
            .map(|i| {
                let (a, b) = (points[i], points[(i + 1) % 4]);
                edge(
                    Curve::Line {
                        origin: a,
                        dir: geom::unit(geom::sub(b, a)).unwrap(),
                    },
                    a,
                    b,
                    i,
                    (i + 1) % 4,
                )
            })
            .collect();
        part(
            Surface::Plane { frame: frame() },
            edges,
            vec![Loop {
                edges: (0..4).map(|i| (i, true)).collect(),
                vertex: None,
            }],
        )
    }

    fn ellipse(major: f64) -> Part {
        let curve = Curve::Ellipse {
            frame: frame(),
            major,
            minor: 1.0,
        };
        // The seam stays inside the rectangle even when the finite curve protrudes.
        let p = curve.value(std::f64::consts::FRAC_PI_2);
        part(
            Surface::Plane { frame: frame() },
            vec![edge(curve, p, p, 0, 0)],
            vec![Loop {
                edges: vec![(0, true)],
                vertex: None,
            }],
        )
    }

    #[test]
    fn convex_containment_checks_the_entire_finite_curve() {
        let support = rectangle();
        assert_eq!(
            assess((&ellipse(1.0), 0), &[(&support, 0)]),
            Coverage::Covered
        );
        assert_ne!(
            assess((&ellipse(2.1), 0), &[(&support, 0)]),
            Coverage::Covered
        );
        let mut holed = rectangle();
        holed.edges.push(ellipse(0.5).edges[0].clone());
        holed.faces[0].loops.push(Loop {
            edges: vec![(4, false)],
            vertex: None,
        });
        let holed = Part::new(holed.faces.clone(), holed.edges.clone(), vec![]);
        assert_ne!(
            assess((&ellipse(1.0), 0), &[(&holed, 0)]),
            Coverage::Covered
        );
    }

    #[test]
    fn containment_cannot_drop_a_conflicting_pcurve() {
        let patch = ellipse(1.0).with_file_pcurves(
            std::collections::BTreeMap::from([(
                0,
                vec![(
                    0,
                    RawCurve2d::Ellipse {
                        point: [0.0; 2],
                        reference: None,
                        semi_axis_1: 3.0,
                        semi_axis_2: 1.0,
                    },
                )],
            )]),
            Default::default(),
            (1.0, 1.0),
            Some(1e-7),
        );
        assert_ne!(assess((&patch, 0), &[(&rectangle(), 0)]), Coverage::Covered);
    }

    fn cylinder(height: f64, with_hole: bool) -> Part {
        let arc = |z| Curve::Circle {
            frame: Frame {
                origin: [0.0, 0.0, z],
                ..frame()
            },
            radius: 2.0,
        };
        let (a, b) = ([2.0, 0.0, 0.0], [2.0, 0.0, height]);
        let mut edges = vec![
            edge(arc(0.0), a, a, 0, 0),
            edge(
                Curve::Line {
                    origin: a,
                    dir: [0.0, 0.0, 1.0],
                },
                a,
                b,
                0,
                1,
            ),
            edge(arc(height), b, b, 1, 1),
        ];
        let mut loops = vec![Loop {
            edges: vec![(0, false), (1, true), (2, true), (1, false)],
            vertex: None,
        }];
        if with_hole {
            let points = vec![
                [2.0, 0.0, 1.0],
                [2.0, 0.25, 1.0],
                [2.0, 0.25, 2.0],
                [2.0, 0.0, 2.0],
                [2.0, 0.0, 1.0],
            ];
            let curve = Curve::Nurbs(
                NurbsCurve::new(
                    1,
                    points.clone(),
                    vec![1.0; 5],
                    vec![0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0],
                )
                .unwrap(),
            );
            edges.push(edge(curve, points[0], points[0], 2, 2));
            loops.push(Loop {
                edges: vec![(3, true)],
                vertex: None,
            });
        }
        part(
            Surface::Cylinder {
                frame: frame(),
                radius: 2.0,
            },
            edges,
            loops,
        )
    }

    #[test]
    fn full_cylindrical_bands_cover_local_holes_but_keep_their_height() {
        let patch = cylinder(3.0, true);
        assert_eq!(
            assess((&patch, 0), &[(&cylinder(3.0, false), 0)]),
            Coverage::Covered
        );
        assert_ne!(
            assess((&patch, 0), &[(&cylinder(2.9, false), 0)]),
            Coverage::Covered
        );
        assert_ne!(
            assess((&cylinder(3.0, false), 0), &[(&patch, 0)]),
            Coverage::Covered
        );
    }

    #[test]
    fn equal_axis_contours_keep_the_circle_proof_without_fitting_other_ellipses() {
        for minor in [2.0, 2.0_f64.next_up(), 3.0] {
            let mut part = cylinder(3.0, false);
            part.edges[0].curve = Curve::Ellipse {
                frame: frame(),
                major: 2.0,
                minor,
            };
            assert_eq!(
                closed_contour((&part, 0), 0, false, &frame()).is_some(),
                minor == 2.0
            );
        }
    }

    #[test]
    fn cylindrical_hole_pcurves_cannot_leave_the_outer_band() {
        let patch = cylinder(3.0, true).with_file_pcurves(
            std::collections::BTreeMap::from([(
                0,
                vec![(
                    3,
                    RawCurve2d::BSpline {
                        degree: 1,
                        points: vec![[0.0, 1.0], [0.1, 1.0], [0.1, 4.0], [0.0, 4.0], [0.0, 1.0]],
                        multiplicities: vec![2, 1, 1, 1, 2],
                        knots: vec![0.0, 1.0, 2.0, 3.0, 4.0],
                        weights: None,
                    },
                )],
            )]),
            Default::default(),
            (1.0, 1.0),
            Some(1e-7),
        );
        assert_ne!(
            assess((&patch, 0), &[(&cylinder(3.0, false), 0)]),
            Coverage::Covered
        );
    }
}
