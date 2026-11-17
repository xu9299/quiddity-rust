//! Planar coverage where every boundary meeting the patch is a straight edge.
//!
//! A straight edge has no sampling chord error. Expanding its scan interval hides real gaps
//! in one parameter direction while retaining them in the other. Here the original boundary
//! segments are subtracted without expansion. Their crossings split the scan into bands on
//! which the uncovered length is affine, so midpoint integration gives its area. The relative
//! area tolerance is unchanged. This path does not turn a chord into a proof of a curved edge:
//! a curved support edge's replacement can change coverage only inside its conservative
//! analytic/control-hull box, including its vertices. Straight edges also account for the
//! discrepancy between their vertices and their analytic line. These bounds bracket the
//! uncovered area: only a decisive bracket uses this path, with the same acceptance budget.
//! Curved bounds above that budget, ambiguous boundaries, cylinders and oversized arrangements
//! remain indeterminate; ineligible curves do not require constructing an arrangement.

use super::representations::{straight_foot, straight_origin};
use super::{FaceRef, intervals, same_surface, subtract};
use crate::brep::Edge;
use crate::geom::{self, Curve, Frame, Surface, V3};
use crate::sampling::{edge_interval, extremes_along};

type Point = (f64, f64);
type Segment = (Point, Point);
type Loops = Vec<Vec<Point>>;

struct Boundary {
    loops: Loops,
    line_error: f64,
    straight_points: Option<Vec<Point>>,
    size: f64,
}

struct AreaBracket {
    gap: f64,
    error: f64,
    tolerance: f64,
    outer_gap_lower: Option<f64>,
}

// Bound the quadratic arrangement work; this is an eligibility limit, not a geometric epsilon.
const MAX_SEGMENTS: usize = 512;

#[derive(Clone, Copy)]
struct Box2 {
    low: Point,
    high: Point,
}

impl Box2 {
    fn empty() -> Self {
        Self {
            low: (f64::INFINITY, f64::INFINITY),
            high: (f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }

    fn add(&mut self, p: Point) -> Option<()> {
        if !p.0.is_finite() || !p.1.is_finite() {
            return None;
        }
        self.low = (self.low.0.min(p.0), self.low.1.min(p.1));
        self.high = (self.high.0.max(p.0), self.high.1.max(p.1));
        Some(())
    }

    fn disjoint(self, other: Self, pad: f64) -> bool {
        self.high.0 + pad < other.low.0
            || other.high.0 + pad < self.low.0
            || self.high.1 + pad < other.low.1
            || other.high.1 + pad < self.low.1
    }

    fn intersection_area(self, other: Self, pad: f64) -> f64 {
        let width = (self.high.0 + pad).min(other.high.0) - (self.low.0 - pad).max(other.low.0);
        let height = (self.high.1 + pad).min(other.high.1) - (self.low.1 - pad).max(other.low.1);
        width.max(0.0) * height.max(0.0)
    }
}

/// `None` means this path cannot establish the boundary representation it needs.
pub(super) fn covered(patch: FaceRef<'_>, supports: &[FaceRef<'_>]) -> Option<bool> {
    let area = finite_area(patch, supports, &[])?;
    if area.gap + area.error <= area.tolerance {
        Some(true)
    } else if area.gap > area.tolerance + area.error
        || area.outer_gap_lower.is_some_and(|g| g > area.tolerance)
    {
        Some(false)
    } else {
        None
    }
}

/// An acceptance-only union proof. A matched one-loop support can fill a loop removed
/// from another support's parity region. Charge the complete boundary correspondence;
/// integrate that reduced region alone, so the filler never pays a coarse curved box.
pub(super) fn covered_filled(patch: FaceRef<'_>, supports: &[FaceRef<'_>]) -> Option<bool> {
    for &owner in supports {
        if let Some(area) = finite_area(patch, &[owner], supports)
            && area.gap + area.error <= area.tolerance
        {
            return Some(true);
        }
    }
    None
}

/// Integrate finite analytic segments, retaining uncertainty only at their common vertices
/// and between the file's independent curve representations. No infinite lines are intersected.
pub(super) fn canonical_frame((part, face): FaceRef<'_>) -> Option<Frame> {
    let Surface::Plane { frame } = &part.faces[face].surface else {
        return None;
    };
    let mut z = geom::unit(frame.z)?;
    let k = geom::dominant_axis(z);
    if z[k] < 0.0 {
        z = geom::scale(z, -1.0);
    }
    let vertex = part.faces[face]
        .loops
        .iter()
        .flat_map(|lp| &lp.edges)
        .flat_map(|&(e, _)| [part.edges[e].start, part.edges[e].end])
        .min_by(|a, b| {
            a[0].total_cmp(&b[0])
                .then(a[1].total_cmp(&b[1]))
                .then(a[2].total_cmp(&b[2]))
        })?;
    let origin = geom::sub(
        vertex,
        geom::scale(z, geom::dot(geom::sub(vertex, frame.origin), z)),
    );
    let x = geom::default_x_axis(z);
    Some(Frame {
        origin,
        x,
        y: geom::cross(z, x),
        z,
    })
}

fn finite_area(
    patch: FaceRef<'_>,
    supports: &[FaceRef<'_>],
    fillers: &[FaceRef<'_>],
) -> Option<AreaBracket> {
    let frame = canonical_frame(patch)?;
    let target = finite_loops(patch, &frame, None, &[])?;
    let mut bounds = Box2::empty();
    for &p in target.loops.iter().flatten() {
        bounds.add(p)?;
    }
    let area = uncovered_area(&target.loops, &[], bounds)?;
    if !area.is_finite() || area <= 0.0 {
        return None;
    }
    // Apply the same relative tolerance to the bounded target area. The extra target term
    // uses the lower area bound for acceptance and the upper one for rejection.
    let mut error = target.line_error * (1.0 + 1e-9);
    let mut theirs = Vec::new();
    let mut boundaries = Vec::new();
    for &support in supports {
        if same_surface(
            &patch.0.faces[patch.1].surface,
            &support.0.faces[support.1].surface,
        ) {
            let boundary = finite_loops(support, &frame, Some(bounds), fillers)?;
            error += boundary.line_error;
            theirs.push(boundary.loops.clone());
            boundaries.push(boundary);
        }
    }
    let gap = uncovered_area(&target.loops, &theirs, bounds)?;
    let tolerance = 1e-9 * area;
    let outer_gap_lower =
        if fillers.is_empty() && gap + error > tolerance && gap <= tolerance + error {
            outside_straight_envelopes(&target, &boundaries, bounds)
        } else {
            None
        };
    Some(AreaBracket {
        gap,
        error,
        tolerance,
        outer_gap_lower,
    })
}

/// Every finite straight boundary path, including its vertex joins, lies in this hull.
/// The parity region of any closed path also lies there, even with holes or crossings.
/// Enlarging all supports can only decrease the uncovered area. This is rejection-only:
/// no area inside an envelope establishes actual support or restores a removed loop.
fn outside_straight_envelopes(
    target: &Boundary,
    supports: &[Boundary],
    bounds: Box2,
) -> Option<f64> {
    let mut envelopes = Vec::new();
    let mut error = target.line_error * (1.0 + 1e-9);
    for support in supports {
        let polygon = straight_outer(support.straight_points.as_ref()?, support.size)?;
        for ends in polygon.windows(2) {
            error += segment_error([ends[0], ends[1]], [0.0; 2], support.size, Some(bounds))?;
        }
        envelopes.push(vec![polygon]);
    }
    Some(uncovered_area(&target.loops, &envelopes, bounds)? - error)
}

fn straight_outer(points: &[Point], size: f64) -> Option<Vec<Point>> {
    // Keep adaptive orientation away from overflow and subnormal products. This only
    // limits eligibility; it is neither a geometric tolerance nor an acceptance margin.
    if points.len() > 4096
        || !size.is_finite()
        || size > 1e75
        || points.iter().any(|p| {
            [p.0, p.1]
                .into_iter()
                .any(|x| !x.is_finite() || x.abs() > 1e75 || (x != 0.0 && x.abs() < 1e-75))
        })
    {
        return None;
    }
    let pad = 64.0 * f64::EPSILON * size;
    let mut expanded = Vec::new();
    for &(x, y) in points {
        for a in [(x - pad).next_down(), (x + pad).next_up()] {
            for b in [(y - pad).next_down(), (y + pad).next_up()] {
                expanded.push((a, b));
            }
        }
    }
    expanded.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    expanded.dedup();
    let turn = |a: Point, b: Point, c: Point| {
        let coord = |p: Point| robust::Coord { x: p.0, y: p.1 };
        robust::orient2d(coord(a), coord(b), coord(c))
    };
    let half = |ordered: Vec<Point>| {
        let mut hull = Vec::new();
        for p in ordered {
            while hull.len() >= 2 && turn(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0 {
                hull.pop();
            }
            hull.push(p);
        }
        hull.pop();
        hull
    };
    let mut hull = half(expanded.clone());
    hull.extend(half(expanded.into_iter().rev().collect()));
    if hull.len() < 3 || hull.len() > MAX_SEGMENTS {
        return None;
    }
    hull.push(hull[0]);
    Some(hull)
}

fn finite_loops(
    (part, face): FaceRef<'_>,
    frame: &Frame,
    clip: Option<Box2>,
    fillers: &[FaceRef<'_>],
) -> Option<Boundary> {
    let mut contours = Vec::new();
    let mut error = 0.0;
    let mut straight_points = Some(Vec::new());
    let surface = &part.faces[face].surface;
    let mut base_size = frame.origin.into_iter().map(f64::abs).fold(1.0, f64::max);
    if let Surface::Plane { frame } = surface {
        base_size = frame
            .origin
            .into_iter()
            .map(f64::abs)
            .fold(base_size, f64::max);
    }
    for &(e, _) in part.faces[face].loops.iter().flat_map(|lp| &lp.edges) {
        for p in [part.edges[e].start, part.edges[e].end] {
            base_size = p.into_iter().map(f64::abs).fold(base_size, f64::max);
        }
    }
    let mut outer_size = base_size;
    for (loop_index, lp) in part.faces[face].loops.iter().enumerate() {
        if lp.edges.is_empty() {
            return None;
        }
        if clip.is_some()
            && !fillers.is_empty()
            && lp
                .edges
                .iter()
                .any(|&(e, _)| straight_origin(&part.edges[e].curve).is_none())
            && let Some(cost) =
                super::corresponding::loop_match_error((part, face), loop_index, fillers)
        {
            straight_points = None;
            error += cost;
            continue;
        }
        let mut points = Vec::new();
        let mut first = None;
        let mut last = None;
        for &(index, forward) in &lp.edges {
            if part.unresolved_edges().contains(&index) {
                return None;
            }
            let edge = &part.edges[index];
            let (a, b) = if forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            if last.is_some_and(|p| p != a) {
                return None;
            }
            first.get_or_insert(a);
            last = Some(b);
            let vertices = [a, b].map(|p| project(frame, p));
            use super::representations::{Representation, bounded_alternatives};
            let representations = bounded_alternatives((part, face), index)?;
            if representations
                .iter()
                .any(|r| matches!(r, Representation::Hull(_)))
            {
                straight_points = None;
                let patch = clip?;
                error += curve_error(edge, frame, patch)?;
                for r in representations {
                    match r {
                        Representation::Exact(curve) => {
                            error += curve_error(
                                &Edge {
                                    curve,
                                    ..edge.clone()
                                },
                                frame,
                                patch,
                            )?;
                        }
                        Representation::Hull(hull) => {
                            let mut bounds = Box2::empty();
                            let mut size = base_size;
                            for p in [a, b].into_iter().chain(hull) {
                                size = p.into_iter().map(f64::abs).fold(size, f64::max);
                                bounds.add(project(frame, p))?;
                            }
                            error += bounds.intersection_area(patch, 256.0 * f64::EPSILON * size);
                        }
                    }
                }
                points.extend(vertices);
                continue;
            }
            let alternatives: Vec<_> = representations
                .into_iter()
                .filter_map(|r| match r {
                    Representation::Exact(curve) => Some(curve),
                    Representation::Hull(_) => None,
                })
                .collect();
            if straight_origin(&edge.curve).is_none() {
                straight_points = None;
                let patch = clip?;
                error += curve_error(edge, frame, patch)?;
                for curve in alternatives {
                    let alternate = Edge {
                        curve,
                        ..edge.clone()
                    };
                    error += curve_error(&alternate, frame, patch)?;
                }
                points.extend(vertices);
                continue;
            }
            let origin = straight_origin(&edge.curve)?;
            let mut size = origin.into_iter().map(f64::abs).fold(base_size, f64::max);
            for curve in &alternatives {
                let origin = straight_origin(curve)?;
                size = origin.iter().copied().map(f64::abs).fold(size, f64::max);
            }
            let pad = 64.0 * f64::EPSILON * size;
            let world = [
                straight_foot(&edge.curve, a)?,
                straight_foot(&edge.curve, b)?,
            ];
            let feet = world.map(|p| project(frame, p));
            outer_size = outer_size.max(size);
            if let Some(points) = &mut straight_points {
                points.extend(vertices);
                points.extend(feet);
            }
            let mut radii = std::array::from_fn::<_, 2, _>(|i| {
                (vertices[i].0 - feet[i].0).hypot(vertices[i].1 - feet[i].1) + pad
            });
            let uncertainty = part
                .file_uncertainty()
                .filter(|x| x.is_finite() && *x >= 0.0);
            let limit = uncertainty.map_or(f64::INFINITY, |x| x + pad);
            if [a, b]
                .into_iter()
                .zip(world)
                .any(|(p, q)| geom::dist(p, q) > limit)
            {
                return None;
            }
            if uncertainty.is_none() {
                // Constructed parts need not declare a modelling precision. Retain both
                // the vertex path and finite analytic segment in their area bracket; an
                // absent declaration is neither zero uncertainty nor a licence to heal.
                error += paired_error(vertices, feet, size, clip)?;
                for curve in &alternatives {
                    let projected = [straight_foot(curve, a)?, straight_foot(curve, b)?]
                        .map(|p| project(frame, p));
                    if let Some(points) = &mut straight_points {
                        points.extend(projected);
                    }
                    error += paired_error(vertices, projected, size, clip)?;
                }
                points.extend(vertices);
                continue;
            }
            error += segment_error(feet, [0.0; 2], size, clip)?;
            for curve in &alternatives {
                let param_feet = [a, b].map(|p| {
                    let q = straight_foot(curve, p)?;
                    (geom::dist(p, q) <= limit).then_some(project(frame, q))
                });
                let [Some(c), Some(d)] = param_feet else {
                    return None;
                };
                let alternate = [c, d];
                if let Some(points) = &mut straight_points {
                    points.extend(alternate);
                }
                for i in 0..2 {
                    radii[i] = radii[i].max(
                        (vertices[i].0 - alternate[i].0).hypot(vertices[i].1 - alternate[i].1)
                            + pad,
                    );
                }
                error += paired_error(feet, alternate, size, clip)?;
            }
            for (vertex, radius) in vertices.into_iter().zip(radii) {
                let local = Box2 {
                    low: vertex,
                    high: vertex,
                };
                error += clip.map_or(4.0 * radius * radius, |patch| {
                    local.intersection_area(patch, radius)
                });
            }
            // These joins stay inside the endpoint boxes charged above. They describe no
            // additional support outside those boxes and do not modify the input B-Rep.
            points.extend([vertices[0], feet[0], feet[1], vertices[1]]);
        }
        if last != first {
            return None;
        }
        points.push(points[0]);
        contours.push(points);
    }
    Some(Boundary {
        loops: contours,
        line_error: error,
        straight_points,
        size: outer_size,
    })
}

fn project(frame: &Frame, p: V3) -> Point {
    let delta = geom::sub(p, frame.origin);
    (geom::dot(delta, frame.x), geom::dot(delta, frame.y))
}

#[cfg(test)]
fn line_error(edge: &Edge, frame: &Frame, patch: Option<Box2>) -> Option<f64> {
    let Curve::Line { origin, dir } = edge.curve else {
        return None;
    };
    if !geom::norm(dir).is_finite() || geom::norm(dir) == 0.0 {
        return None;
    }
    let mut size = frame.origin.into_iter().map(f64::abs).fold(1.0, f64::max);
    let mut errors = [0.0_f64; 2];
    for p in [origin, edge.start, edge.end] {
        size = p.into_iter().map(f64::abs).fold(size, f64::max);
    }
    let ends = [edge.start, edge.end].map(|p| project(frame, p));
    for (i, (p, uv)) in [edge.start, edge.end].into_iter().zip(ends).enumerate() {
        let q = project(frame, edge.curve.value(edge.curve.parameter(p)));
        if !q.0.is_finite() || !q.1.is_finite() {
            return None;
        }
        errors[i] = (q.0 - uv.0).hypot(q.1 - uv.1);
    }
    segment_error(ends, errors, size, patch)
}

fn segment_error(
    ends: [Point; 2],
    errors: [f64; 2],
    size: f64,
    patch: Option<Box2>,
) -> Option<f64> {
    let mut bounds = Box2::empty();
    for point in ends {
        bounds.add(point)?;
    }
    // The swept segment lies inside this capsule and its expanded axis-aligned box.
    // Neither bound enlarges the support: an indecisive area bracket falls back.
    let radius = errors[0].max(errors[1]) + 64.0 * f64::EPSILON * size;
    let length = (ends[1].0 - ends[0].0).hypot(ends[1].1 - ends[0].1);
    let capsule = 2.0 * radius * length + std::f64::consts::PI * radius * radius;
    // A target boundary can move outside its nominal box, so only support error is clipped.
    let Some(patch) = patch else {
        return Some(capsule);
    };
    let full = capsule.min(bounds.intersection_area(patch, radius));
    // The displacement between two straight segments is affine. On each half its norm is
    // bounded by the larger interpolated endpoint bound. Their two capsules therefore still
    // contain the entire boundary sweep, without charging the largest error along all of it.
    // The extra eight epsilons cover midpoint interpolation; target uncertainty stays full.
    let middle = (
        ends[0].0 + 0.5 * (ends[1].0 - ends[0].0),
        ends[0].1 + 0.5 * (ends[1].1 - ends[0].1),
    );
    let points = [ends[0], middle, ends[1]];
    let deviations = [errors[0], 0.5 * errors[0] + 0.5 * errors[1], errors[1]];
    let mut split = 0.0;
    for i in 0..2 {
        let radius = deviations[i].max(deviations[i + 1]) + 72.0 * f64::EPSILON * size;
        let mut bounds = Box2::empty();
        bounds.add(points[i])?;
        bounds.add(points[i + 1])?;
        let length = (points[i + 1].0 - points[i].0).hypot(points[i + 1].1 - points[i].1);
        let capsule = 2.0 * radius * length + std::f64::consts::PI * radius * radius;
        split += capsule.min(bounds.intersection_area(patch, radius));
    }
    Some(full.min(split))
}

/// A straight homotopy F(t,s)=p(t)+s(q(t)-p(t)) has an affine Jacobian determinant.
/// Convexity of absolute value bounds its integral by the mean of the four corner values.
/// The area formula bounds the image even when the sweep folds or overlaps itself.
fn paired_error(
    ends: [Point; 2],
    other: [Point; 2],
    size: f64,
    patch: Option<Box2>,
) -> Option<f64> {
    let delta = |a: Point, b: Point| (a.0 - b.0, a.1 - b.1);
    let norm = |p: Point| p.0.hypot(p.1);
    let directions = [delta(ends[1], ends[0]), delta(other[1], other[0])];
    let deviations = [delta(other[0], ends[0]), delta(other[1], ends[1])];
    let mut area = 0.0;
    for a in directions {
        for b in deviations {
            area += (a.0 * b.1 - a.1 * b.0).abs() * 0.25;
        }
    }
    let pad = 64.0 * f64::EPSILON * size;
    area += 8.0
        * pad
        * (directions.into_iter().map(norm).sum::<f64>()
            + deviations.into_iter().map(norm).sum::<f64>())
        + 32.0 * pad * pad;
    if !area.is_finite() {
        return None;
    }
    Some(segment_error(ends, deviations.map(norm), size, patch)?.min(area))
}

fn curve_error(edge: &Edge, frame: &Frame, patch: Box2) -> Option<f64> {
    let mut bounds = Box2::empty();
    let mut size = frame.origin.into_iter().map(f64::abs).fold(1.0, f64::max);
    match &edge.curve {
        Curve::Circle { frame, radius } => {
            size = frame.origin.into_iter().map(f64::abs).fold(size, f64::max);
            size = size.max(radius.abs());
        }
        Curve::Ellipse {
            frame,
            major,
            minor,
        } => {
            size = frame.origin.into_iter().map(f64::abs).fold(size, f64::max);
            size = size.max(major.abs()).max(minor.abs());
        }
        _ => {}
    }
    let mut add = |p: V3| {
        size = p.into_iter().map(f64::abs).fold(size, f64::max);
        bounds.add(project(frame, p))
    };
    for p in [edge.start, edge.end] {
        add(p)?;
    }
    match &edge.curve {
        Curve::Line { .. } => {
            for p in [edge.start, edge.end] {
                add(edge.curve.value(edge.curve.parameter(p)))?;
            }
        }
        Curve::Nurbs(curve) => {
            if !curve.weights().iter().all(|w| w.is_finite() && *w > 0.0) {
                return None;
            }
            for &p in curve.control_points() {
                add(p)?;
            }
            if bounds.disjoint(patch, 64.0 * f64::EPSILON * size) {
                return Some(0.0);
            }
            if let Some(mut trimmed) = finite_spline_bounds(edge, frame, patch) {
                for p in [edge.start, edge.end] {
                    trimmed.add(project(frame, p))?;
                }
                return Some(trimmed.intersection_area(patch, 64.0 * f64::EPSILON * size));
            }
        }
        Curve::Circle { .. } | Curve::Ellipse { .. } => {
            let interval = edge_interval(
                &edge.curve,
                edge.start,
                edge.end,
                edge.same_sense,
                edge.vertices.0 == edge.vertices.1,
            );
            if !interval.0.is_finite() || !interval.1.is_finite() {
                return None;
            }
            add(edge.curve.value(interval.0))?;
            add(edge.curve.value(interval.1))?;
            for p in extremes_along(&edge.curve, interval, &[frame.x, frame.y]) {
                add(p)?;
            }
        }
    }
    // Projection roundoff enlarges the error bound, never the support or acceptance budget.
    Some(bounds.intersection_area(patch, 64.0 * f64::EPSILON * size))
}

/// Refined positive Bezier hulls enclose the finite spline, including periodic wraps.
/// Every final hull is retained even when refinement exhausts its work budget.
fn finite_spline_bounds(edge: &Edge, frame: &Frame, patch: Box2) -> Option<Box2> {
    let Curve::Nurbs(curve) = &edge.curve else {
        return None;
    };
    let interval = edge_interval(
        &edge.curve,
        edge.start,
        edge.end,
        edge.same_sense,
        edge.vertices.0 == edge.vertices.1,
    );
    let accuracy = (patch.high.0 - patch.low.0).min(patch.high.1 - patch.low.1) * 1e-13;
    let (low, high) =
        super::spline_bounds::bounds(curve, interval, |p| project(frame, p), accuracy)?;
    Some(Box2 { low, high })
}

fn uncovered_area(mine: &Loops, theirs: &[Loops], bounds: Box2) -> Option<f64> {
    let loops = std::iter::once(mine).chain(theirs).flatten();
    let mut heights = vec![bounds.low.1, bounds.high.1];
    let mut segments = Vec::new();
    for lp in loops {
        for w in lp.windows(2) {
            let mut box2 = Box2::empty();
            box2.add(w[0])?;
            box2.add(w[1])?;
            for p in w {
                if bounds.low.1 < p.1 && p.1 < bounds.high.1 {
                    heights.push(p.1);
                }
            }
            if !box2.disjoint(bounds, 0.0) && w[0].1 != w[1].1 {
                segments.push((w[0], w[1]));
                if segments.len() > MAX_SEGMENTS {
                    return None;
                }
            }
        }
    }
    for (i, &a) in segments.iter().enumerate() {
        for &b in &segments[..i] {
            if let Some(v) = crossing_height(a, b)
                && bounds.low.1 < v
                && v < bounds.high.1
            {
                heights.push(v);
            }
        }
    }
    heights.sort_by(f64::total_cmp);
    heights.dedup();
    let mut area = 0.0;
    for band in heights.windows(2) {
        let v = band[0] + 0.5 * (band[1] - band[0]);
        let mut left = intervals(mine, v);
        for support in theirs {
            for cut in intervals(support, v) {
                left = subtract(&left, cut);
            }
        }
        let length = left.iter().map(|(a, b)| (b - a).max(0.0)).sum::<f64>();
        area += (band[1] - band[0]) * length;
    }
    area.is_finite().then_some(area)
}

pub(super) fn polygon_area(loops: &Vec<Vec<(f64, f64)>>) -> Option<f64> {
    let mut bounds = Box2::empty();
    for &p in loops.iter().flatten() {
        bounds.add(p)?;
    }
    uncovered_area(loops, &[], bounds)
}

fn crossing_height((a, b): Segment, (c, d): Segment) -> Option<f64> {
    let cross = |p: Point, q: Point| p.0 * q.1 - p.1 * q.0;
    let delta = |p: Point, q: Point| (p.0 - q.0, p.1 - q.1);
    let (r, s, t) = (delta(b, a), delta(d, c), delta(c, a));
    let determinant = cross(r, s);
    if determinant == 0.0 {
        return None;
    }
    let (p, q) = (cross(t, s) / determinant, cross(t, r) / determinant);
    ((0.0..=1.0).contains(&p) && (0.0..=1.0).contains(&q)).then_some(a.1 + p * r.1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::{Edge, Face, Loop, Part};

    fn face(contours: &[Vec<Point>]) -> Part {
        let mut edges = Vec::new();
        let mut loops = Vec::new();
        for points in contours {
            let base = edges.len();
            for (i, &p) in points.iter().enumerate() {
                let q = points[(i + 1) % points.len()];
                let (start, end) = ([p.0, p.1, 0.0], [q.0, q.1, 0.0]);
                edges.push(Edge {
                    curve: Curve::Line {
                        origin: start,
                        dir: geom::unit(geom::sub(end, start)).unwrap(),
                    },
                    start,
                    end,
                    vertices: (base + i, base + (i + 1) % points.len()),
                    same_sense: true,
                    samples: vec![start, end],
                });
            }
            loops.push(Loop {
                edges: (base..edges.len()).map(|i| (i, true)).collect(),
                vertex: None,
            });
        }
        let a = contours[0][0];
        let b = contours[0][1];
        let x = geom::unit([b.0 - a.0, b.1 - a.1, 0.0]).unwrap();
        Part::new(
            vec![Face {
                surface: Surface::Plane {
                    frame: Frame {
                        origin: [a.0, a.1, 0.0],
                        x,
                        y: [-x[1], x[0], 0.0],
                        z: [0.0, 0.0, 1.0],
                    },
                },
                reversed: false,
                loops,
                solid: None,
                pcurves: Vec::new(),
            }],
            edges,
            Vec::new(),
        )
    }

    fn rectangle(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point> {
        vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    }

    fn rotate(points: &[Point], angle: f64) -> Vec<Point> {
        let (s, c) = angle.sin_cos();
        points
            .iter()
            .map(|&(x, y)| (c * x - s * y, s * x + c * y))
            .collect()
    }

    #[test]
    fn a_separate_finite_disk_fills_a_support_loop_but_retains_annular_gaps() {
        use crate::cover::{Coverage, assess};
        let patch = face(&[rectangle(-2.0, -2.0, 2.0, 2.0)]);
        let disk = |radius: f64, offset: f64| {
            let mut p = face(&[rectangle(-2.0, -2.0, 2.0, 2.0)]);
            let curve = Curve::Circle {
                frame: Frame {
                    origin: [offset, 0.0, 0.0],
                    x: [1.0, 0.0, 0.0],
                    y: [0.0, 1.0, 0.0],
                    z: [0.0, 0.0, 1.0],
                },
                radius,
            };
            let seam = curve.value(0.0);
            p.edges = vec![Edge {
                curve,
                start: seam,
                end: seam,
                vertices: (0, 0),
                same_sense: true,
                samples: vec![],
            }];
            p.faces[0].loops = vec![Loop {
                edges: vec![(0, true)],
                vertex: None,
            }];
            p
        };
        let filler = disk(1.0, 0.0);
        let mut wall = face(&[rectangle(-2.0, -2.0, 2.0, 2.0)]);
        let hole_index = wall.edges.len();
        wall.edges.push(filler.edges[0].clone());
        wall.faces[0].loops.push(Loop {
            edges: vec![(hole_index, true)],
            vertex: None,
        });
        assert_ne!(assess((&patch, 0), &[(&wall, 0)]), Coverage::Covered);
        assert_eq!(
            assess((&patch, 0), &[(&wall, 0), (&filler, 0)]),
            Coverage::Covered
        );
        // The analytic annular gap is pi*(1-r^2), far above the original area budget.
        for bad in [disk(0.999, 0.0), disk(1.0, 0.001)] {
            assert_ne!(
                assess((&patch, 0), &[(&wall, 0), (&bad, 0)]),
                Coverage::Covered
            );
        }
        let mut perforated = disk(1.0, 0.0);
        perforated.edges.push(disk(0.1, 0.0).edges.remove(0));
        perforated.faces[0].loops.push(Loop {
            edges: vec![(1, true)],
            vertex: None,
        });
        assert_ne!(
            assess((&patch, 0), &[(&wall, 0), (&perforated, 0)]),
            Coverage::Covered
        );
    }

    #[test]
    fn a_remote_native_surface_pcurve_keeps_its_complete_control_hull() {
        use crate::brep::RawCurve2d;
        use crate::nurbs::{NurbsCurve, NurbsSurface};
        let patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        for near in [false, true] {
            let mut support = face(&[rectangle(-1.0, -1.0, 20.0, 20.0)]);
            let curve = Curve::Nurbs(
                NurbsCurve::new(
                    2,
                    vec![
                        [10.0, 10.0, 0.0],
                        [10.0, 11.0, 0.0],
                        [11.0, 10.0, 0.0],
                        [10.0, 10.0, 0.0],
                    ],
                    vec![1.0; 4],
                    vec![0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0],
                )
                .unwrap(),
            );
            support.edges.push(Edge {
                curve,
                start: [10.0, 10.0, 0.0],
                end: [10.0, 10.0, 0.0],
                vertices: (4, 4),
                same_sense: true,
                samples: vec![],
            });
            support.faces[0].loops.push(Loop {
                edges: vec![(4, false)],
                vertex: None,
            });
            let lower = if near { 0.0 } else { 9.0 };
            support.faces.push(Face {
                surface: Surface::Freeform {
                    kind: "BSPLINE",
                    surface: Box::new(
                        NurbsSurface::new(
                            1,
                            1,
                            vec![
                                vec![[lower, lower, 0.0], [lower, 11.0, 0.0]],
                                vec![[11.0, lower, 0.0], [11.0, 11.0, 0.0]],
                            ],
                            vec![vec![1.0; 2]; 2],
                            vec![0.0, 0.0, 1.0, 1.0],
                            vec![0.0, 0.0, 1.0, 1.0],
                        )
                        .unwrap(),
                    ),
                },
                reversed: false,
                solid: None,
                pcurves: vec![],
                loops: vec![Loop {
                    edges: vec![(4, true)],
                    vertex: None,
                }],
            });
            let support = Part::new(support.faces.clone(), support.edges.clone(), vec![])
                .with_file_pcurves(
                    std::collections::BTreeMap::from([(
                        1,
                        vec![(
                            4,
                            RawCurve2d::BSpline {
                                degree: 2,
                                points: vec![[0.5, 0.5], [0.5, 1.0 + 1e-8], [1.0, 0.5], [0.5, 0.5]],
                                multiplicities: vec![3, 1, 3],
                                knots: vec![0.0, 1.0, 2.0],
                                weights: None,
                            },
                        )],
                    )]),
                    std::collections::BTreeMap::from([(1, vec![0])]),
                    (1.0, 1.0),
                    Some(1e-7),
                );
            assert_eq!(super::super::covered((&patch, 0), &[(&support, 0)]), !near);
        }
    }

    #[test]
    fn finite_straight_sweeps_keep_their_transverse_area() {
        let ends = [(0.0, 0.0), (1.0, 0.0)];
        for (other, expected) in [
            ([(0.0, 0.1), (1.0, 0.1)], 0.1),
            ([(0.0, 0.0), (1.0, 0.1)], 0.05),
        ] {
            let bound = paired_error(ends, other, 1.0, None).unwrap();
            assert!(bound >= expected && bound <= expected + 1e-10);
        }
        // A longitudinal slide has no swept area, although its Euclidean displacement is large.
        assert!(paired_error(ends, [(0.1, 0.0), (1.1, 0.0)], 2.0, None).unwrap() < 1e-10);
    }

    #[test]
    fn straight_strip_gaps_do_not_depend_on_the_patch_frame() {
        let square = rectangle(0.0, 0.0, 1.0, 1.0);
        for degrees in [0.0f64, 17.0, 45.0, 89.0, 123.0, 179.0, 271.0] {
            for gap in [
                0.0, 1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 3.99e-4, 4.01e-4, 0.1,
            ] {
                let support = face(&[rotate(
                    &rectangle(0.0, 0.0, 1.0, 1.0 - gap),
                    degrees.to_radians(),
                )]);
                let mut points = rotate(&square, degrees.to_radians());
                for _ in 0..4 {
                    let patch = face(&[points.clone()]);
                    assert_eq!(
                        covered((&patch, 0), &[(&support, 0)]),
                        Some(gap <= 1e-9),
                        "angle={degrees}, gap={gap}"
                    );
                    assert_eq!(
                        super::super::covered((&patch, 0), &[(&support, 0)]),
                        gap <= 1e-9
                    );
                    points.rotate_left(1);
                }
            }
        }
    }

    #[test]
    fn holes_are_not_filled_and_separate_supports_can_fill_them() {
        let outer = rectangle(0.0, 0.0, 1.0, 1.0);
        let hole = rectangle(0.4, 0.4, 0.6, 0.6);
        let patch = face(std::slice::from_ref(&outer));
        let ring = face(&[outer, hole.clone()]);
        let filler = face(&[hole]);
        assert_eq!(covered((&patch, 0), &[(&ring, 0)]), Some(false));
        assert_eq!(
            covered((&patch, 0), &[(&ring, 0), (&filler, 0)]),
            Some(true)
        );
        assert_eq!(
            covered((&patch, 0), &[(&filler, 0), (&ring, 0)]),
            Some(true)
        );
    }

    #[test]
    fn crossing_boundaries_are_integration_events() {
        let mut patch = rectangle(0.0, 0.0, 1.0, 1.0);
        patch.push(patch[0]);
        let left = vec![(0.0, 0.0), (0.2, 0.0), (0.8, 1.0), (0.0, 1.0), (0.0, 0.0)];
        let right = vec![(0.8, 0.0), (1.0, 0.0), (1.0, 1.0), (0.2, 1.0), (0.8, 0.0)];
        for degrees in [0.0f64, 17.0, 45.0, 89.0, 123.0] {
            let mine = vec![rotate(&patch, degrees.to_radians())];
            let theirs = vec![
                vec![rotate(&left, degrees.to_radians())],
                vec![rotate(&right, degrees.to_radians())],
            ];
            let mut bounds = Box2::empty();
            for &p in mine.iter().flatten() {
                bounds.add(p).unwrap();
            }
            let area = uncovered_area(&mine, &theirs, bounds).unwrap();
            assert!((area - 0.15).abs() < 1e-12, "angle={degrees}, area={area}");
        }
    }

    #[test]
    fn touching_diagonal_supports_cover_without_a_scan_margin() {
        let patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        let a = face(&[vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]]);
        let b = face(&[vec![(1.0, 1.0), (0.0, 1.0), (1.0, 0.0)]]);
        assert_eq!(covered((&patch, 0), &[(&a, 0), (&b, 0)]), Some(true));
    }

    #[test]
    fn analytic_curve_boxes_do_not_discard_a_curve_near_the_patch() {
        let frame = Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        };
        let patch = Box2 {
            low: (0.0, 0.0),
            high: (1.0, 1.0),
        };
        for centre in [0.5, 10.0] {
            let curve = Curve::Circle {
                frame: Frame {
                    origin: [centre, centre, 0.0],
                    ..frame
                },
                radius: 0.1,
            };
            let point = [centre + 0.1, centre, 0.0];
            let edge = Edge {
                curve,
                start: point,
                end: point,
                vertices: (0, 0),
                same_sense: true,
                samples: vec![],
            };
            let error = curve_error(&edge, &frame, patch).unwrap();
            if centre == 10.0 {
                assert_eq!(error, 0.0);
            } else {
                assert!((error - 0.04).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn discontinuous_loops_are_not_closed_by_guessing() {
        let patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        let mut support = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        support.edges[0].start[0] = 0.01;
        assert_eq!(covered((&patch, 0), &[(&support, 0)]), None);
    }

    #[test]
    fn a_trimmed_arc_does_not_use_the_whole_circles_box() {
        let frame = Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        };
        let edge = Edge {
            curve: Curve::Circle {
                frame: Frame {
                    origin: [0.5, 2.0, 0.0],
                    ..frame
                },
                radius: 1.5,
            },
            start: [2.0, 2.0, 0.0],
            end: [-1.0, 2.0, 0.0],
            vertices: (0, 1),
            same_sense: true,
            samples: vec![],
        };
        let patch = Box2 {
            low: (0.0, 0.0),
            high: (1.0, 1.0),
        };
        assert_eq!(curve_error(&edge, &frame, patch), Some(0.0));
        let opposite = Edge {
            same_sense: false,
            ..edge
        };
        assert!(curve_error(&opposite, &frame, patch).unwrap() > 0.49);
    }

    #[test]
    fn a_tangent_circle_is_charged_to_the_original_area_budget() {
        let mut support = face(&[rectangle(-4.0, -4.0, 4.0, 4.0)]);
        let index = support.edges.len();
        support.edges.push(Edge {
            curve: Curve::Circle {
                frame: Frame {
                    origin: [-1.0, 0.5, 0.0],
                    x: [1.0, 0.0, 0.0],
                    y: [0.0, 1.0, 0.0],
                    z: [0.0, 0.0, 1.0],
                },
                radius: 1.0,
            },
            start: [0.0, 0.5, 0.0],
            end: [0.0, 0.5, 0.0],
            vertices: (index, index),
            same_sense: true,
            samples: vec![],
        });
        support.faces[0].loops.push(Loop {
            edges: vec![(index, true)],
            vertex: None,
        });
        let support = Part::new(support.faces.clone(), support.edges.clone(), vec![]);
        let mut square = rectangle(0.0, 0.0, 1.0, 1.0);
        for _ in 0..4 {
            let patch = face(&[square.clone()]);
            assert_eq!(covered((&patch, 0), &[(&support, 0)]), Some(true));
            square.rotate_left(1);
        }
        let patch = face(&[rectangle(-0.1, 0.0, 1.0, 1.0)]);
        assert_eq!(
            covered((&patch, 0), &[(&support, 0)]),
            None,
            "a material circle segment cannot be discarded as roundoff"
        );
    }

    #[test]
    fn declared_strip_gaps_remain_gaps_under_unit_scaling() {
        for scale in [0.001, 1.0, 1000.0] {
            let patch = face(&[rectangle(0.0, 0.0, scale, scale)]);
            let support = face(&[rectangle(0.0, 0.0, scale, scale * (1.0 - 1e-4))]);
            assert_eq!(covered((&patch, 0), &[(&support, 0)]), Some(false));
        }
    }

    #[test]
    fn rational_curves_use_the_positive_weight_control_hull() {
        let frame = Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        };
        let patch = Box2 {
            low: (0.0, 0.0),
            high: (1.0, 1.0),
        };
        for offset in [0.0, 10.0] {
            let points = vec![
                [offset, offset, 0.0],
                [offset + 1.0, offset + 2.0, 0.0],
                [offset + 2.0, offset, 0.0],
            ];
            let curve = crate::nurbs::NurbsCurve::new(
                2,
                points.clone(),
                vec![1.0, 3.0, 1.0],
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            )
            .unwrap();
            let edge = Edge {
                curve: Curve::Nurbs(curve),
                start: points[0],
                end: points[2],
                vertices: (0, 1),
                same_sense: true,
                samples: vec![],
            };
            let error = curve_error(&edge, &frame, patch).unwrap();
            if offset == 10.0 {
                assert_eq!(error, 0.0);
            } else {
                assert_eq!(error, 1.0);
            }
        }
    }

    #[test]
    fn off_line_vertices_do_not_establish_a_tiny_gap() {
        let mut square = rectangle(0.0, 0.0, 1.0, 1.0);
        for gap in [1e-8, 1e-4] {
            let mut support = face(&[rectangle(0.0, 0.0, 1.0, 1.0 - gap)]);
            let Curve::Line { origin, .. } = &mut support.edges[2].curve else {
                unreachable!()
            };
            origin[1] += 2e-8;
            for _ in 0..4 {
                let patch = face(&[square.clone()]);
                assert_eq!(
                    covered((&patch, 0), &[(&support, 0)]),
                    if gap == 1e-8 { None } else { Some(false) }
                );
                square.rotate_left(1);
            }
        }
    }

    fn outward_alternate(gap: f64, angle: f64) -> Part {
        use crate::brep::RawCurve2d;
        let support = face(&[rotate(&rectangle(gap, -1.0, 1.0, 2.0), angle)]);
        support.with_file_pcurves(
            std::collections::BTreeMap::from([(
                0,
                vec![(
                    1,
                    RawCurve2d::Line {
                        point: [2.0 - gap, 0.0],
                        direction: [0.0, 1.0],
                        magnitude: 1.0,
                    },
                )],
            )]),
            Default::default(),
            (1.0, 1.0),
            None,
        )
    }

    #[test]
    fn an_outward_alternate_cannot_hide_a_gap_on_the_opposite_side() {
        use crate::cover::{Coverage, assess};
        for angle in [0.0, 0.37, 1.2] {
            let support = outward_alternate(5e-9, angle);
            let mut square = rotate(&rectangle(0.0, 0.0, 1.0, 1.0), angle);
            for _ in 0..4 {
                let patch = face(&[square.clone()]);
                let bracket = finite_area((&patch, 0), &[(&support, 0)], &[]).unwrap();
                assert!(bracket.gap <= bracket.tolerance + bracket.error);
                assert_eq!(assess((&patch, 0), &[(&support, 0)]), Coverage::Uncovered);
                square.rotate_left(1);
            }
        }
    }

    #[test]
    fn outer_envelopes_retain_other_supports_and_target_uncertainty() {
        use crate::cover::{Coverage, assess};
        let gap = 5e-9;
        let support = outward_alternate(gap, 0.0);
        let patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        let filler = face(&[rectangle(-1.0, -1.0, 2.0 * gap, 2.0)]);
        assert_ne!(
            assess((&patch, 0), &[(&support, 0), (&filler, 0)]),
            Coverage::Uncovered
        );
        let complete = face(&[rectangle(-1.0, -1.0, 2.0, 2.0)]);
        assert_eq!(
            assess((&patch, 0), &[(&support, 0), (&complete, 0)]),
            Coverage::Covered
        );
        let mut uncertain = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        let Curve::Line { origin, .. } = &mut uncertain.edges[3].curve else {
            unreachable!()
        };
        origin[0] += 2.0 * gap;
        assert_eq!(
            assess((&uncertain, 0), &[(&support, 0)]),
            Coverage::Indeterminate
        );
    }

    #[test]
    fn an_alternate_that_can_cover_the_gap_stays_indeterminate() {
        use crate::brep::RawCurve2d;
        use crate::cover::{Coverage, assess};
        let gap = 5e-9;
        let base = outward_alternate(gap, 0.0);
        let support = base.with_file_pcurves(
            std::collections::BTreeMap::from([(
                0,
                vec![(
                    3,
                    RawCurve2d::Line {
                        point: [-gap, 0.0],
                        direction: [0.0, -1.0],
                        magnitude: 1.0,
                    },
                )],
            )]),
            Default::default(),
            (1.0, 1.0),
            None,
        );
        let patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        assert_eq!(
            assess((&patch, 0), &[(&support, 0)]),
            Coverage::Indeterminate
        );
    }

    #[test]
    fn target_line_discrepancies_are_also_uncertainty() {
        let support = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        let mut patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        let Curve::Line { origin, .. } = &mut patch.edges[2].curve else {
            unreachable!()
        };
        origin[1] += 1e-7;
        assert_eq!(covered((&patch, 0), &[(&support, 0)]), None);
    }

    #[test]
    fn target_uncertainty_includes_the_outside_of_its_nominal_box() {
        let patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        let Surface::Plane { frame } = &patch.faces[0].surface else {
            unreachable!()
        };
        let mut edge = patch.edges[2].clone();
        let Curve::Line { origin, .. } = &mut edge.curve else {
            unreachable!()
        };
        origin[1] += 0.1;
        let full = line_error(&edge, frame, None).unwrap();
        let clipped = line_error(
            &edge,
            frame,
            Some(Box2 {
                low: (0.0, 0.0),
                high: (1.0, 1.0),
            }),
        )
        .unwrap();
        let capsule = 0.2 + std::f64::consts::PI * 0.01;
        assert!((full - capsule).abs() < 1e-12);
        assert!(clipped < full - 0.1);
    }

    #[test]
    fn endpoint_varying_error_can_resolve_a_near_budget_gap() {
        let mut support = face(&[rectangle(0.0, 0.0, 1.0, 1.0 - 1.036e-9)]);
        support.edges[2].curve = Curve::Line {
            origin: [0.0, 1.0 - 1.036e-9 + 0.8e-11, 0.0],
            dir: geom::unit([1.0, 1.2e-11, 0.0]).unwrap(),
        };
        let mut square = rectangle(0.0, 0.0, 1.0, 1.0);
        for _ in 0..4 {
            let patch = face(&[square.clone()]);
            assert_eq!(covered((&patch, 0), &[(&support, 0)]), Some(false));
            square.rotate_left(1);
        }
    }

    #[test]
    fn an_indeterminate_boundary_does_not_reenter_sampled_acceptance() {
        use super::super::{Coverage, assess, covered as established, sampled_legacy};
        let mut support = face(&[rectangle(0.0, 0.0, 1.0, 1.0 - 1e-8)]);
        let Curve::Line { origin, .. } = &mut support.edges[2].curve else {
            unreachable!()
        };
        origin[1] += 2e-8;
        let mut square = rectangle(0.0, 0.0, 1.0, 1.0);
        let mut old = Vec::new();
        for _ in 0..4 {
            let patch = face(&[square.clone()]);
            assert_eq!(
                assess((&patch, 0), &[(&support, 0)]),
                Coverage::Indeterminate
            );
            assert!(!established((&patch, 0), &[(&support, 0)]));
            old.push(sampled_legacy((&patch, 0), &[(&support, 0)]));
            square.rotate_left(1);
        }
        assert!(old.contains(&true) && old.contains(&false));
    }

    #[test]
    fn the_patch_parameter_origin_and_axes_do_not_select_the_answer() {
        use super::super::{Coverage, assess};
        let supports = [
            face(&[rectangle(-1.0, -1.0, 2.0, 2.0)]),
            face(&[rectangle(0.0, 0.0, 1.0, 1.0 - 1e-4)]),
        ];
        for angle in [0.0_f64, 0.37, 1.2] {
            for (support, expected) in supports
                .iter()
                .zip([Coverage::Covered, Coverage::Uncovered])
            {
                let mut patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
                let Surface::Plane { frame } = &mut patch.faces[0].surface else {
                    unreachable!()
                };
                frame.origin = [100.0, -70.0, 0.0];
                frame.x = [angle.cos(), angle.sin(), 0.0];
                frame.y = [-angle.sin(), angle.cos(), 0.0];
                assert_eq!(assess((&patch, 0), &[(support, 0)]), expected);
            }
        }
    }

    #[test]
    fn a_far_3d_curve_does_not_hide_an_interior_pcurve() {
        use crate::brep::RawCurve2d;
        let patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        let mut support = face(&[rectangle(-2.0, -2.0, 2.0, 2.0)]);
        let index = support.edges.len();
        let vertex = [10.25, 10.0, 0.0];
        support.edges.push(Edge {
            curve: Curve::Circle {
                frame: Frame {
                    origin: [10.0, 10.0, 0.0],
                    x: [1.0, 0.0, 0.0],
                    y: [0.0, 1.0, 0.0],
                    z: [0.0, 0.0, 1.0],
                },
                radius: 0.25,
            },
            start: vertex,
            end: vertex,
            vertices: (index, index),
            same_sense: true,
            samples: vec![],
        });
        support.faces[0].loops.push(Loop {
            edges: vec![(index, true)],
            vertex: None,
        });
        for pcurve in [
            RawCurve2d::Circle {
                point: [2.5, 2.5],
                reference: None,
                radius: 0.25,
            },
            RawCurve2d::Line {
                point: [2.5, 2.5],
                direction: [1.0, -1.0],
                magnitude: 1.0,
            },
        ] {
            let bounded =
                crate::brep::Part::new(support.faces.clone(), support.edges.clone(), vec![])
                    .with_file_pcurves(
                        std::collections::BTreeMap::from([(0, vec![(index, pcurve)])]),
                        Default::default(),
                        (1.0, 1.0),
                        Some(1e-7),
                    );
            assert_eq!(covered((&patch, 0), &[(&bounded, 0)]), None);
            assert_ne!(
                super::super::assess((&patch, 0), &[(&bounded, 0)]),
                super::super::Coverage::Covered
            );
        }
    }

    #[test]
    fn a_finite_spline_box_does_not_use_other_parameter_spans() {
        let frame = Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        };
        let curve = Curve::Nurbs(
            crate::nurbs::NurbsCurve::new(
                2,
                vec![[0.0, 10.0, 0.0], [2.0, 10.0, 0.0], [2.0, 0.0, 0.0]],
                vec![1.0; 3],
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            )
            .unwrap(),
        );
        let mut edge = Edge {
            start: curve.value(0.0),
            end: curve.value(0.1),
            curve,
            vertices: (0, 1),
            same_sense: true,
            samples: vec![],
        };
        let patch = Box2 {
            low: (0.0, 0.0),
            high: (1.0, 1.0),
        };
        assert_eq!(curve_error(&edge, &frame, patch), Some(0.0));
        edge.end = edge.curve.value(1.0);
        assert!(curve_error(&edge, &frame, patch).unwrap() > 0.9);
    }

    #[test]
    fn a_rotated_rational_circle_keeps_its_tangent_bound() {
        use crate::brep::RawCurve2d;
        let mut support = face(&[rectangle(-4.0, -4.0, 4.0, 4.0)]);
        let index = support.edges.len();
        let frame = Frame {
            origin: [-1.0, 0.5, 0.0],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        };
        support.edges.push(Edge {
            curve: Curve::Circle { frame, radius: 1.0 },
            start: [0.0, 0.5, 0.0],
            end: [0.0, 0.5, 0.0],
            vertices: (index, index),
            same_sense: true,
            samples: vec![],
        });
        support.faces[0].loops.push(Loop {
            edges: vec![(index, true)],
            vertex: None,
        });
        let (s, c) = 0.37_f64.sin_cos();
        let poles = [
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (-1.0, 1.0),
            (-1.0, 0.0),
            (-1.0, -1.0),
            (0.0, -1.0),
            (1.0, -1.0),
            (1.0, 0.0),
        ]
        .map(|(x, y)| {
            let world = [-1.0 + c * x - s * y, 0.5 + s * x + c * y, 0.0];
            let uv = support.faces[0].surface.parameters(world, None).unwrap();
            [uv.0, uv.1]
        })
        .to_vec();
        let pcurve = RawCurve2d::BSpline {
            degree: 2,
            points: poles,
            multiplicities: vec![3, 2, 2, 2, 3],
            knots: vec![0.0, 0.25, 0.5, 0.75, 1.0],
            weights: Some(
                (0..9)
                    .map(|i| {
                        if i % 2 == 0 {
                            1.0
                        } else {
                            std::f64::consts::FRAC_1_SQRT_2
                        }
                    })
                    .collect(),
            ),
        };
        let support = Part::new(support.faces.clone(), support.edges.clone(), vec![])
            .with_file_pcurves(
                std::collections::BTreeMap::from([(0, vec![(index, pcurve)])]),
                Default::default(),
                (1.0, 1.0),
                Some(1e-7),
            );
        let patch = face(&[rectangle(0.0, 0.0, 1.0, 1.0)]);
        assert_eq!(covered((&patch, 0), &[(&support, 0)]), Some(true));
        let overlapping = face(&[rectangle(-0.1, 0.0, 1.0, 1.0)]);
        assert_ne!(
            super::super::assess((&overlapping, 0), &[(&support, 0)]),
            super::super::Coverage::Covered
        );
    }
}
