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
    let (gap, error, tolerance) = finite_area(patch, supports)?;
    if gap + error <= tolerance {
        Some(true)
    } else if gap > tolerance + error {
        Some(false)
    } else {
        None
    }
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

fn finite_area(patch: FaceRef<'_>, supports: &[FaceRef<'_>]) -> Option<(f64, f64, f64)> {
    let frame = canonical_frame(patch)?;
    let target = finite_loops(patch, &frame, None)?;
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
    for &support in supports {
        if same_surface(
            &patch.0.faces[patch.1].surface,
            &support.0.faces[support.1].surface,
        ) {
            let boundary = finite_loops(support, &frame, Some(bounds))?;
            error += boundary.line_error;
            theirs.push(boundary.loops);
        }
    }
    Some((
        uncovered_area(&target.loops, &theirs, bounds)?,
        error,
        1e-9 * area,
    ))
}

fn finite_loops((part, face): FaceRef<'_>, frame: &Frame, clip: Option<Box2>) -> Option<Boundary> {
    let mut contours = Vec::new();
    let mut error = 0.0;
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
    for lp in &part.faces[face].loops {
        if lp.edges.is_empty() {
            return None;
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
            let alternatives = super::representations::alternatives((part, face), index)?;
            if !matches!(edge.curve, Curve::Line { .. }) {
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
            let Curve::Line { origin, .. } = edge.curve else {
                return None;
            };
            let mut size = origin.into_iter().map(f64::abs).fold(base_size, f64::max);
            for curve in &alternatives {
                let Curve::Line { origin, .. } = curve else {
                    return None;
                };
                size = origin.iter().copied().map(f64::abs).fold(size, f64::max);
            }
            let pad = 64.0 * f64::EPSILON * size;
            let world = [a, b].map(|p| edge.curve.value(edge.curve.parameter(p)));
            let feet = world.map(|p| project(frame, p));
            let mut radii = std::array::from_fn::<_, 2, _>(|i| {
                (vertices[i].0 - feet[i].0).hypot(vertices[i].1 - feet[i].1) + pad
            });
            let limit = part
                .file_uncertainty()
                .filter(|x| x.is_finite() && *x >= 0.0)
                .unwrap_or(0.0)
                + pad;
            if [a, b]
                .into_iter()
                .zip(world)
                .any(|(p, q)| geom::dist(p, q) > limit)
            {
                return None;
            }
            error += segment_error(feet, [0.0; 2], size, clip)?;
            for curve in &alternatives {
                let param_feet = [a, b].map(|p| {
                    let q = curve.value(curve.parameter(p));
                    (geom::dist(p, q) <= limit).then_some(project(frame, q))
                });
                let [Some(c), Some(d)] = param_feet else {
                    return None;
                };
                let alternate = [c, d];
                let deviations = std::array::from_fn(|i| {
                    radii[i] = radii[i].max(
                        (vertices[i].0 - alternate[i].0).hypot(vertices[i].1 - alternate[i].1)
                            + pad,
                    );
                    (feet[i].0 - alternate[i].0).hypot(feet[i].1 - alternate[i].1)
                });
                error += segment_error(feet, deviations, size, clip)?;
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
        Curve::Line { .. } => return Some(0.0),
        Curve::Nurbs(curve) => {
            if !curve.weights().iter().all(|w| w.is_finite() && *w > 0.0) {
                return None;
            }
            for &p in curve.control_points() {
                add(p)?;
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
        let support = crate::brep::Part::new(support.faces.clone(), support.edges.clone(), vec![])
            .with_file_pcurves(
                std::collections::BTreeMap::from([(
                    0,
                    vec![(
                        index,
                        RawCurve2d::Circle {
                            point: [0.5, 0.5],
                            reference: None,
                            radius: 0.25,
                        },
                    )],
                )]),
                Default::default(),
                (1.0, 1.0),
                Some(1e-7),
            );
        assert_eq!(covered((&patch, 0), &[(&support, 0)]), None);
        assert_ne!(
            super::super::assess((&patch, 0), &[(&support, 0)]),
            super::super::Coverage::Covered
        );
    }
}
