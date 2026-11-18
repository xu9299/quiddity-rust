//! A bounded straight planar query needs no sampling chord margin. Literal finite
//! polygons provide the integration reference; complete edge homotopies enclose
//! the differences from their analytic and available file pcurve representations.
//! The relative area budget is unchanged. An indecisive or ineligible query keeps
//! the existing sampled path, including every curved and cylindrical query.

use super::interval::Interval as I;
use super::{FaceRef, same_surface, straight_intervals as straight};
use crate::brep::{Pcurve, RawCurve2d};
use crate::geom::{self, Curve, Frame, Surface, V3};

type Point = (f64, f64);
type Loops = Vec<Vec<Point>>;
type Exact = [[I; 2]; 2];
type Bounds = (Point, Point);

struct Boundary {
    loops: Loops,
    error: f64,
}

pub(super) fn covered(patch: FaceRef<'_>, supports: &[FaceRef<'_>]) -> Option<bool> {
    if !matches!(patch.0.faces[patch.1].surface, Surface::Plane { .. }) {
        return None;
    }
    let relevant: Vec<_> = supports
        .iter()
        .copied()
        .filter(|s| same_surface(&patch.0.faces[patch.1].surface, &s.0.faces[s.1].surface))
        .collect();
    // Bound the new proof's work before constructing any pcurve image. These
    // limits select a proof path; they do not change a geometric tolerance.
    let mut count = 0usize;
    for (part, face) in std::iter::once(patch).chain(relevant.iter().copied()) {
        for &(index, _) in part.faces[face].loops.iter().flat_map(|lp| &lp.edges) {
            count += 1;
            if count > 512 || !eligible((part, face), index) {
                return None;
            }
        }
    }
    let frame = canonical_frame(patch)?;
    let target = boundary(patch, &frame, None)?;
    let low = target
        .loops
        .iter()
        .flatten()
        .fold((f64::INFINITY, f64::INFINITY), |a, p| {
            (a.0.min(p.0), a.1.min(p.1))
        });
    let high = target
        .loops
        .iter()
        .flatten()
        .fold((f64::NEG_INFINITY, f64::NEG_INFINITY), |a, p| {
            (a.0.max(p.0), a.1.max(p.1))
        });
    let bounds = (low, high);
    let mut theirs = Vec::new();
    let mut error = I::point(target.error);
    for support in relevant {
        // A failed representation cannot be dropped to prove a rejection.
        let other = boundary(support, &frame, Some(bounds))?;
        theirs.push(other.loops);
        error = error + I::point(other.error);
    }
    let area = super::planar_intervals::area(&target.loops, &[], bounds)?;
    let gap = super::planar_intervals::area(&target.loops, &theirs, bounds)?;
    let lower = (area - I::point(target.error)).lo;
    let upper = (area + I::point(target.error)).hi;
    if lower > 0.0 && (gap + error).hi <= (I::point(1e-9) * I::point(lower)).lo {
        Some(true)
    } else if (gap - error).lo > (I::point(1e-9) * I::point(upper)).hi {
        Some(false)
    } else {
        None
    }
}

/// Check representation kinds before paying for any projected interval or
/// homotopy. Most curved queries can therefore retain their old path immediately.
fn eligible((part, face): FaceRef<'_>, index: usize) -> bool {
    let line = match &part.edges[index].curve {
        Curve::Line { .. } => true,
        Curve::Nurbs(c) => {
            c.degree() == 1
                && c.control_points().len() == 2
                && c.weights().iter().all(|w| w.is_finite() && *w > 0.0)
        }
        _ => false,
    };
    if !line {
        return false;
    }
    let straight = |raw: &RawCurve2d, surface: &Surface| {
        if !matches!(surface, Surface::Plane { .. }) {
            return false;
        }
        match raw {
            RawCurve2d::Line { .. } => true,
            RawCurve2d::BSpline {
                degree,
                points,
                multiplicities,
                knots,
                weights,
            } => {
                *degree == 1
                    && points.len() == 2
                    && multiplicities == &[2, 2]
                    && knots.len() == 2
                    && knots.iter().all(|x| x.is_finite())
                    && knots[0] < knots[1]
                    && weights
                        .as_ref()
                        .is_none_or(|w| w.len() == 2 && w.iter().all(|x| x.is_finite() && *x > 0.0))
            }
        }
    };
    if part
        .file_pcurves(face)
        .iter()
        .any(|(e, pc)| *e == index && !straight(pc, &part.faces[face].surface))
    {
        return false;
    }
    if part.file_pcurves(face).is_empty()
        && part.faces[face]
            .pcurves
            .iter()
            .any(|(e, pc)| *e == index && !matches!(pc, Pcurve::Line { .. }))
    {
        return false;
    }
    for &adjacent in &part.edge_faces()[index] {
        if adjacent != face
            && part
                .primary_pcurves(adjacent)
                .any(|(e, pc)| *e == index && !straight(pc, &part.faces[adjacent].surface))
        {
            return false;
        }
    }
    true
}

fn canonical_frame((part, face): FaceRef<'_>) -> Option<Frame> {
    let Surface::Plane { frame } = part.faces[face].surface else {
        return None;
    };
    let mut z = geom::unit(frame.z)?;
    if z[geom::dominant_axis(z)] < 0.0 {
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

fn alternatives(
    (part, face): FaceRef<'_>,
    index: usize,
    vertices: [V3; 2],
) -> Option<Vec<straight::Ends>> {
    let edge = &part.edges[index];
    let mut images = vec![straight::curve_ends(&edge.curve, vertices)?];
    let units = part.file_units().unwrap_or((1.0, 1.0));
    if !units.0.is_finite() || units.0 <= 0.0 {
        return None;
    }
    let mut add = |raw: &RawCurve2d, surface: &Surface, units| {
        images.push(straight::pcurve_ends(raw, surface, units, vertices)?);
        Some(())
    };
    for (e, raw) in part.file_pcurves(face) {
        if *e == index {
            add(raw, &part.faces[face].surface, units)?;
        }
    }
    if part.file_pcurves(face).is_empty() {
        for (e, pc) in &part.faces[face].pcurves {
            if *e != index {
                continue;
            }
            let Pcurve::Line { point, dir } = pc else {
                return None;
            };
            add(
                &RawCurve2d::Line {
                    point: [point.0, point.1],
                    direction: [dir.0, dir.1],
                    magnitude: 1.0,
                },
                &part.faces[face].surface,
                (1.0, 1.0),
            )?;
        }
    }
    for &adjacent in &part.edge_faces()[index] {
        if adjacent == face {
            continue;
        }
        for (e, raw) in part.primary_pcurves(adjacent) {
            if *e == index {
                add(raw, &part.faces[adjacent].surface, units)?;
            }
        }
    }
    Some(images)
}

fn boundary((part, face): FaceRef<'_>, frame: &Frame, clip: Option<Bounds>) -> Option<Boundary> {
    let mut loops = Vec::new();
    let mut error = I::point(0.0);
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
            let vertices = if forward {
                [edge.start, edge.end]
            } else {
                [edge.end, edge.start]
            };
            if last.is_some_and(|v| v != vertices[0]) {
                return None;
            }
            first.get_or_insert(vertices[0]);
            last = Some(vertices[1]);
            let images = alternatives((part, face), index, vertices)?;
            let exact_vertices = vertices.map(straight::point);
            let exact = exact_vertices.map(|p| straight::project(frame, p));
            let rounded = vertices.map(|p| {
                let d = geom::sub(p, frame.origin);
                (geom::dot(d, frame.x), geom::dot(d, frame.y))
            });
            let base = rounded.map(|p| [I::point(p.0), I::point(p.1)]);
            error = error + I::point(swept_area(base, exact, clip)?);
            for ends in images {
                if let Some(limit) = part
                    .file_uncertainty()
                    .filter(|v| v.is_finite() && *v >= 0.0)
                {
                    for (p, q) in exact_vertices.into_iter().zip(ends) {
                        if straight::norm_upper(std::array::from_fn(|k| p[k] - q[k])) > limit {
                            return None;
                        }
                    }
                }
                let mapped = ends.map(|p| straight::project(frame, p));
                error = error + I::point(swept_area(exact, mapped, clip)?);
            }
            points.push(rounded[0]);
        }
        if first != last {
            return None;
        }
        points.push(*points.first()?);
        loops.push(points);
    }
    (!loops.is_empty() && error.hi.is_finite()).then_some(Boundary {
        loops,
        error: error.hi,
    })
}

/// The straight homotopy's affine Jacobian is bounded by the mean absolute
/// determinant at its four corners. Its complete image also lies in the corner
/// box. Only support error is clipped to the nominal target box: target error
/// already pays for every point outside that box.
fn swept_area(a: Exact, b: Exact, clip: Option<Bounds>) -> Option<f64> {
    let delta = |p: [I; 2], q: [I; 2]| [p[0] - q[0], p[1] - q[1]];
    let directions = [delta(a[1], a[0]), delta(b[1], b[0])];
    let deviations = [delta(b[0], a[0]), delta(b[1], a[1])];
    let mut area = I::point(0.0);
    for d in directions {
        for v in deviations {
            area = area + I::point((d[0] * v[1] - d[1] * v[0]).abs_upper()) * I::point(0.25);
        }
    }
    let mut bound = area.hi;
    if let Some((lo, hi)) = clip {
        let mut low = [f64::INFINITY; 2];
        let mut high = [f64::NEG_INFINITY; 2];
        for p in a.into_iter().chain(b) {
            for k in 0..2 {
                low[k] = low[k].min(p[k].lo);
                high[k] = high[k].max(p[k].hi);
            }
        }
        let width = (I::point(high[0].min(hi.0)) - I::point(low[0].max(lo.0)))
            .hi
            .max(0.0);
        let height = (I::point(high[1].min(hi.1)) - I::point(low[1].max(lo.1)))
            .hi
            .max(0.0);
        bound = bound.min((I::point(width) * I::point(height)).hi);
    }
    (bound.is_finite() && bound >= 0.0).then_some(bound)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::{Edge, Face, Loop, Part};

    fn polygon(contours: &[Vec<Point>], angle: f64, offset: f64, start: usize) -> Part {
        let (s, c) = angle.sin_cos();
        let placed = |p: Point| [c * p.0 - s * p.1 + offset, s * p.0 + c * p.1 - offset, 0.0];
        let mut edges = Vec::new();
        let mut loops = Vec::new();
        for contour in contours {
            let mut points = contour.clone();
            points.rotate_left(start % contour.len());
            let base = edges.len();
            for i in 0..points.len() {
                let (a, b) = (placed(points[i]), placed(points[(i + 1) % points.len()]));
                edges.push(Edge {
                    curve: Curve::Line {
                        origin: a,
                        dir: geom::sub(b, a),
                    },
                    start: a,
                    end: b,
                    vertices: (base + i, base + (i + 1) % points.len()),
                    same_sense: true,
                    samples: vec![a, b],
                });
            }
            loops.push(Loop {
                edges: (base..edges.len()).map(|i| (i, true)).collect(),
                vertex: None,
            });
        }
        let first = &edges[0];
        let x = geom::unit(geom::sub(first.end, first.start)).unwrap();
        Part::new(
            vec![Face {
                surface: Surface::Plane {
                    frame: Frame {
                        origin: first.start,
                        x,
                        y: [-x[1], x[0], 0.0],
                        z: [0.0, 0.0, 1.0],
                    },
                },
                reversed: false,
                loops,
                solid: None,
                pcurves: vec![],
            }],
            edges,
            vec![],
        )
    }

    fn rectangle(a: f64, b: f64, c: f64, d: f64) -> Vec<Point> {
        vec![(a, b), (c, b), (c, d), (a, d)]
    }

    #[test]
    fn a_real_strip_gap_cannot_be_hidden_by_restarting_the_patch_loop() {
        for angle in [0.0, 0.37, 1.2] {
            for offset in [0.0, 123.4567] {
                let support = polygon(&[rectangle(0.0, 0.0, 1.0, 0.9999)], angle, offset, 0);
                for start in 0..4 {
                    let target = polygon(&[rectangle(0.0, 0.0, 1.0, 1.0)], angle, offset, start);
                    assert_eq!(covered((&target, 0), &[(&support, 0)]), Some(false));
                    assert!(!super::super::covered((&target, 0), &[(&support, 0)]));
                }
            }
        }
    }

    #[test]
    fn complete_support_holes_and_disjoint_fillers_keep_their_exact_areas() {
        let outer = rectangle(0.0, 0.0, 1.0, 1.0);
        let inner = rectangle(0.25, 0.25, 0.75, 0.75);
        for angle in [0.0, 0.37, 1.2] {
            let target = polygon(std::slice::from_ref(&outer), angle, 0.0, 0);
            let ring = polygon(&[outer.clone(), inner.clone()], angle, 0.0, 0);
            let fill = polygon(std::slice::from_ref(&inner), angle, 0.0, 0);
            assert_eq!(covered((&target, 0), &[(&target, 0)]), Some(true));
            assert_eq!(covered((&target, 0), &[(&ring, 0)]), Some(false));
            assert_eq!(covered((&target, 0), &[(&ring, 0), (&fill, 0)]), Some(true));
        }
    }

    #[test]
    fn crossing_edges_are_integrated_between_their_actual_intersections() {
        let square = polygon(&[rectangle(0.0, 0.0, 1.0, 1.0)], 0.37, 0.0, 0);
        let triangle = polygon(&[vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]], 0.37, 0.0, 0);
        assert_eq!(covered((&square, 0), &[(&triangle, 0)]), Some(false));
        let other = polygon(&[vec![(1.0, 1.0), (0.0, 1.0), (1.0, 0.0)]], 0.37, 0.0, 0);
        assert_eq!(
            covered((&square, 0), &[(&triangle, 0), (&other, 0)]),
            Some(true)
        );
    }

    #[test]
    fn a_pcurve_with_a_shorter_finite_domain_cannot_supply_missing_material() {
        let target = polygon(&[rectangle(0.0, 0.0, 1.0, 1.0)], 0.0, 0.0, 0);
        let support = polygon(&[rectangle(0.0, 0.0, 1.0, 1.0)], 0.0, 0.0, 0).with_file_pcurves(
            std::collections::BTreeMap::from([(
                0,
                vec![(
                    0,
                    RawCurve2d::BSpline {
                        degree: 1,
                        points: vec![[0.0, 0.0], [0.5, 0.0]],
                        multiplicities: vec![2, 2],
                        knots: vec![0.0, 1.0],
                        weights: None,
                    },
                )],
            )]),
            Default::default(),
            (1.0, 1.0),
            Some(1e-7),
        );
        assert_eq!(covered((&target, 0), &[(&support, 0)]), None);
    }

    #[test]
    fn a_mismatched_file_line_is_not_an_exact_topological_chord() {
        let target = polygon(&[rectangle(0.0, 0.0, 1.0, 1.0)], 0.0, 0.0, 0);
        let mut support = polygon(&[rectangle(0.0, 0.0, 1.0, 1.0)], 0.0, 0.0, 0);
        if let Curve::Line { origin, .. } = &mut support.edges[2].curve {
            origin[1] -= 0.001;
        }
        assert_ne!(covered((&target, 0), &[(&support, 0)]), Some(true));
    }

    #[test]
    fn curved_edges_keep_the_existing_query_path() {
        let target = polygon(&[rectangle(0.0, 0.0, 1.0, 1.0)], 0.0, 0.0, 0);
        let mut support = polygon(&[rectangle(0.0, 0.0, 1.0, 1.0)], 0.0, 0.0, 0);
        support.edges[0].curve = Curve::Circle {
            frame: Frame {
                origin: [0.5, 0.0, 0.0],
                x: [1.0, 0.0, 0.0],
                y: [0.0, 1.0, 0.0],
                z: [0.0, 0.0, 1.0],
            },
            radius: 0.5,
        };
        assert_eq!(covered((&target, 0), &[(&support, 0)]), None);
    }
}
