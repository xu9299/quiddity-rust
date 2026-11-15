//! Independent geometric expectations for cached edge bounds, including periodic seams.
use super::*;
use crate::geom::{Curve, Frame};
use std::sync::OnceLock;

fn edge(curve: Curve, start: V3, end: V3, same_sense: bool, closed: bool) -> Edge {
    Edge {
        curve,
        start,
        end,
        same_sense,
        vertices: (0, usize::from(!closed)),
        samples: vec![start, end],
    }
}

fn check(edge: &Edge, cache: &OnceLock<EdgeInterval>, q: V3, want: V3, between: bool) {
    let (point, actual_between) = nearest_on_edge(edge, q, cache);
    assert!(geom::dist(point, want) < 1e-10, "{point:?} != {want:?}");
    assert_eq!(actual_between, between);
}

#[test]
fn trimmed_line_keeps_interior_and_endpoint_decisions() {
    for reverse in [false, true] {
        let (a, b) = if reverse { (4.0, 2.0) } else { (2.0, 4.0) };
        let edge = edge(
            Curve::Line {
                origin: [0.0; 3],
                dir: [1.0, 0.0, 0.0],
            },
            [a, 0.0, 0.0],
            [b, 0.0, 0.0],
            !reverse,
            false,
        );
        let cache = OnceLock::new();
        for _ in 0..2 {
            check(&edge, &cache, [3.0, 5.0, 0.0], [3.0, 0.0, 0.0], true);
            check(&edge, &cache, [0.0, 1.0, 0.0], [2.0, 0.0, 0.0], false);
            check(&edge, &cache, [6.0, 1.0, 0.0], [4.0, 0.0, 0.0], false);
        }
    }
}

#[test]
fn periodic_arc_crosses_seam_and_closed_edge_skips_bounds() {
    let curve = Curve::Circle {
        frame: Frame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        },
        radius: 1.0,
    };
    for reverse in [false, true] {
        let (start, end) = if reverse {
            ([0.0, 1.0, 0.0], [0.0, -1.0, 0.0])
        } else {
            ([0.0, -1.0, 0.0], [0.0, 1.0, 0.0])
        };
        let arc = edge(curve.clone(), start, end, !reverse, false);
        let cache = OnceLock::new();
        for _ in 0..2 {
            check(&arc, &cache, [3.0, 0.0, 0.0], [1.0, 0.0, 0.0], true);
            check(&arc, &cache, [-3.0, 0.2, 0.0], [0.0, 1.0, 0.0], false);
        }
    }
    let closed = edge(curve, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], true, true);
    let cache = OnceLock::new();
    check(&closed, &cache, [-2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], true);
    assert!(cache.get().is_none());
}

#[test]
fn trimmed_quadratic_spline_retains_original_curve_geometry() {
    // Quadratic Bezier with x=2t and y=2t(1-t), trimmed to t in [1/4, 3/4].
    // The point above the symmetry axis is nearest the apex; remote side points
    // have an unconstrained projection outside the trim and select its endpoints.
    let curve = Curve::Nurbs(
        crate::nurbs::NurbsCurve::new(
            2,
            vec![[0.0; 3], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0; 3],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        )
        .unwrap(),
    );
    let edge = edge(curve, [0.5, 0.375, 0.0], [1.5, 0.375, 0.0], true, false);
    let cache = OnceLock::new();
    for _ in 0..2 {
        check(&edge, &cache, [1.0, 2.0, 0.0], [1.0, 0.5, 0.0], true);
        check(&edge, &cache, [-2.0, 0.0, 0.0], edge.start, false);
        check(&edge, &cache, [4.0, 0.0, 0.0], edge.end, false);
    }
}

#[test]
fn populated_intervals_follow_part_transform_and_repair_lifecycle() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/two_solids.step");
    let loaded = crate::read_step_file(&path).unwrap();
    let mut part = Part::new(
        loaded.faces.clone(),
        loaded.edges.clone(),
        loaded.solids.clone(),
    );
    let index = part
        .edges
        .iter()
        .position(|e| matches!(e.curve, Curve::Line { .. }) && !e.is_closed())
        .unwrap();
    let rays = RayCaster::for_solid(&part, 0);
    nearest_on_edge(
        &part.edges[index],
        part.edges[index].start,
        &rays.intervals[index],
    );
    assert!(rays.intervals[index].get().is_some());
    let moved = part
        .transformed(&[
            [1.0, 0.0, 0.0, 10.0],
            [0.0, 1.0, 0.0, -7.0],
            [0.0, 0.0, 1.0, 4.0],
        ])
        .unwrap();
    let moved_rays = RayCaster::for_solid(&moved, 0);
    assert!(moved_rays.intervals[index].get().is_none());
    check(
        &moved.edges[index],
        &moved_rays.intervals[index],
        moved.edges[index].start,
        moved.edges[index].start,
        true,
    );
    drop(rays);
    part.record_inverted_shell_repair(0, true);
    let repaired = RayCaster::for_solid(&part, 0);
    assert!(repaired.intervals[index].get().is_none());
    check(
        &part.edges[index],
        &repaired.intervals[index],
        part.edges[index].start,
        part.edges[index].start,
        true,
    );
}

#[test]
#[ignore = "manual fixed-bound inversion timing; no wall-clock assertion"]
fn compare_repeated_endpoint_inversion_with_cached_bounds() {
    let curve = Curve::Nurbs(
        crate::nurbs::NurbsCurve::new(
            2,
            vec![[0.0; 3], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0; 3],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        )
        .unwrap(),
    );
    let edge = edge(curve, [0.5, 0.375, 0.0], [1.5, 0.375, 0.0], true, false);
    let queries: Vec<V3> = (0..20_000)
        .map(|i| [-1.0 + 4.0 * i as f64 / 19_999.0, 1.5, 0.0])
        .collect();
    // The previous operation sequence, including each query's two endpoint inversions.
    let original = |q: V3| {
        let t = edge.curve.parameter(q);
        let (a, b) = edge_interval(&edge.curve, edge.start, edge.end, edge.same_sense, false);
        let (lo, hi) = (a.min(b), a.max(b));
        let t = edge
            .curve
            .period()
            .map_or(t, |period| lo + (t - lo).rem_euclid(period));
        if (lo..=hi).contains(&t) {
            (edge.curve.value(t), true)
        } else if geom::dist(edge.start, q) <= geom::dist(edge.end, q) {
            (edge.start, false)
        } else {
            (edge.end, false)
        }
    };
    let bits = |(point, interior): (V3, bool)| (point.map(f64::to_bits), interior);
    let started = std::time::Instant::now();
    let expected: Vec<_> = queries.iter().map(|&q| bits(original(q))).collect();
    let original_seconds = started.elapsed().as_secs_f64();
    let cache = OnceLock::new();
    let started = std::time::Instant::now();
    let actual: Vec<_> = queries
        .iter()
        .map(|&q| bits(nearest_on_edge(&edge, q, &cache)))
        .collect();
    let cached_seconds = started.elapsed().as_secs_f64();
    assert_eq!(actual, expected);
    eprintln!(
        "Synthetic quadratic spline: {} queries; original={original_seconds:.6}s; cached={cached_seconds:.6}s",
        queries.len()
    );
}
