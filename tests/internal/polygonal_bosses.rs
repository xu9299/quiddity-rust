use super::common;
use crate::features::Context;
use crate::features::experimental_geometry::GeometryGraph;
use crate::features::polygonal_bosses::{PolygonalOptions, analyse};
use crate::read_step_file;

#[test]
fn face_inventory_defers_blend_tables_until_a_surface_query() {
    let part = read_step_file(&common::fixtures().join("golden_polygonal_boss.step")).unwrap();
    let graph = GeometryGraph::for_solid(&part, 0);
    let face = graph.faces()[0];
    assert!(graph.is_planar(face));
    assert!(!graph.neighbours(face).is_empty());
    assert!(!graph.has_blend_tables());
    assert!(graph.surface_fact(face).is_some());
    assert!(graph.has_blend_tables());
}

#[test]
fn coincident_bosses_keep_their_own_sides_and_free_end_cap_in_solid_order() {
    let path = common::fixtures().join("captured/662cf4f54c3b63bc.step.gz");
    let part = read_step_file(&path).unwrap();
    let result = analyse(&Context::new(&part), &PolygonalOptions::default());
    assert!(result.refusals.is_empty());
    let [first, second] = result.accepted.as_slice() else {
        panic!("expected one boss on each of the two coincident solids")
    };
    assert_eq!(first.record, second.record);
    assert_eq!(first.defining, [6, 7, 8, 9, 10, 11]);
    assert_eq!(first.context, [12]);
    assert_eq!(second.defining, [19, 20, 21, 22, 23, 24]);
    assert_eq!(second.context, [25]);
}
