//! Public corpus faces exercising STEP boundary selection at the geometry API.
mod common;

use haecceity::{
    geom, read_step_file,
    step::{IDENTITY, read_step_file_placed},
};

#[test]
fn conical_rims_use_the_geometry_of_the_selected_pcurves() {
    let Some(dir) = common::corpus_dir() else {
        assert!(std::env::var_os("QUIDDITY_CORPUS_REQUIRED").is_none());
        return;
    };
    let part = read_step_file(&dir.join("cadgenbench_inputs/cgb243.step.gz")).unwrap();
    // A nearly level generator coordinate on a cone, a planar rational
    // quadratic, and a nearly level generator coordinate on a cylinder.
    for (entity, expected) in [
        (5444, [3.5, 4.]),
        (56072, [32., 32.5]),
        (57220, [13.5, 14.]),
    ] {
        let face = (0..part.faces.len())
            .find(|&f| part.face_source(f).is_some_and(|s| s.entity == entity))
            .unwrap();
        let mut radii: Vec<_> = part
            .face_edges(face)
            .into_iter()
            .filter_map(|e| match part.boundary_edge_support(e)? {
                geom::Curve::Circle { radius, .. } => Some(*radius),
                _ => None,
            })
            .collect();
        radii.sort_by(f64::total_cmp);
        assert_eq!(radii.len(), 2, "face #{entity}");
        assert!(
            radii
                .iter()
                .zip(expected)
                .all(|(r, w)| (r - w).abs() < geom::COORD_FLOOR),
            "face #{entity}: {radii:?}"
        );
    }
}

#[test]
fn circular_feature_support_keeps_the_selected_spline_and_its_interval() {
    let captured = common::fixtures().join("captured");
    for (file, entities) in [
        ("0d84b2b05e79d24b.step.gz", [265, 328, 465, 518]),
        ("fa0600fbcc4940d4.step.gz", [398, 452, 632, 685]),
    ] {
        let part = read_step_file(&captured.join(file)).unwrap();
        for entity in entities {
            let edge = (0..part.edges.len())
                .find(|&e| part.edge_source(e).is_some_and(|s| s.entity == entity))
                .unwrap();
            let Some(geom::Curve::Nurbs(spline)) = part.boundary_edge_curve(edge) else {
                panic!("{file} #{entity}: expected the selected rational spline");
            };
            let Some(geom::Curve::Circle { frame, radius }) = part.boundary_edge_support(edge)
            else {
                panic!("{file} #{entity}: missing circular rim");
            };
            assert!((radius - 5.).abs() < 1e-10);
            let (a, b) = part.boundary_edge_interval(edge).unwrap();
            assert!(geom::dist(spline.value(a), part.edges[edge].start) < 1e-8);
            assert!(geom::dist(spline.value(b), part.edges[edge].end) < 1e-8);
            let q = frame.to_local(spline.value(0.5 * (a + b)));
            assert!((q[0].hypot(q[1]) - radius).abs() < 1e-10 && q[2].abs() < 1e-10);
        }
    }
}

#[test]
fn sweeping_a_selected_quarter_arc_keeps_its_finite_side_and_empty_cut() {
    use haecceity::{
        classify::Classifier,
        rays::RayCaster,
        sweep::extrude_face,
        volume::{Probe, common_volume},
    };
    let file = common::fixtures().join("captured/85e5cdbaee3f051b.step.gz");
    let part = read_step_file(&file).unwrap();
    assert_eq!(part.face_source(6).unwrap().entity, 473);
    let sweep = extrude_face(&part, 6, [0.; 3], [-25., 0., 0.]).unwrap();
    let rays = RayCaster::for_solid(&sweep, 0);
    let origin = [-19.914436919, 11.004339173, -11.];
    let hits = rays.trimmed_hits(origin, [0., 0., 1.], 24.).unwrap();
    let dy: f64 = origin[1] - 15.;
    let expected = 21. - (16. - dy * dy).sqrt();
    assert_eq!(hits.len(), 2);
    assert!((hits[0].t - expected).abs() < 1e-8 && (hits[1].t - 21.).abs() < 1e-10);
    let (volume, _) = sweep.solid_mass(0).unwrap();
    assert!((volume - 100. * std::f64::consts::PI).abs() < 1e-8);
    assert_eq!(
        common_volume(&Classifier::for_solid(&part, 0), &Probe::Solid(rays)),
        Some(0.)
    );
}

#[test]
fn a_closed_spline_hole_sweeps_to_a_complete_band_at_its_knot_seam() {
    use haecceity::{
        classify::Classifier,
        geom::Surface,
        rays::RayCaster,
        sweep::extrude_face,
        volume::{Probe, common_volume},
    };
    let Some(dir) = common::corpus_dir() else {
        assert!(std::env::var_os("QUIDDITY_CORPUS_REQUIRED").is_none());
        return;
    };
    let part = read_step_file(&dir.join("mfcadpp/10155.step")).unwrap();
    assert_eq!(part.face_source(24).unwrap().entity, 2387);
    let depth = 7.98651;
    let inset = 1e-6 * (1. + depth);
    let height = depth - 2. * inset;
    let probe = extrude_face(&part, 24, [0., 0., -inset], [0., 0., -height]).unwrap();
    let side = probe
        .faces
        .iter()
        .position(|f| matches!(f.surface, Surface::Freeform { .. }) && f.loops.len() == 2)
        .unwrap();
    // Each rim ends at its starting STEP vertex after one full turn.
    // The side area is perimeter times height, and volume is cap area times height.
    let rim = part
        .face_edges(24)
        .into_iter()
        .find(|&e| part.edges[e].is_closed())
        .unwrap();
    let geom::Curve::Circle { radius, .. } = part.edges[rim].curve else {
        panic!("the file supplies an independent 3D circle for the hole");
    };
    let expected_side = std::f64::consts::TAU * radius * height;
    assert!((probe.face_mass(side).unwrap()[0] - expected_side).abs() < 1e-8);
    let expected_volume = part.face_mass(24).unwrap()[0] * height;
    assert!((probe.solid_mass(0).unwrap().0 - expected_volume).abs() < 1e-8);
    assert_eq!(
        common_volume(
            &Classifier::for_solid(&part, 0),
            &Probe::Solid(RayCaster::for_solid(&probe, 0))
        ),
        Some(0.)
    );
}

#[test]
fn trim_moments_survive_placement() {
    let Some(dir) = common::corpus_dir() else {
        assert!(std::env::var_os("QUIDDITY_CORPUS_REQUIRED").is_none());
        return;
    };
    let mut moved = IDENTITY;
    for (row, offset) in moved.iter_mut().zip([123.456, -78.9, 41.3]) {
        row[3] = offset;
    }
    for (name, faces) in [
        ("cadgenbench/threaded_connector_109.step", &[17, 139][..]),
        ("cadgenbench_inputs/cgb242.step.gz", &[298, 304][..]),
    ] {
        let native = read_step_file(&dir.join(name)).unwrap();
        let placed = read_step_file_placed(&dir.join(name), &moved).unwrap();
        for &f in faces {
            let a = native
                .face_moments(f)
                .unwrap_or_else(|| panic!("{name} face {f}, native"));
            let b = placed
                .face_moments(f)
                .unwrap_or_else(|| panic!("{name} face {f}, translated"));
            assert!(
                (a.area - b.area).abs() < 1e-9 * a.area + 1e-9,
                "{name} face {f}"
            );
            assert!(
                geom::dist(geom::add(a.centroid, moved.map(|r| r[3])), b.centroid) < 1e-6,
                "{name} face {f}"
            );
        }
    }
}

#[test]
fn an_oblique_cylinder_cut_has_the_same_domain_on_both_sides_of_its_seam() {
    let Some(dir) = common::corpus_dir() else {
        assert!(std::env::var_os("QUIDDITY_CORPUS_REQUIRED").is_none());
        return;
    };
    // Two planes bound each cylindrical wall. At u = 0 their intersections
    // place the probe below outside the oblique end, and v = -5 inside.
    for (file, entity, below) in [("10096", 1989, -17.), ("10138", 2636, -20.33)] {
        let part = read_step_file(&dir.join(format!("mfcadpp/{file}.step"))).unwrap();
        let f = (0..part.faces.len())
            .find(|&f| part.face_source(f).unwrap().entity == entity)
            .unwrap();
        let domain = part.domain(f).unwrap();
        for u in [-1e-10, 0., 1e-10, std::f64::consts::TAU] {
            assert!(!domain.contains(u, below), "{file} #{entity} {u}");
            assert!(domain.contains(u, -5.), "{file} #{entity} {u}");
        }
        if file == "10096" {
            assert!(domain.contains(std::f64::consts::PI, below));
        }
    }
}

#[test]
fn a_shared_vertex_closes_the_declared_mesh_loop() {
    let Some(dir) = common::corpus_dir() else {
        assert!(std::env::var_os("QUIDDITY_CORPUS_REQUIRED").is_none());
        return;
    };
    let mut moved = IDENTITY;
    for (row, offset) in moved.iter_mut().zip([123.456, -78.9, 41.3]) {
        row[3] = offset;
    }
    let file = dir.join("cadgenbench/flanged_spool_132.step");
    for m in [IDENTITY, moved] {
        let part = read_step_file_placed(&file, &m).unwrap();
        let f = 2;
        assert_eq!(part.face_source(f).unwrap().entity, 182);
        let anchor = part.face_anchor(f).unwrap();
        assert!(part.domain(f).unwrap().contains(anchor.uv.0, anchor.uv.1));
        let mesh = part
            .triangulate(f, 1e-3 * part.bounds().diagonal(), 0.3)
            .unwrap();
        assert!(!mesh.triangles.is_empty());
    }
}

#[test]
fn the_winding_ray_through_a_closed_pcurve_end_reads_the_annular_groove() {
    let Some(dir) = common::corpus_dir() else {
        assert!(std::env::var_os("QUIDDITY_CORPUS_REQUIRED").is_none());
        return;
    };
    let part = read_step_file(&dir.join("mfcadpp_holdout/419.step")).unwrap();
    let floor = 28;
    assert_eq!(part.face_source(floor).unwrap().entity, 2005);
    let domain = part.domain(floor).unwrap();
    // The two concentric STEP cylinders bound radii 1.350424756018 and
    // 2.071143801424 on this plane. Straddle the stored pcurve's endpoint
    // and test material, the central hole and exterior on both sides.
    for u in [-1e-8_f64, 0., 1e-13, 3e-13, 1e-8] {
        for v in [-2.08_f64, -2.07, -1.36, -1.34, 0., 1.34, 1.36, 2.07, 2.08] {
            let r = u.hypot(v);
            assert_eq!(
                domain.contains(u, v),
                (1.350424756018..2.071143801424).contains(&r),
                "{u} {v}"
            );
        }
    }
    for wall in [26, 30] {
        assert_eq!(part.arc(floor, wall), Some(haecceity::brep::Arc::Concave));
        assert_eq!(part.arc(wall, floor), Some(haecceity::brep::Arc::Concave));
    }
}

#[test]
fn a_closed_master_repeats_its_knots_on_the_measured_turn() {
    let file =
        common::fixtures().join("captured/effective_surfaces/parts/201ebceb4462a612.step.gz");
    let part = read_step_file(&file).unwrap();
    let e = (0..part.edges.len())
        .find(|&e| part.edge_source(e).unwrap().entity == 374)
        .unwrap();
    let Some(geom::Curve::Nurbs(curve)) = part.boundary_edge_curve(e) else {
        panic!("the master on the plane is a rational quadratic")
    };
    let (a, b) = part.boundary_edge_interval(e).unwrap();
    assert!((a - curve.domain().1).abs() < 1e-10 && b > curve.domain().1);
    // The three 120-degree rational arcs describe the circle of radius 8.
    let length = part.declared_edge_length(e).unwrap();
    assert!(
        (length - 16. * std::f64::consts::PI).abs() < 1e-9,
        "{length}"
    );
    for fraction in [1. / 6., 0.25, 0.5, 0.875] {
        let (p, _) = part.declared_edge_at_length(e, fraction).unwrap();
        let chord = geom::dist(part.edges[e].start, p);
        assert!((chord - 16. * (fraction * std::f64::consts::PI).sin()).abs() < 1e-9);
    }
}
