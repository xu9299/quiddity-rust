//! Public corpus faces exercising STEP boundary selection at the geometry API.
mod common;

use haecceity::{
    geom, read_step_file,
    step::{IDENTITY, read_step_file_placed},
};

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
