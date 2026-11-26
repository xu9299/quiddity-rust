//! Reproduce the narrow-face report against unchanged, frozen library builds.

use std::hash::{Hash, Hasher};

const FACES: &[usize] = &[3, 16, 55, 112, 136, 138, 141, 142, 148];

fn signature(mesh: &haecceity::mesh::FaceMesh) -> u64 {
    let mut hash = std::collections::hash_map::DefaultHasher::new();
    for p in &mesh.points {
        p.map(f64::to_bits).hash(&mut hash);
    }
    for &(u, v) in &mesh.uv {
        (u.to_bits(), v.to_bits()).hash(&mut hash);
    }
    mesh.triangles.hash(&mut hash);
    mesh.folded_dropped.0.hash(&mut hash);
    mesh.folded_dropped.1.to_bits().hash(&mut hash);
    hash.finish()
}

fn run(path: &std::path::Path) {
    let part = haecceity::read_step_file(path).unwrap();
    let diagonal = part.bounds().diagonal();
    assert!(diagonal > 0.0 && diagonal.is_finite());
    let deflection = 0.001 * diagonal;
    let meshes = part.mesh(deflection, 0.3);
    println!(
        "{{\"kind\":\"whole_part\",\"faces\":{},\"deflection\":{:.17e},\"angular\":0.3,\"triangles\":{},\"refusals\":{}}}",
        part.faces.len(),
        deflection,
        meshes
            .iter()
            .filter_map(|m| m.as_ref().ok())
            .map(|m| m.triangles.len())
            .sum::<usize>(),
        meshes.iter().filter(|m| m.is_err()).count()
    );
    for &face in FACES {
        let area = part.face_moments(face).unwrap().area;
        assert!(area > 0.0 && area.is_finite());
        let entity = part.face_source(face).unwrap().entity;
        let mesh = meshes[face].as_ref().unwrap();
        println!(
            "{{\"kind\":\"face\",\"face\":{},\"entity\":{},\"reference_area\":{:.17e},\"mesh_area\":{:.17e},\"area_error_percent\":{:.17e},\"triangles\":{},\"folded_dropped\":{},\"folded_area\":{:.17e},\"converged\":{},\"fingerprint\":\"{:016X}\"}}",
            face,
            entity,
            area,
            mesh.area(),
            100.0 * (mesh.area() / area - 1.0),
            mesh.triangles.len(),
            mesh.folded_dropped.0,
            mesh.folded_dropped.1,
            mesh.converged,
            signature(mesh)
        );
    }
}

fn main() {
    let path = std::path::PathBuf::from(std::env::args_os().nth(1).expect("STEP input"));
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(move || run(&path))
        .unwrap()
        .join()
        .unwrap();
}
