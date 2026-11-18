
// Read-only route census compiled as an extra test on the validation branch.
#[cfg(test)]
#[test]
fn survey_master_boundary_routes() {
    use serde_json::{Value, json};
    use std::{fs, path::PathBuf};

    fn inspect(path: &std::path::Path) -> Value {
        let part = crate::step::read_step_file(path).unwrap();
        let mut flagged = 0;
        let mut native = 0;
        let mut usable_loops = 0;
        let mut master = 0;
        let mut projected_also_valid = 0;
        let mut same_area_bits = 0;
        for (i, face) in part.faces.iter().enumerate() {
            if part.primary_pcurves(i).next().is_none() {
                continue;
            }
            flagged += 1;
            if !matches!(&face.surface, Surface::Freeform { kind, .. } if matches!(*kind, "BSPLINE" | "BEZIER")) {
                continue;
            }
            native += 1;
            if part.integration_pcurves(i).is_some() {
                usable_loops += 1;
            }
            let area = part.face_mass(i);
            if !part.cache[i].mass.get().unwrap().preferred {
                continue;
            }
            master += 1;
            if let Some(projected) = part.projected_face_mass(i).terms {
                projected_also_valid += 1;
                if projected[0].to_bits() == area.unwrap()[0].to_bits() {
                    same_area_bits += 1;
                }
            }
        }
        json!({"faces":part.faces.len(),"flagged":flagged,"native":native,
            "usable_loops_before_quadrature":usable_loops,"master":master,
            "valid_projection_among_master":projected_also_valid,
            "identical_area_among_both":same_area_bits})
    }

    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let control = inspect(&source.join("tests/fixtures/master-trim-square.step"));
    assert_eq!(control["master"], 1);
    assert_eq!(control["valid_projection_among_master"], 1);
    assert_eq!(control["identical_area_among_both"], 0);
    let manifest: Value = serde_json::from_slice(
        &fs::read(source.join("../../tests/fixtures/corpus.json")).unwrap()
    ).unwrap();
    let corpus = PathBuf::from(std::env::var_os("QUIDDITY_CORPUS").unwrap());
    let mut rows = std::collections::BTreeMap::new();
    for row in manifest["files"].as_array().unwrap() {
        let name = row["file"].as_str().unwrap();
        assert!(rows.insert(name, inspect(&corpus.join(name))).is_none());
    }
    let report = json!({"control":control,"files":rows});
    fs::write(std::env::var_os("MASTER_TRIM_ROUTES").unwrap(),
        serde_json::to_vec_pretty(&report).unwrap()).unwrap();
}
