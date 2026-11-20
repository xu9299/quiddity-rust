//! Validation driver for one corpus input. Build the same driver against both revisions.
use quiddity::kernel::step::read_step_file;
use serde_json::{Value, json};
use std::{error::Error, io, path::Path};

fn run() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let [path, digest] = args.as_slice() else {
        return Err("expected a file and its SHA-256".into());
    };
    let part = match read_step_file(Path::new(path)) {
        Ok(part) => part,
        Err(error) => {
            println!("{}", json!({"read_error": error.to_string()}));
            return Ok(());
        }
    };
    // Exercise recognition first, preserving its normal cold-cache query order.
    let recognition: Value = serde_json::from_str(
        &quiddity::file_api::recognise_part(&part, digest).map_err(io::Error::other)?,
    )?;
    let faces: Vec<_> = (0..part.faces.len())
        .map(|i| {
            let area = part.face_mass(i).map(|m| m[0]);
            json!({"index": i, "entity": part.face_source(i).map(|s| s.entity), "area": area})
        })
        .collect();
    // Preserve the recognition and legacy area survey's normal query order. A
    // complete moment query can fail separately on a very small finite face.
    let moments: Vec<_> = (0..part.faces.len())
        .map(|i| part.face_moments(i).map(|m| (m.area, m.centroid, m.second_moments)))
        .collect();
    println!("{}", json!({"faces": faces, "recognition": recognition, "moments": moments}));
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
