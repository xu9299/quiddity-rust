//! Copy into crates/haecceity/tests/face_reproducer.rs in the pinned upstream checkout.
//! Supply FACE_REPRO_INPUT and FACE_REPRO_OPERATION=mass or moments. No input is rewritten.

use std::path::Path;
use std::time::Instant;

#[test]
#[ignore = "requires an explicitly selected extracted STEP sample"]
fn measure_extracted_faces() {
    let input = std::env::var("FACE_REPRO_INPUT").expect("set FACE_REPRO_INPUT");
    let operation = std::env::var("FACE_REPRO_OPERATION").unwrap_or_else(|_| "mass".to_owned());
    assert!(matches!(operation.as_str(), "mass" | "moments"));
    let started = Instant::now();
    let part = haecceity::step::read_step_file(Path::new(&input)).expect("read STEP sample");
    eprintln!(
        "read_seconds={} faces={}",
        started.elapsed().as_secs_f64(),
        part.faces.len()
    );
    assert!(!part.faces.is_empty());
    for face in 0..part.faces.len() {
        let source = part.face_source(face).map(|s| s.entity);
        eprintln!("starting face={face} source={source:?} operation={operation}");
        let started = Instant::now();
        if operation == "mass" {
            eprintln!("result={:?}", part.face_mass(face));
        } else {
            eprintln!("result={:?}", part.face_moments(face));
        }
        eprintln!("seconds={}", started.elapsed().as_secs_f64());
    }
}
