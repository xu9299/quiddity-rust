#[test]
#[ignore = "explicit local comparison against pinned upstream main"]
fn compare_upstream_integral_cost() {
    #[link(name = "kernel32")]
    unsafe extern "system" {
        fn GetCurrentThread() -> *mut std::ffi::c_void;
        fn GetThreadTimes(
            thread: *mut std::ffi::c_void,
            creation: *mut u64,
            exit: *mut u64,
            kernel: *mut u64,
            user: *mut u64,
        ) -> i32;
    }
    let cpu = || {
        let (mut creation, mut exit, mut kernel, mut user) = (0, 0, 0, 0);
        // FILETIME is an eight-byte little-endian integer on this Windows x64 probe.
        unsafe {
            assert_ne!(
                GetThreadTimes(
                    GetCurrentThread(),
                    &mut creation,
                    &mut exit,
                    &mut kernel,
                    &mut user
                ),
                0
            );
        }
        (kernel + user) as f64 * 1e-7
    };
    let names = [
        "held_side_knot.step",
        "late_turn.step",
        "rotated_seam_band.step",
    ];
    let mut parts = Vec::new();
    for name in names {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures")
            .join(name);
        let part = crate::read_step_file(&path).unwrap();
        let faces: Vec<_> = (0..part.faces.len())
            .filter(|&f| matches!(part.faces[f].surface, Surface::Freeform { .. }))
            .collect();
        assert!(!faces.is_empty());
        for &face in &faces {
            part.face_mass(face);
            part.face_moments(face);
        }
        parts.push((name, part, faces));
    }
    let mut rows = Vec::new();
    let mut expected = std::collections::BTreeMap::new();
    for round in 0..5 {
        for offset in 0..parts.len() {
            let (name, part, faces) = &parts[(round + offset) % parts.len()];
            for &face in faces {
                for operation in ["mass", "moments"] {
                    let start_cpu = cpu();
                    let start_wall = std::time::Instant::now();
                    for _ in 0..128 {
                        let bits = if operation == "mass" {
                            part.compute_face_mass(face)
                                .map(|v| v.map(f64::to_bits).to_vec())
                        } else {
                            part.compute_face_moments(face).map(moment_bits)
                        };
                        let key = format!("{name}/{face}/{operation}");
                        if let Some(previous) = expected.get(&key) {
                            assert_eq!(previous, &bits);
                        } else {
                            expected.insert(key, bits);
                        }
                    }
                    rows.push(serde_json::json!({
                        "fixture": name, "face": face, "operation": operation,
                        "round": round, "repeats": 128,
                        "wall_seconds": start_wall.elapsed().as_secs_f64(),
                        "thread_cpu_seconds": cpu() - start_cpu,
                    }));
                }
            }
        }
    }
    let output = std::env::var_os("QUIDDITY_SPLIT_BENCH_OUTPUT").expect("explicit output");
    std::fs::write(output, serde_json::to_vec_pretty(&serde_json::json!({
            "expected": expected, "rows": rows,
            "scope": "Three public trimmed-surface fixtures; warmed geometry caches and direct mass/moment computations; excludes STEP reads, builds and full recognition",
        })).unwrap()).unwrap();
}
