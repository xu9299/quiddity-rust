//! Research-only process probe linked to an unchanged haecceity library.
use std::hint::black_box;
use std::time::Instant;

fn measure(path: &str, calls: usize) {
    let fd = unsafe { mesh_pmu_begin() };
    assert!(fd >= 0);
    let start = Instant::now();
    let part = haecceity::read_step_file(std::path::Path::new(path)).unwrap();
    let read_seconds = start.elapsed().as_secs_f64();
    let diagonal = part.bounds().diagonal();
    let deflection = if diagonal > 0. && diagonal.is_finite() {
        0.001 * diagonal
    } else {
        0.001
    };
    let mut mesh_seconds = Vec::with_capacity(calls);
    let mut expected = None;
    for _ in 0..calls {
        let start = Instant::now();
        let meshes = part.mesh(deflection, 0.3);
        mesh_seconds.push(start.elapsed().as_secs_f64());
        let mut counts = [meshes.len(), 0usize, 0usize, 0usize];
        for mesh in &meshes {
            match mesh {
                Ok(mesh) => {
                    counts[1] += mesh.points.len();
                    counts[2] += mesh.triangles.len();
                }
                Err(_) => counts[3] += 1,
            }
        }
        assert_eq!(*expected.get_or_insert(counts), counts);
        black_box(&meshes);
        drop(meshes);
    }
    drop(part);
    assert_eq!(unsafe { mesh_pmu_end(fd) }, 0);
    let [faces, points, triangles, refusals] = expected.unwrap();
    println!(
        "{{\"calls\":{calls},\"faces\":{faces},\"points\":{points},\"triangles\":{triangles},\"refusals\":{refusals},\"read_seconds\":{read_seconds},\"mesh_seconds\":{mesh_seconds:?}}}"
    );
}

fn main() {
    let path = std::env::args().nth(1).expect("input STEP path");
    let calls: usize = std::env::args()
        .nth(2)
        .expect("mesh call count")
        .parse()
        .unwrap();
    assert!(calls == 1 || calls == 10);
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(move || measure(&path, calls))
        .unwrap()
        .join()
        .unwrap();
}

unsafe extern "C" { fn mesh_pmu_begin() -> i32; fn mesh_pmu_end(fd: i32) -> i32; }
