//! Single-fixture probe for the three motion-sensitive section comparisons.
#[path = "../src/engine.rs"]
mod engine;
pub use engine::*;
extern crate self as quiddity;
use features::Context;
use features::effective_surfaces::EffectiveFaces;
use kernel::step::{IDENTITY, read_step_file_placed};
use serde_json::json;
use std::path::Path;

fn main() {
    let [file, axis]: [String; 2] = std::env::args().skip(1).collect::<Vec<_>>().try_into().unwrap();
    let placed = match axis.as_str() {
        "x" => [[1.,0.,0.,0.],[0.,0.,-1.,0.],[0.,1.,0.,0.]],
        "z" => [[0.,-1.,0.,0.],[1.,0.,0.,0.],[0.,0.,1.,0.]],
        _ => panic!("unknown placement"),
    };
    let mut measurements = Vec::new();
    for (name, placement) in [("identity", IDENTITY), ("quarter_turn", placed)] {
        eprintln!("PLACEMENT {name}");
        let part = read_step_file_placed(Path::new(&file), &placement).unwrap();
        let ctx = Context::new(&part);
        let surfaces = EffectiveFaces::new(&ctx);
        let pockets: Vec<_> = features::cylindrical_pockets::cylindrical_pocket_proofs(&ctx, &surfaces)
            .iter().map(|p| json!({"floor":p.floor,"walls":p.walls,"stock":p.stock,"volume":p.volume})).collect();
        let candidates = features::section_recess_geometry::candidates(&ctx, &surfaces).unwrap();
        let faces: Vec<_> = candidates.iter().map(|c| &c.constituent_faces).collect();
        let found = features::recognise(&part);
        let sections: Vec<_> = found.defining.iter()
            .filter(|(family,_)| ["section_passages","section_recesses"].contains(family))
            .map(|(family, sets)| (family, sets)).collect();
        let areas: Vec<_> = (0..part.faces.len()).map(|f| part.face_mass(f).map(|v| v[0])).collect();
        measurements.push(json!({"placement":name,"pocket_proofs":pockets,"candidate_face_sets":faces,"physical_section_families":sections,"areas":areas}));
    }
    println!("{}", json!(measurements));
}
