//! Read the committed probe inputs and compare every result with the normal snapshot predicate.
use haecceity::{
    classify::Classifier,
    geom::Bounds,
    rays::RayCaster,
    read_step, read_step_file,
    volume::{Probe, common_volume, probe_volume},
};
use serde_json::{Value, json};
use std::{collections::BTreeMap, env, fs, path::PathBuf};
fn main() {
    let args: Vec<PathBuf> = env::args_os().skip(1).map(PathBuf::from).collect();
    let inputs: Value = serde_json::from_slice(&fs::read(&args[0]).unwrap()).unwrap();
    let snapshots: Vec<Value> = serde_json::from_slice(&fs::read(&args[1]).unwrap()).unwrap();
    let mut asked = BTreeMap::<&str, Vec<&Value>>::new();
    for p in inputs["probes"].as_array().unwrap() {
        asked
            .entry(p["file"].as_str().unwrap())
            .or_default()
            .push(p);
    }
    assert_eq!(asked.len(), snapshots.len());
    let mut examined = 0;
    let mut changed = 0;
    let mut failed = Vec::new();
    for expected in &snapshots {
        let file = expected["file"].as_str().unwrap();
        let part = read_step_file(&args[2].join(file)).unwrap();
        let classifiers: Vec<_> = (0..part.solids.len())
            .map(|s| Classifier::for_solid(&part, s))
            .collect();
        let cases = &asked[file];
        let reference = expected["probes"].as_array().unwrap();
        assert_eq!(cases.len(), reference.len());
        for (index, (p, r)) in cases.iter().zip(reference).enumerate() {
            let other;
            let probe = if let Some(b) = p.get("box") {
                let c: [f64; 6] = serde_json::from_value(b.clone()).unwrap();
                Probe::Box(Bounds {
                    min: [c[0], c[1], c[2]],
                    max: [c[3], c[4], c[5]],
                })
            } else {
                other = read_step(p["step"].as_str().unwrap().as_bytes()).unwrap();
                Probe::Solid(RayCaster::for_solid(&other, 0))
            };
            let actual: Option<f64> = p["solids"]
                .as_array()
                .unwrap()
                .iter()
                .map(|s| common_volume(&classifiers[s.as_u64().unwrap() as usize], &probe))
                .sum();
            let pv = probe_volume(&probe);
            let snapshot = r["volume"].as_f64();
            let sv = r["probe_volume"].as_f64();
            let same = match (actual, pv, snapshot, sv) {
                (Some(g), Some(_), Some(0.0), Some(_)) => g == 0.0,
                (Some(g), Some(gv), Some(w), Some(wv)) if (w - wv).abs() <= 1e-9 * wv => {
                    (g - gv).abs() <= 1e-9 * gv
                }
                (Some(g), Some(gv), Some(w), Some(wv)) => (g / gv - w / wv).abs() <= 1e-6,
                _ => actual == snapshot && pv == sv,
            };
            if !same {
                failed.push(json!({"file":file,"index":index,"caller":p["caller"],"actual":actual,"probe_volume":pv,"snapshot":snapshot,"snapshot_probe":sv,"occt":p["volume"],"occt_probe":p["probe_volume"]}));
            }
            examined += 1;
            changed += usize::from(actual != snapshot || pv != sv);
        }
    }
    let result = json!({"files":asked.len(),"probes":examined,"raw_values_changed":changed,"ordinary_failures":failed});
    fs::write(&args[3], serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    println!(
        "{examined} probes, {changed} raw changes, {} beyond the ordinary snapshot comparator",
        failed.len()
    );
}
