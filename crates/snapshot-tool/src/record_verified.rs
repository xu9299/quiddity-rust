//! Explicit local recording of reviewed, already-compared remote producer output.
//! This is an intentional expectation change, never a check or a CI acceptance mode.

use super::*;

fn metadata_only(path: &str) -> bool {
    path == "CHANGELOG.md"
        || path == "CONTRIBUTING.md"
        || path.starts_with("docs/")
        || path.starts_with("crates/snapshot-tool/")
        || path.starts_with("tests/fixtures/rust-snapshot/")
}

fn validate_receipt(
    receipt: &Value,
    before_inputs: &Digests,
    expected: &Digests,
    produced: &Digests,
) -> Result<()> {
    if receipt["schema"] != 1 || receipt["kind"] != "quiddity-verified-snapshots" {
        return Err("unsupported verified snapshot receipt".into());
    }
    for field in ["source_commit", "source_tree"] {
        if !receipt[field]
            .as_str()
            .is_some_and(|v| v.len() == 40 && v.bytes().all(|b| b.is_ascii_hexdigit()))
        {
            return Err(format!("verified receipt requires a full {field}").into());
        }
    }
    for field in ["ordinary_comparisons", "coverage", "workspace"] {
        if receipt[field] != "passed" {
            return Err(format!("verified receipt does not establish {field}").into());
        }
    }
    if receipt["origin"].as_object().is_none_or(|o| o.is_empty()) {
        return Err("verified receipt requires producer origin information".into());
    }
    compare_inventory(produced, expected, true)?;
    if receipt["inputs"] != json!(before_inputs)
        || receipt["expected_snapshots"] != json!(expected)
        || receipt["verified_snapshots"] != json!(produced)
    {
        return Err(
            "verified receipt input, previous expectation, or output digests differ".into(),
        );
    }
    Ok(())
}

pub(super) fn run(
    root: &Path,
    dir: &Path,
    corpus: &Path,
    before_inputs: &Digests,
    head: &str,
) -> Result<()> {
    // Keep this invariant here as well as in the CLI dispatcher.
    if std::env::var_os("CI").is_some() {
        return Err("snapshot recording is forbidden in CI".into());
    }
    let receipt: Value = serde_json::from_slice(&fs::read(dir.join("receipt.json"))?)?;
    let expected = snapshot_digests(&root.join(PREFIX))?;
    let produced = snapshot_digests(dir)?;
    validate_receipt(&receipt, before_inputs, &expected, &produced)?;
    check(root, dir, before_inputs, true)?;
    let source = receipt["source_commit"].as_str().unwrap();
    if !output(root, &["merge-base", "--is-ancestor", source, "HEAD"])?
        .status
        .success()
        || receipt["source_tree"] != git(root, &["rev-parse", &format!("{source}^{{tree}}")])?
    {
        return Err("verified producer source is not an identified ancestor of HEAD".into());
    }
    let changed = git(root, &["diff", "--name-only", source, "HEAD"])?;
    if changed.lines().any(|path| !metadata_only(path)) {
        return Err(
            "recognition, test, dependency or input source changed since verification".into(),
        );
    }
    let staging = Temporary::new()?;
    for name in produced.keys() {
        let target = staging.0.join(name);
        fs::create_dir_all(target.parent().unwrap())?;
        fs::copy(dir.join(name), target)?;
    }
    if snapshot_digests(&staging.0)? != produced || snapshot_digests(dir)? != produced {
        return Err("verified output changed during recording".into());
    }
    clean(root)?;
    if git(root, &["rev-parse", "HEAD"])? != head
        || snapshot_digests(&root.join(PREFIX))? != expected
        || inputs(root, corpus)? != *before_inputs
    {
        return Err("checkout or expectations changed during recording".into());
    }
    let provenance = json!({
        "schema": 1,
        "source_commit": source,
        "command": format!("{INVOCATION} --record-verified DIR"),
        "inputs": before_inputs,
        "snapshots": produced,
        "recorded_from": receipt["origin"],
    });
    fs::write(
        staging.0.join("provenance.json"),
        format!("{}\n", serde_json::to_string_pretty(&provenance)?),
    )?;
    install(root, &staging.0)?;
    println!(
        "Recorded {} reviewed snapshots from Rust {source}; review the diff, commit, and run strict --check",
        produced.len()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> (Value, Digests, Digests, Digests) {
        let inputs = Digests::from([("fixture.step".into(), "input".into())]);
        let expected = Digests::from([
            ("corpus.json.gz".into(), "old-corpus".into()),
            ("pmi.json.gz".into(), "same-pmi".into()),
        ]);
        let mut produced = expected.clone();
        produced.insert("corpus.json.gz".into(), "reviewed-corpus".into());
        let receipt = json!({
            "schema": 1, "kind": "quiddity-verified-snapshots",
            "source_commit": "1".repeat(40), "source_tree": "2".repeat(40),
            "ordinary_comparisons": "passed", "coverage": "passed", "workspace": "passed",
            "origin": {"job": "reviewed producer"},
            "inputs": inputs, "expected_snapshots": expected, "verified_snapshots": produced,
        });
        (receipt, inputs, expected, produced)
    }

    #[test]
    fn recording_requires_complete_original_and_produced_digest_inventories() {
        let (receipt, inputs, expected, produced) = fixture();
        assert!(validate_receipt(&receipt, &inputs, &expected, &produced).is_ok());
        for field in ["inputs", "expected_snapshots", "verified_snapshots"] {
            let mut damaged = receipt.clone();
            damaged[field] = json!({});
            assert!(validate_receipt(&damaged, &inputs, &expected, &produced).is_err());
        }
        let mut orphan = produced.clone();
        orphan.insert("extra.json.gz".into(), "extra".into());
        let mut damaged = receipt;
        damaged["verified_snapshots"] = json!(orphan);
        assert!(validate_receipt(&damaged, &inputs, &expected, &orphan).is_err());
    }

    #[test]
    fn recording_refuses_failed_checks_or_missing_source_identity() {
        let (receipt, inputs, expected, produced) = fixture();
        for field in [
            "ordinary_comparisons",
            "coverage",
            "workspace",
            "source_commit",
            "source_tree",
            "origin",
        ] {
            let mut damaged = receipt.clone();
            damaged[field] = Value::Null;
            assert!(validate_receipt(&damaged, &inputs, &expected, &produced).is_err());
        }
        assert!(metadata_only("crates/snapshot-tool/src/main.rs"));
        for path in [
            "src/recognition.rs",
            "crates/haecceity/src/mass.rs",
            "Cargo.lock",
            "tests/internal/mass.rs",
            "tests/fixtures/input.step",
        ] {
            assert!(!metadata_only(path));
        }
    }

    #[test]
    fn recording_is_explicit_and_cannot_be_a_check_mode() {
        for flags in [
            vec!["--record-verified"],
            vec!["--record-verified", "dir", "--check"],
            vec!["--record-verified", "dir", "--diff", "HEAD"],
            vec!["--record-verified", "dir", "--coverage-only"],
        ] {
            assert!(Args::parse(flags.into_iter().map(str::to_owned)).is_err());
        }
        assert!(Args::parse(["--record-verified", "dir"].into_iter().map(str::to_owned)).is_ok());
    }

    #[test]
    fn ci_recording_refuses_before_any_file_access() {
        if std::env::var_os("QUIDDITY_RECORD_GUARD_CHILD").is_some() {
            let missing = Path::new("nonexistent-verified-snapshot-test");
            let error = run(missing, missing, missing, &Digests::new(), "").unwrap_err();
            assert_eq!(error.to_string(), "snapshot recording is forbidden in CI");
        } else {
            assert!(
                Command::new(std::env::current_exe().unwrap())
                    .args([
                        "record_verified::tests::ci_recording_refuses_before_any_file_access",
                        "--exact"
                    ])
                    .env("QUIDDITY_RECORD_GUARD_CHILD", "1")
                    .env("CI", "1")
                    .status()
                    .unwrap()
                    .success()
            );
        }
    }
}
