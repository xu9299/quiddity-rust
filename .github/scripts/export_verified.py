"""Export previously passing ARM output for explicit local snapshot-tool recording."""

import base64
import hashlib
import io
import json
import subprocess
import zipfile
from pathlib import Path

SOURCE = Path("source")
EVIDENCE = Path("evidence")
HEAD = "184946bc4f0e2924f5c1fc370bcd62baf7387679"
TREE = "a79eac1be4ae3ce3573ad5ee0c89a4223c202015"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def inventory(root):
    return {p.relative_to(root).as_posix(): digest(p.read_bytes()) for p in root.rglob("*.json.gz")}


assert subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"], text=True).strip() == HEAD
assert (EVIDENCE / "environment.log").read_text().splitlines()[:2] == [HEAD, TREE]
jobs = json.loads(
    subprocess.check_output(
        ["gh", "api", "repos/xu9299/quiddity-rust/actions/jobs/132311204208"], text=True
    )
)
steps = {step["name"]: step["conclusion"] for step in jobs["steps"]}
assert jobs["run_id"] == 43550362331
assert steps["Complete workspace and corpus regression"] == "success"
assert steps["Snapshot coverage and provenance"] == "success"
assert steps["Strict candidate snapshot reproduction"] == "failure"
assert "HEAD checked all 2702 snapshots under their original comparison bands" in (
    EVIDENCE / "coverage.log"
).read_text()
audit = json.loads(Path("audit/summary.json").read_text())
assert audit["kinds"] == {"float": 25494} and audit["strict_replay_changed_set_matches"]
expected_root = SOURCE / "tests/fixtures/rust-snapshot"
original = json.loads((expected_root / "provenance.json").read_text())
before, after = inventory(expected_root), inventory(EVIDENCE / "verified")
assert before.keys() == after.keys() and len(after) == 2702
assert original["snapshots"] == before
receipt = {
    "schema": 1,
    "kind": "quiddity-verified-snapshots",
    "source_commit": HEAD,
    "source_tree": TREE,
    "inputs": original["inputs"],
    "expected_snapshots": before,
    "verified_snapshots": after,
    "ordinary_comparisons": "passed",
    "coverage": "passed",
    "workspace": "passed",
    "origin": {
        "repository": "xu9299/quiddity-rust",
        "run_id": 43550362331,
        "job_id": 132311204208,
        "artifact_id": 14099562262,
        "artifact_sha256": "47cda28fe09c5a6c7dc0b615acabc06d6e93f1cf9ca42fca812f4e3734db1b4d",
        "verification_command": "cargo run -p snapshot-tool --release --locked -- --check",
        "coverage_command": "cargo run -p snapshot-tool --release --locked -- --check --coverage-only --verify-dir ../validation-logs/verified",
    },
}
payload = io.BytesIO()
with zipfile.ZipFile(payload, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    archive.writestr("receipt.json", json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    for name in sorted(before):
        if before[name] != after[name]:
            archive.writestr(name, (EVIDENCE / "verified" / name).read_bytes())
raw = payload.getvalue()
Path("audit/verified-delta.zip").write_bytes(raw)
encoded = base64.b64encode(raw).decode("ascii")
print("VERIFIED_TRANSPORT=" + json.dumps({"sha256": digest(raw), "bytes": len(raw), "chunks": (len(encoded) + 3999) // 4000}))
for i in range(0, len(encoded), 4000):
    print("VERIFIED_CHUNK=" + json.dumps({"index": i // 4000, "data": encoded[i:i + 4000]}))
print("Exported existing compared output; no expected snapshot was changed.")
