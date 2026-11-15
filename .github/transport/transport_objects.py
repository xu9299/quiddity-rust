"""Transport exact, already-recorded Git objects; never produce or record snapshots in CI."""

import base64
import hashlib
import json
import subprocess
from pathlib import Path

spec = json.loads(Path("maintenance/.github/transport/objects.json").read_text())
source = Path("source")
parent = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
assert parent == spec["parent"]
assert spec["commit"] == "c294e9b1f5ab93b68b8c26a319ccf299889ed000"
for item in spec["objects"]:
    if "artifact_path" in item:
        name = item["artifact_path"]
        assert ".." not in Path(name).parts and not Path(name).is_absolute()
        payload = (Path("evidence/verified") / name).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == item["sha256"]
    else:
        payload = base64.b64decode(item["base64"], validate=True)
    actual = (
        subprocess.check_output(
            ["git", "hash-object", "-w", "-t", item["type"], "--stdin"],
            cwd=source,
            input=payload,
        )
        .decode()
        .strip()
    )
    assert actual == item["oid"], (actual, item["oid"])
    print(f"Verified existing {item['type']} object {actual}")
assert (
    subprocess.check_output(
        ["git", "rev-parse", spec["commit"] + "^{tree}"],
        cwd=source,
        text=True,
    ).strip()
    == spec["tree"]
)
assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=source).strip()
subprocess.run(
    ["git", "merge-base", "--is-ancestor", parent, spec["commit"]],
    cwd=source,
    check=True,
)
subprocess.run(
    ["git", "push", "origin", spec["commit"] + ":refs/heads/perf/nurbs-isocurve-ready-052"],
    cwd=source,
    check=True,
)
print("Transported the exact local commit; no expectation was generated or recorded here.")
