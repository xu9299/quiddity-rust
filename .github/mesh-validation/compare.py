"""Compare the emitted effective-surface bytes with the committed snapshot and the upstream control."""

import gzip
import hashlib
import json
from pathlib import Path
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar_changes(left, right, path=()):
    if isinstance(left, dict) and isinstance(right, dict) and left.keys() == right.keys():
        for key in left:
            yield from scalar_changes(left[key], right[key], (*path, key))
    elif isinstance(left, list) and isinstance(right, list) and len(left) == len(right):
        for i, (a, b) in enumerate(zip(left, right)):
            yield from scalar_changes(a, b, (*path, i))
    elif type(left) != type(right) or repr(left) != repr(right):
        yield {"path": path, "expected": left, "actual": right}


def main():
    label, commit, tree = sys.argv[1:]
    emitted = Path("../emitted") / label / "effective_surfaces"
    committed = Path("tests/fixtures/rust-snapshot/effective_surfaces")
    out = Path("../evidence") / label
    observed = {p.relative_to(emitted).as_posix(): digest(p) for p in sorted(emitted.rglob("*.json.gz"))}
    official = {p.relative_to(committed).as_posix(): digest(p) for p in sorted(committed.rglob("*.json.gz"))}
    assert observed and observed.keys() == official.keys()
    changes = {}
    for name, signature in observed.items():
        if signature == official[name]:
            continue
        with gzip.open(committed/name, "rt") as f:
            old = json.load(f)
        with gzip.open(emitted/name, "rt") as f:
            new = json.load(f)
        changes[name] = list(scalar_changes(old, new))
    result = dict(source=commit, tree=tree, total=len(observed), different_from_committed=changes, sha256=observed)
    if label != "main":
        base = json.loads(Path("../evidence/main/comparison.json").read_text())["sha256"]
        assert observed.keys() == base.keys()
        result["different_from_main"] = {name:sha for name,sha in observed.items() if sha != base[name]}
    (out/"comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    print(label, len(observed), "changed from snapshot", len(changes),
          "changed from upstream", len(result.get("different_from_main",{})), flush=True)


if __name__ == "__main__":
    main()
