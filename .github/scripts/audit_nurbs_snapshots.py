"""Read-only comparison of committed and already-verified snapshot artifacts."""

import argparse
import gzip
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

SOURCE = "184946bc4f0e2924f5c1fc370bcd62baf7387679"
TREE = "a79eac1be4ae3ce3573ad5ee0c89a4223c202015"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def child(path, key):
    return path + "/" + str(key).replace("~", "~0").replace("/", "~1")


def differences(before, after, path=""):
    numeric = (int, float)
    if (
        isinstance(before, numeric)
        and not isinstance(before, bool)
        and isinstance(after, numeric)
        and not isinstance(after, bool)
    ):
        if not (math.isfinite(before) and math.isfinite(after)):
            if before != after:
                yield dict(kind="nonfinite", path=path, before=before, after=after)
        elif before != after:
            delta = abs(after - before)
            yield dict(
                kind="integer" if type(before) is int and type(after) is int else "float",
                path=path,
                before=before,
                after=after,
                absolute=delta,
                normalized=delta / (1 + abs(before)),
            )
        elif type(before) is not type(after):
            yield dict(kind="number_type", path=path, before=before, after=after)
        elif before == 0 and math.copysign(1, before) != math.copysign(1, after):
            yield dict(kind="zero_sign", path=path, before=before, after=after)
        return
    if type(before) is not type(after):
        yield dict(kind="type", path=path, before=type(before).__name__, after=type(after).__name__)
    elif isinstance(before, dict):
        for key in sorted(before.keys() - after.keys()):
            yield dict(kind="key_removed", path=child(path, key))
        for key in sorted(after.keys() - before.keys()):
            yield dict(kind="key_added", path=child(path, key))
        for key in sorted(before.keys() & after.keys()):
            yield from differences(before[key], after[key], child(path, key))
    elif isinstance(before, list):
        if len(before) != len(after):
            yield dict(kind="list_length", path=path, before=len(before), after=len(after))
        for index, (a, b) in enumerate(zip(before, after, strict=False)):
            yield from differences(a, b, child(path, index))
    elif before != after:
        yield dict(kind=type(before).__name__, path=path, before=before, after=after)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    environment = (args.evidence / "environment.log").read_text(encoding="utf-8-sig")
    assert environment.splitlines()[:2] == [SOURCE, TREE], "artifact source identity differs"
    actual_root = args.evidence / "verified"
    expected = {
        p.relative_to(args.expected).as_posix(): p for p in args.expected.rglob("*.json.gz")
    }
    actual = {p.relative_to(actual_root).as_posix(): p for p in actual_root.rglob("*.json.gz")}
    missing, added = (
        sorted(expected.keys() - actual.keys()),
        sorted(actual.keys() - expected.keys()),
    )
    assert not missing and not added and len(actual) == 2702, (missing, added, len(actual))
    strict_log = (args.evidence / "strict.log").read_text(encoding="utf-8-sig")
    match = re.search(
        r"Rust snapshot differs: missing=\[\], added=\[\], changed=(\[.*\])", strict_log
    )
    assert match, "strict failure was not the expected complete-inventory byte comparison"
    strict_changed = json.loads(match.group(1))
    files, all_deltas = [], []
    for name in sorted(actual):
        before, after = expected[name].read_bytes(), actual[name].read_bytes()
        if before == after:
            continue
        plain_before, plain_after = gzip.decompress(before), gzip.decompress(after)
        deltas = [
            dict(file=name, **d)
            for d in differences(json.loads(plain_before), json.loads(plain_after))
        ]
        all_deltas.extend(deltas)
        numbers = [d for d in deltas if d["kind"] == "float"]
        files.append(
            dict(
                file=name,
                expected_sha256=sha(before),
                actual_sha256=sha(after),
                decoded_bytes_equal=plain_before == plain_after,
                kinds=dict(Counter(d["kind"] for d in deltas)),
                max_absolute=max((d["absolute"] for d in numbers), default=0.0),
                max_normalized=max((d["normalized"] for d in numbers), default=0.0),
            )
        )
    kinds = dict(Counter(d["kind"] for d in all_deltas))
    numbers = [d for d in all_deltas if d["kind"] == "float"]
    other = [d for d in all_deltas if d["kind"] not in {"float", "zero_sign"}]
    summary = dict(
        source_commit=SOURCE,
        source_tree=TREE,
        artifact_id=14099562262,
        count=len(actual),
        compressed_changed=len(files),
        strict_replay_changed_count=len(strict_changed),
        strict_replay_changed_set_matches=sorted(strict_changed)
        == sorted(f["file"] for f in files),
        kinds=kinds,
        normalized_formula="abs(after-before)/(1+abs(before))",
        max_absolute=max((d["absolute"] for d in numbers), default=0.0),
        max_normalized=max((d["normalized"] for d in numbers), default=0.0),
        largest_absolute=sorted(numbers, key=lambda d: d["absolute"], reverse=True)[:12],
        largest_normalized=sorted(numbers, key=lambda d: d["normalized"], reverse=True)[:12],
        other_changes=other[:100],
        files=files,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    for name, data in [("summary.json", summary), ("deltas.json", all_deltas)]:
        (args.out / name).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print("SNAPSHOT_AUDIT_JSON=" + json.dumps(summary, separators=(",", ":")))
    print("Read-only audit complete; no expected or produced snapshot was changed.")


if __name__ == "__main__":
    main()
