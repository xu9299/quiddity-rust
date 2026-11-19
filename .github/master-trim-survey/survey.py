"""Compare raw areas, every published fingerprint and the complete recognition envelope."""
import argparse
import collections
import concurrent.futures
import hashlib
import json
import math
from pathlib import Path
import subprocess


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def equal(a, b):
    return canonical(a) == canonical(b)


def changed_keys(a, b):
    return sorted(k for k in a.keys() | b.keys() if k not in a or k not in b or not equal(a[k], b[k]))


def changed_rows(a, b, key):
    left = {r[key]: r for r in a}
    right = {r[key]: r for r in b}
    assert len(a) == len(left) and len(b) == len(right)
    return changed_keys(left, right)


def compare(name, pair):
    old, new = pair
    if "read_error" in old or "read_error" in new:
        return {"file": name, "read_errors": [old.get("read_error"), new.get("read_error")]}
    a, b = old["faces"], new["faces"]
    assert [(f["index"], f["entity"]) for f in a] == [(f["index"], f["entity"]) for f in b], name
    changed, gained, lost = [], [], []
    for x, y in zip(a, b, strict=True):
        before, after = x["area"], y["area"]
        ref = {"index": x["index"], "entity": x["entity"], "before": before, "after": after}
        for value in (before, after):
            assert value is None or (value > 0 and math.isfinite(value)), (name, ref)
        if before is None and after is not None:
            gained.append(ref)
        elif before is not None and after is None:
            lost.append(ref)
        elif before != after:
            changed.append({**ref, "relative": abs(after - before) / before})
    a, b = old["recognition"].copy(), new["recognition"].copy()
    fa, fb = a.pop("fingerprints"), b.pop("fingerprints")
    da, db = a.pop("document"), b.pop("document")
    assert len(fa["faces"]) == len(old["faces"]) and len(fb["faces"]) == len(new["faces"])
    face_fp = changed_rows(fa["faces"], fb["faces"], "index")
    feat_fp = changed_rows(fa["features"], fb["features"], "id")
    return dict(file=name, faces=len(old["faces"]),
                measured=[sum(f["area"] is not None for f in side["faces"]) for side in pair],
                areas=changed, gained=gained, lost=lost,
                fingerprint_fields=changed_keys(fa, fb), face_fingerprints=face_fp,
                feature_fingerprints=feat_fp, document_fields=changed_keys(da, db),
                record_ids=changed_rows(da["records"], db["records"], "id"),
                record_counts=[len(da["records"]), len(db["records"])],
                published_families=changed_keys(a, b))


def main():
    parser = argparse.ArgumentParser()
    for key in ("manifest", "corpus", "baseline", "candidate", "control", "out"):
        parser.add_argument("--" + key, required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    entries = json.loads(args.manifest.read_bytes())["files"]
    names = {x["file"] for x in entries}
    assert len(names) == len(entries)
    on_disk = {p.relative_to(args.corpus).as_posix() for p in args.corpus.rglob("*")
               if p.is_file() and p.name.lower().endswith((".step", ".stp", ".step.gz", ".stp.gz"))}
    assert names == on_disk, (sorted(names - on_disk), sorted(on_disk - names))
    def capture(binary, path, digest):
        cp = subprocess.run([str(binary.resolve()), str(path), digest], check=True, stdout=subprocess.PIPE)
        return json.loads(cp.stdout)
    control_digest = hashlib.sha256(args.control.read_bytes()).hexdigest()
    witness = [capture(binary, args.control, control_digest)["faces"][0]["area"]
               for binary in (args.baseline, args.candidate)]
    assert abs(witness[0] - 1) < 1e-12 and abs(witness[1] - (1 - 0.0004 / 3)) < 1e-12, witness
    def trial(entry):
        name = entry["file"]
        path = args.corpus / name
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == entry["sha256"], name
        results = [capture(binary, path, digest) for binary in (args.baseline, args.candidate)]
        row = compare(name, results)
        print(name, row.get("faces"), len(row.get("areas", [])), flush=True)
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows = sorted(pool.map(trial, entries), key=lambda r: r["file"])
    assessed = [r for r in rows if "read_errors" not in r]
    changes = [{"file": r["file"], **v} for r in assessed for v in r["areas"]]
    def affected(key):
        return sum(bool(r[key]) for r in assessed)
    report = dict(inputs=len(rows), read_failures=[r for r in rows if "read_errors" in r],
                  faces=sum(r["faces"] for r in assessed),
                  measured=[sum(r["measured"][i] for r in assessed) for i in range(2)],
                  changed_previously_measured=len(changes), max_relative_change=max(changes, key=lambda c:c["relative"], default=None),
                  gained=sum(len(r["gained"]) for r in assessed), lost=sum(len(r["lost"]) for r in assessed),
                  area_changed_files=affected("areas"), fingerprint_changed_files=affected("fingerprint_fields"),
                  face_fingerprints_changed=sum(len(r["face_fingerprints"]) for r in assessed),
                  feature_fingerprints_changed=sum(len(r["feature_fingerprints"]) for r in assessed),
                  document_changed_files=affected("document_fields"), record_changed_files=affected("record_ids"),
                  records_changed=sum(len(r["record_ids"]) for r in assessed), published_family_changed_files=affected("published_families"),
                  record_counts=[sum(r["record_counts"][i] for r in assessed) for i in range(2)],
                  counts_by_document_field=dict(collections.Counter(k for r in assessed for k in r["document_fields"])),
                  control_areas=witness)
    (args.out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.out / "per-file.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
