"""Run the complement of the official PR filters using the upstream distribution logic."""
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

source = Path.cwd()
work = source.parent
inputs = work / "workload"
plan = inputs / "tier-plan"
shard = int(os.environ["SHARD"])
sys.path.insert(0, str(source / ".github"))
from shards import assign, filterset, listed

def inventory(name):
    lines = (plan / (name + ".txt")).read_text().splitlines()
    cases = [tuple(line.split("\t", 1)[0].split(" ", 1)) for line in lines if line and not line.startswith("#")]
    assert len(cases) == len(set(cases)), name
    return set(cases)

full = inventory("full")
pr = inventory("pr")
remaining = inventory("main-only")
assert full == pr | remaining and pr.isdisjoint(remaining)
assert (len(full), len(pr), len(remaining)) == (1629, 1497, 132)
weights = json.loads((source / ".config/test-times.json").read_text())["Linux"]
chunks, loads = assign(remaining, weights, 4)
expected = set(chunks[shard-1])
assert expected and set().union(*map(set,chunks)) == remaining
expr = filterset(expected)
common = ["--archive-file",str(inputs / "tests-ubuntu-latest.tar.zst"),"--extract-to",str(source),"--extract-overwrite","--workspace-remap",str(source),"--profile","ci","-E",expr]
selected = subprocess.check_output(["cargo-nextest","nextest","list",*common,"--message-format","json"])
assert set(listed(json.loads(selected))) == expected
dest = work / "receipt"
dest.mkdir(exist_ok=True)
print(f"Complement {shard} of 4: {len(expected)} main-only cases, about {loads[shard-1]:.0f} s", flush=True)
command = ["cargo-nextest","nextest","run",*common,"--status-level","pass"]
with subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1) as proc, (dest/"nextest.log").open("w") as out:
    for line in proc.stdout:
        print(line, end="", flush=True)
        out.write(line)
    result = proc.wait()
if result:
    raise SystemExit(result)
log = re.sub(r"\x1b\[[0-9;]*[a-zA-Z]","",(dest/"nextest.log").read_text())
passed = re.findall(r"^\s*(?:PASS|FLAKY)\s+\[[^\]]+\]\s+(?:\(\s*\d+\s*/\s*\d+\)\s+)?(\S+)\s+(\S+)\s*$",log,re.M)
summary = re.search(r"Summary\s+\[[^\]]+\]\s+(\d+) tests? run:\s*(\d+) passed",log)
assert summary and int(summary[1]) == int(summary[2]) == len(expected)
assert set(passed) == expected and len(passed) == len(expected), (len(expected),len(passed))
receipt = {"shard":shard,"source":os.environ["MERGE"],"tree":os.environ["TREE"],"corpus":os.environ["CORPUS"],
           "origin":(inputs/"origin.txt").read_text().strip(),"digests":(inputs/"SHA256SUMS").read_text().splitlines(),
           "upstream_artifact":15078986237,"upstream_zip_sha256":os.environ["DIGEST"],"pr_selected":len(pr),"full_inventory":len(full),
           "complement":len(remaining),"passed":[" ".join(t) for t in sorted(expected)]}
(dest/"verified.json").write_text(json.dumps(receipt,indent=2)+"\n")
