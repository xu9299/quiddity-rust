"""在能提供硬件事件的同一 Linux 运行器上实测一次读取和一次整件网格。"""

import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
from itertools import pairwise
from pathlib import Path

STAGES = ("base", "pr", "scoped", "flow-scoped")
PUBLIC_STAGES = ("base", "flow", "prep", "stack")
HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(value, f, ensure_ascii=False, indent=2)
        f.write("\n")


def orders(stages):
    assert len(stages) == 4
    forward = [[stages[(o + i) % 4] for i in (0, 1, 3, 2)] for o in range(4)]
    return forward + [row[::-1] for row in forward]


def schedule(plan):
    for i, order in enumerate(plan["orders"]):
        cases = plan["cases"]
        split = i % len(cases)
        for name, path in cases[split:] + cases[:split]:
            for j, stage in enumerate(order):
                yield i, j, name, path, stage


def materialize_flow(source, plan):
    def git(*args, data=None, env=None):
        return subprocess.check_output(["git", *args], cwd=source, input=data, env=env)

    # Only previously public sources are needed for the first performance layer.
    patch = git(
        "diff",
        "--binary",
        plan["sources"]["pr"]["source"],
        plan["sources"]["scoped"]["source"],
        "--",
        "crates/haecceity/src",
    )
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(os.environ, GIT_INDEX_FILE=str(Path(tmp) / "index"))
        git("read-tree", plan["flow_parent"], env=env)
        git("apply", "--cached", "-", data=patch, env=env)
        tree = git("write-tree", env=env).decode().strip()
    expected = plan["sources"]["flow-scoped"]
    assert tree == expected["tree"]
    assert (
        git(
            "hash-object", "-w", "-t", "commit", "--stdin", data=(HERE / "flow.commit").read_bytes()
        )
        .decode()
        .strip()
        == expected["source"]
    )


def measure(source, corpus, out):
    import arm_resources_run as common

    plan = json.loads((HERE / "plan.json").read_bytes())
    out.mkdir()
    (out / "plan.json").write_bytes((HERE / "plan.json").read_bytes())
    stages = tuple(plan["sources"])
    assert platform.system() == "Linux" and plan["orders"] == orders(stages)
    assert digest(Path(__file__)) == plan["controller_sha256"]
    assert digest(HERE / "arm_resources_run.py") == plan["build_driver_sha256"]
    assert digest(HERE / "probe.rs") == plan["probe_sha256"]
    assert digest(HERE / "counter.c") == plan["counter_sha256"]
    for name, path in plan["cases"]:
        assert digest(corpus / path) == plan["input_sha256"][name]
    counter = out / "calibrate"
    obj = out / "counter.o"
    cc = ["cc", "-std=gnu17", "-O2", "-Wall", "-Wextra", "-Werror", str(HERE / "counter.c")]
    subprocess.run([*cc, "-DSTUDY_CALIBRATION", "-pthread", "-o", str(counter)], check=True)
    subprocess.run([*cc, "-c", "-o", str(obj)], check=True)
    allowed = sorted(os.sched_getaffinity(0))
    write(
        out / "runner.json",
        {
            "allowed_cpus": allowed,
            "model": subprocess.check_output(["lscpu"], text=True).splitlines(),
            "kernel": platform.release(),
            "paranoid": Path("/proc/sys/kernel/perf_event_paranoid").read_text().strip(),
        },
    )
    trials = []
    prefix = None
    for cpu in [*allowed, None]:
        trial = ["sudo", "-n", "--"] + ([] if cpu is None else ["taskset", "-c", str(cpu)])
        vals = []
        for n in (2000000, 20000000):
            r = subprocess.run([*trial, str(counter), str(n)], capture_output=True)
            vals.append(
                {
                    "iterations": n,
                    "exit": r.returncode,
                    "stdout": r.stdout.decode().strip(),
                    "stderr": r.stderr.decode().strip(),
                }
            )
            if r.returncode:
                break
        trials.append({"cpu": cpu, "values": vals})
        if len(vals) == 2 and all(v["exit"] == 0 for v in vals):
            a, b = (json.loads(v["stderr"]) for v in vals)
            if (
                all(v["instructions"] > 0 and v["running"] / v["enabled"] >= 0.999 for v in (a, b))
                and 5.5 < (b["instructions"] - a["instructions"]) / 18000000 < 8.5
            ):
                prefix = trial
                break
    write(
        out / "capability.json",
        {
            "execution_prefix": prefix,
            "trials": trials,
            "counter_sha256": digest(obj),
            "worker_instructions_confirmed": prefix is not None,
        },
    )
    if prefix is None:
        raise RuntimeError("No usable worker-thread retired-instruction counter on this runner")

    def counted(cmd):
        result = subprocess.run([*prefix, *map(str, cmd)], capture_output=True)
        if result.returncode:
            print(result.stderr.decode(errors="replace"), file=sys.stderr)
        result.check_returncode()
        event = json.loads(result.stderr)
        assert event["instructions"] > 0 and 0.999 <= event["running"] / event["enabled"] <= 1.0
        return result.stdout, event

    if "flow_parent" in plan:
        materialize_flow(source, plan)
    common.STAGES = stages
    binaries = common.build(source, corpus, out, plan, HERE / "probe.rs")
    text = (HERE / "probe.rs").read_text()
    begin = "    let start = Instant::now();"
    end = "    drop(part);"
    assert text.count(end) == 1 and begin in text
    text = text.replace(
        begin, "    let fd = unsafe { mesh_pmu_begin() };\n    assert!(fd >= 0);\n" + begin, 1
    )
    text = text.replace(end, end + "\n    assert_eq!(unsafe { mesh_pmu_end(fd) }, 0);")
    text += '\nunsafe extern "C" { fn mesh_pmu_begin() -> i32; fn mesh_pmu_end(fd: i32) -> i32; }\n'
    worker = out / "worker-probe.rs"
    worker.write_text(text, encoding="utf8", newline="\n")
    measured = {}
    for stage, record in binaries.items():
        cmd = record["rustc_command"].copy()
        exe = out / ("worker-" + stage)
        cmd[cmd.index(str(HERE / "probe.rs"))] = str(worker)
        cmd[1:1] = ["--crate-name", "selected_mesh_worker", "-C", "link-arg=" + str(obj)]
        cmd[-1] = str(exe)
        subprocess.run(cmd, check=True)
        measured[stage] = {
            "path": str(exe),
            "sha256": digest(exe),
            "library_sha256": record["library_sha256"],
        }
    write(
        out / "measured-binaries.json",
        {"stages": measured, "worker_probe_sha256": digest(worker), "object_sha256": digest(obj)},
    )
    expected = {}
    with (out / "samples.jsonl").open("x", encoding="utf8", newline="\n") as f:
        for i, j, name, path, stage in schedule(plan):
            raw, event = counted([measured[stage]["path"], corpus / path, "1"])
            shape = json.loads(raw)
            counts = tuple(shape[k] for k in ("faces", "points", "triangles", "refusals"))
            assert expected.setdefault((name, stage), counts) == counts
            print(
                json.dumps(
                    {"round": i, "position": j, "case": name, "stage": stage, **event, **shape}
                ),
                file=f,
                flush=True,
            )
    write(
        out / "completion.json",
        {
            "samples_sha256": digest(out / "samples.jsonl"),
            "processes": len(list(schedule(plan))),
            "source_trees": plan["sources"],
            "plan_sha256": digest(out / "plan.json"),
        },
    )


def prepare(public_layers=False):
    import cloud
    import reproduce as shared
    import scoped

    cloud.REPO, cloud.OUT, cloud.SOURCE = scoped.REPO, scoped.OUT / "cloud", scoped.SOURCE
    cloud.identity()
    folder = cloud.OUT / ("linux-public-layers" if public_layers else "linux-worker-instructions")
    prereq = shared.read(cloud.OUT / "linux-pmu-self/audit.json")
    assert prereq["counter_usable"]
    corpus = shared.INPUT.parent.parent
    ref = shared.read(scoped.PLAN)
    focus = ["string", "cgb202", "cgb243", "narrow", "cgb245"]
    nist = sorted((corpus / "nist").glob("nist_*.stp"))
    assert len(nist) == 10
    paths = dict(ref["cases"])
    cases = [[name, paths[name]] for name in focus] + [
        [p.stem, p.relative_to(corpus).as_posix()] for p in nist
    ]
    flow = shared.read(scoped.OUT / "layers/flow-scoped.json")
    stages = PUBLIC_STAGES if public_layers else STAGES
    heads = (
        (shared.BASE, flow["parent"], scoped.LAYERS["prep-scoped"], scoped.STACK)
        if public_layers
        else (shared.BASE, shared.PR, scoped.SOURCE, flow["source"])
    )
    build_driver = shared.ROOT / "native/quiddity_mesh_merged_317/arm_resources_run.py"
    probe = build_driver.with_name("arm_resource_probe.rs")
    counter = HERE / "scoped_worker_counter.c"
    plan = {
        "sources": {
            s: {"source": h, "tree": cloud.git("rev-parse", h + "^{tree}").decode().strip()}
            for s, h in zip(stages, heads, strict=True)
        },
        "orders": orders(stages),
        "cases": cases,
        "input_sha256": {n: digest(corpus / p) for n, p in cases},
        "corpus_commit": ref["corpus_commit"],
        "counter_sha256": digest(counter),
        "controller_sha256": digest(Path(__file__)),
        "probe_sha256": digest(probe),
        "build_driver_sha256": digest(build_driver),
        "processes": len(cases) * 8 * 4,
        "rounds": 8,
        "calls_per_fresh_process": 1,
        "scope": (
            "User-mode hardware instructions on the worker thread, immediately before STEP read "
            "through release of the Part after bounds, one public whole mesh and count checks. "
            "Record enabled and running time. Compare paired rounds without scaling."
        ),
    }
    if public_layers:
        for h in heads:
            assert cloud.git(
                "for-each-ref", "--format=%(refname)", "--contains", h, "refs/remotes/origin"
            ).strip()
    else:
        plan["flow_parent"] = flow["parent"]
    write(folder / "plan.json", plan)
    track = "public-layer" if public_layers else "worker"
    branch = f"validation/selected-trim-{track}-instructions-349-20261125"
    workflow = ".github/workflows/selected-trim-instructions-349.yml"
    files = {
        workflow: (folder / "workflow.yml").read_bytes(),
        ".mesh-pmu/counter.c": counter.read_bytes(),
        ".mesh-pmu/probe.rs": probe.read_bytes(),
        ".mesh-pmu/arm_resources_run.py": build_driver.read_bytes(),
        ".mesh-pmu/run.py": Path(__file__).read_bytes(),
        ".mesh-pmu/plan.json": (folder / "plan.json").read_bytes(),
    }
    if not public_layers:
        files[".mesh-pmu/flow.commit"] = cloud.git("cat-file", "commit", flow["source"])
    with tempfile.TemporaryDirectory(dir=shared.ROOT / "development/workspaces") as tmp:
        env = dict(os.environ, GIT_INDEX_FILE=str(Path(tmp) / "index"))
        cloud.git("read-tree", scoped.SOURCE, env=env)
        for path, data in files.items():
            blob = cloud.git("hash-object", "-w", "--stdin", data=data).decode().strip()
            cloud.git("update-index", "--add", "--cacheinfo", f"100644,{blob},{path}", env=env)
        tree = cloud.git("write-tree", env=env).decode().strip()
    msg = (
        "Measure retired instructions for fixed selected mesh sources\n\n"
        f"Signed-off-by: {cloud.PERSON}\n"
    )
    commit = cloud.git("commit-tree", tree, "-p", scoped.SOURCE, data=msg.encode()).decode().strip()
    cloud.git("update-ref", "refs/heads/" + branch, commit, "0" * 40)
    write(
        folder / "prepared.json",
        {
            "source": scoped.SOURCE,
            "workflow_source": commit,
            "branch": branch,
            "workflow": workflow,
            "repository": cloud.REMOTE,
        },
    )
    cloud.push(folder.name)


def audit(public_layers=False):
    import scoped
    from resource_report import paired

    folder = (
        scoped.OUT
        / "cloud"
        / ("linux-public-layers" if public_layers else "linux-worker-instructions")
    )
    evidence = folder / "evidence"
    plan = json.loads((folder / "plan.json").read_bytes())
    stages = PUBLIC_STAGES if public_layers else STAGES
    assert tuple(plan["sources"]) == stages and plan["orders"] == orders(stages)
    completed = json.loads((evidence / "completion.json").read_bytes())
    assert (
        digest(folder / "plan.json") == digest(evidence / "plan.json") == completed["plan_sha256"]
    )
    assert digest(evidence / "samples.jsonl") == completed["samples_sha256"]
    assert json.loads((evidence / "capability.json").read_bytes())["worker_instructions_confirmed"]
    build = json.loads((evidence / "build-identities.json").read_bytes())
    measured = json.loads((evidence / "measured-binaries.json").read_bytes())
    assert measured["worker_probe_sha256"] == digest(evidence / "worker-probe.rs")
    assert (
        measured["object_sha256"]
        == json.loads((evidence / "capability.json").read_bytes())["counter_sha256"]
    )
    assert (
        build["plan_sha256"] == completed["plan_sha256"]
        and build["controller_sha256"] == plan["build_driver_sha256"]
    )
    for stage, identity in plan["sources"].items():
        assert all(build["binaries"][stage][k] == v for k, v in identity.items())
        assert (
            measured["stages"][stage]["library_sha256"]
            == build["binaries"][stage]["library_sha256"]
        )
    samples = [json.loads(line) for line in (evidence / "samples.jsonl").read_text().splitlines()]
    assert [(r["round"], r["position"], r["case"], r["stage"]) for r in samples] == [
        (i, j, c, s) for i, j, c, _, s in schedule(plan)
    ]
    assert len(samples) == plan["processes"] == completed["processes"] == 480
    summary = {}
    for case, _ in plan["cases"]:
        rows = {(r["round"], r["stage"]): r for r in samples if r["case"] == case}
        for r in rows.values():
            assert (
                r["calls"] == 1
                and r["instructions"] > 0
                and 0.999 <= r["running"] / r["enabled"] <= 1
            )
        shapes = {
            s: {
                tuple(rows[i, s][k] for k in ("faces", "points", "triangles", "refusals"))
                for i in range(8)
            }
            for s in stages
        }
        assert all(len(v) == 1 for v in shapes.values())
        medians = {
            s: statistics.median(rows[i, s]["instructions"] for i in range(8)) for s in stages
        }
        pairs = [(b, a) for a, b in pairwise(stages)] + [(s, stages[0]) for s in stages[2:]]
        compare = {
            a + "/" + b: paired(
                [(rows[i, a]["instructions"], rows[i, b]["instructions"]) for i in range(8)]
            )
            for a, b in pairs
        }
        summary[case] = {
            "medians": medians,
            "ratios": compare,
            "shapes": {s: sorted(v) for s, v in shapes.items()},
        }
    write(folder / "summary.json", summary)
    write(
        folder / "audit.json",
        {
            "passed": True,
            "processes": len(samples),
            "mesh_calls": len(samples),
            "minimum_running_fraction": min(r["running"] / r["enabled"] for r in samples),
            "samples_sha256": completed["samples_sha256"],
            "summary_sha256": digest(folder / "summary.json"),
        },
    )


if __name__ == "__main__":
    if sys.argv[1] == "measure":
        measure(*(Path(s).resolve() for s in sys.argv[2:]))
    else:
        {"prepare": prepare, "audit": audit}[sys.argv[1]](sys.argv[2:] == ["public-layers"])
