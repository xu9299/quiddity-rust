"""核对一份 Linux 二进制在原生环境及 Callgrind 下生成的逐面网格。"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from ablation import BASE, HEAD

HERE = Path(__file__).resolve().parent


def sha(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def main(source, corpus, evidence):
    evidence.mkdir(exist_ok=True)
    step = corpus / "cadgenbench_inputs/cgb243.step.gz"
    report = {"corpus": subprocess.check_output(["git", "-C", corpus, "rev-parse", "HEAD"]).decode().strip(),
              "step_sha256": sha(step), "probe_sha256": sha(HERE / "probe.rs"),
              "valgrind": subprocess.check_output(["valgrind", "--version"]).decode().strip(),
              "rustc": subprocess.check_output(["rustc", "-Vv"]).decode(),
              "lscpu": subprocess.check_output(["lscpu"]).decode(), "stages": {}}
    for stage, commit in (("pr349", BASE), ("pr364", HEAD)):
        subprocess.run(["git", "checkout", "--detach", commit], cwd=source, check=True)
        with (evidence / f"{stage}-build.log").open("wb") as log:
            raw = subprocess.check_output(["cargo", "build", "--locked", "--release", "-p", "haecceity",
                                          "--message-format=json"], cwd=source, stderr=log)
        events = [json.loads(row) for row in raw.splitlines()]
        def library(name):
            return next(Path(f) for e in events if e.get("target", {}).get("name") == name
                        for f in e["filenames"] if f.endswith(".rlib"))
        lib, serde = library("haecceity"), library("serde_json")
        binary = Path(os.environ["RUNNER_TEMP"]) / f"face-counts-{stage}"
        subprocess.run(["rustc", "--edition=2024", "-C", "opt-level=3", "-L", f"dependency={serde.parent}",
                        "--extern", f"haecceity={lib}", "--extern", f"serde_json={serde}", HERE / "probe.rs", "-o", binary], check=True)
        stats = {}
        for mode in ("native", "callgrind", "native-repeat", "native-glibc-no-fma"):
            output = evidence / f"{stage}-{mode}.json"
            cmd = [binary, "counts", step, "default", output]
            env = os.environ.copy()
            if mode == "callgrind":
                cmd = ["valgrind", "--tool=callgrind", "--cache-sim=no", "--branch-sim=no",
                       "--callgrind-out-file=" + str(binary) + ".cg", *cmd]
            if mode == "native-glibc-no-fma":
                env["GLIBC_TUNABLES"] = "glibc.cpu.hwcaps=-FMA,-FMA4"
            subprocess.run(cmd, env=env, check=True)
            faces = json.loads(output.read_bytes())
            stats[mode] = {"sha256": sha(output), "faces": len(faces),
                           "triangles": sum(f.get("triangles", 0) for f in faces),
                           "refusals": sum("refusal" in f for f in faces)}
            print(stage, mode, stats[mode], flush=True)
        assert stats["native"] == stats["native-repeat"]
        report["stages"][stage] = {"commit": commit, "library_sha256": sha(lib), "probe_sha256": sha(binary),
                                   "modes": stats}
        (evidence / "linux-geometry.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf8")


if __name__ == "__main__":
    main(*(Path(s).resolve() for s in sys.argv[1:]))
