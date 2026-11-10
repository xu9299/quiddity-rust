"""Compare native-module instruction counts with an optional DynamoRIO client."""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import threading
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[2]
CASES = (
    "golden_simple_through_hole.step",
    "golden_bolt_circle_and_rectangular_grid.step",
    "golden_blind_pockets_and_pocket_patterns.step",
    "turned_body_disc.step",
    "@corpus/inventory_refusal/13975.step.gz",
)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def make_marker(counter):
    emit = ctypes.WinDLL(str(counter.parent / "icount_calibration.dll")).icount_marker
    emit.argtypes = [ctypes.c_int, ctypes.c_char_p]
    emit.restype = None

    def mark(kind, label):
        emit(int(kind == "BEGIN"), label.encode("ascii"))

    return mark


def calibrate(args):
    mark = make_marker(args.counter)
    library = ctypes.WinDLL(str(args.counter.parent / "icount_calibration.dll"))
    loop = library.icount_loop
    loop.argtypes = [ctypes.c_uint32]
    loop.restype = ctypes.c_uint32
    assert loop(1000) == 1000  # Outside every measured region.
    expected = {}
    mark("BEGIN", "empty")
    assert sum(range(1000)) == 499500
    mark("END", "empty")
    expected["empty"] = 0
    for n in (0, 1, 10, 1000):
        label = f"calibrate_{n}"
        mark("BEGIN", label)
        assert loop(n) == n
        mark("END", label)
        expected[label] = 3 * n + 4
    mark("BEGIN", "two_threads")
    assert loop(7) == 7
    thread = threading.Thread(target=loop, args=(11,))
    thread.start()
    thread.join()
    mark("END", "two_threads")
    expected["two_threads"] = (3 * 7 + 4) + (3 * 11 + 4)
    # Exercise more boundary transitions than the complete five-case benchmark.
    for index in range(128):
        label = f"boundary_{index}"
        n = index % 17
        mark("BEGIN", label)
        assert loop(n) == n
        mark("END", label)
        expected[label] = 3 * n + 4
    write_json(args.output / "calibration-expected.json", expected)


def worker(args):
    import quiddity
    from quiddity import _native

    mark = make_marker(args.counter)
    native = Path(_native.__file__).resolve()
    assert native.is_relative_to(SOURCE)
    native_sha = digest(native)
    report = {
        "source_commit": args.source_commit,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "native_sha256": native_sha,
        "repetitions": args.repeats,
        "scope": (
            "Executed application instructions inside quiddity's _native.pyd, "
            "summed across threads within explicit API-call regions. Imported "
            "libraries, kernel code and Python interpreter instructions are excluded."
        ),
        "cases": [],
    }
    for case_index, name in enumerate(CASES):
        path = (
            (args.corpus / name.removeprefix("@corpus/"))
            if name.startswith("@corpus/")
            else SOURCE / "tests/fixtures" / name
        )
        source_sha = digest(path)
        recognition = json.loads(quiddity.recognise(str(path)))
        framed = recognition["document"]["frame"]["status"] == "inferred"

        def separate(path=path, framed=framed):
            return quiddity.recognise(str(path)), quiddity.evidence(str(path), framed=framed)

        def shared(path=path):
            return quiddity.recognise_with_evidence(str(path))

        def compare(left, right, name=name, source_sha=source_sha):
            result = json.loads(right)
            assert result["recognition"] == json.loads(left[0]), name
            assert result["evidence"] == json.loads(left[1]), name
            assert result["evidence_error"] is None, name
            assert result["recognition"]["document"]["source_sha256"] == source_sha
            return result

        expected = compare(separate(), shared())  # Both paths warm outside regions.
        pairs = []
        for repetition in range(args.repeats):
            order = ("separate", "shared") if repetition % 2 == 0 else ("shared", "separate")
            values = {}
            regions = {}
            for mode in order:
                label = f"case{case_index}_{mode}_{repetition}"
                mark("BEGIN", label)
                values[mode] = separate() if mode == "separate" else shared()
                mark("END", label)
                regions[mode] = label
            assert compare(values["separate"], values["shared"]) == expected
            pairs.append({"order": list(order), "regions": regions})
        assert digest(path) == source_sha
        report["cases"].append(
            {
                "input": name,
                "sha256": source_sha,
                "framed": framed,
                "complete_outputs_equal": True,
                "pairs": pairs,
            }
        )
        write_json(args.output / f"case{case_index}-output.json", expected)
        write_json(args.output / "worker.json", report)
        print(f"Completed {name}", flush=True)
    assert digest(native) == native_sha
    report["complete"] = True
    write_json(args.output / "worker.json", report)


def run_profile(args, mode, module):
    raw = args.output / f"{mode}-counts.jsonl"
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(SOURCE / "crates/quiddity-py/python")
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PATH"] = (
        str(args.dynamorio / "ext/lib64/release") + os.pathsep + environment["PATH"]
    )
    config = args.output / "dynamorio-config"
    temporary = args.output / "tmp"
    config.mkdir(exist_ok=True)
    temporary.mkdir(exist_ok=True)
    environment["DYNAMORIO_CONFIGDIR"] = str(config)
    environment["TEMP"] = environment["TMP"] = str(temporary)
    command = [
        str(args.dynamorio / "bin64/drrun.exe"),
        "-quiet",
        "-thread_private",
        "-c",
        str(args.counter),
        str(raw),
        module,
        "--",
        sys._base_executable,
        str(Path(__file__).resolve()),
        "--worker",
        mode,
        "--counter",
        str(args.counter),
        "--corpus",
        str(args.corpus),
        "--output",
        str(args.output),
        "--repeats",
        str(args.repeats),
        "--source-commit",
        args.source_commit,
    ]
    with (args.output / f"{mode}.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
    rows = [json.loads(line) for line in raw.read_text().splitlines()]
    assert rows[-1]["complete"] and rows[-1]["regions"] == len(rows) - 1
    counts = {row["region"]: row for row in rows[:-1]}
    assert len(counts) == len(rows) - 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dynamorio", type=Path)
    parser.add_argument("--counter", required=True, type=Path)
    parser.add_argument("--corpus", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=8)
    parser.add_argument("--worker", choices=("calibrate", "measure"))
    parser.add_argument("--calibrate-only", action="store_true")
    parser.add_argument("--source-commit", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("dynamorio", "counter", "corpus", "output"):
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.resolve())
    if args.worker:
        calibrate(args) if args.worker == "calibrate" else worker(args)
        return
    if args.dynamorio is None:
        parser.error("dynamorio is required unless running a worker")
    if args.repeats < 4 or args.repeats % 2:
        parser.error("repeats must be even and at least four")
    args.source_commit = subprocess.check_output(
        ["git", "-C", str(SOURCE), "rev-parse", "HEAD"], text=True
    ).strip()
    tool_sources = (
        "instruction_counter.c",
        "calibration.asm",
        "build_counter.ps1",
        Path(__file__).name,
    )
    source_hashes = {name: digest(Path(__file__).parent / name) for name in tool_sources}
    args.output.mkdir(parents=True, exist_ok=False)
    calibration = run_profile(args, "calibrate", "icount_calibration.dll")
    expected = json.loads((args.output / "calibration-expected.json").read_text())
    assert {key: row["instructions"] for key, row in calibration.items()} == expected
    assert calibration["two_threads"]["threads"] == 2
    write_json(args.output / "calibration.json", calibration)
    print("Instruction counter calibration passed", flush=True)
    if args.calibrate_only:
        return
    reference = args.output / "reference"
    reference.mkdir()
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(SOURCE / "crates/quiddity-py/python")
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    with (args.output / "reference.log").open("w", encoding="utf-8") as log:
        subprocess.run(
            [
                sys._base_executable,
                str(Path(__file__).resolve()),
                "--worker",
                "measure",
                "--counter",
                str(args.counter),
                "--corpus",
                str(args.corpus),
                "--output",
                str(reference),
                "--repeats",
                "0",
                "--source-commit",
                args.source_commit,
            ],
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    counts = run_profile(args, "measure", "_native.pyd")
    report = json.loads((args.output / "worker.json").read_text())
    assert report["complete"] and len(counts) == len(CASES) * args.repeats * 2
    reference_report = json.loads((reference / "worker.json").read_text())
    assert reference_report["complete"]
    assert reference_report["native_sha256"] == report["native_sha256"]
    for index in range(len(CASES)):
        name = f"case{index}-output.json"
        assert json.loads((reference / name).read_text()) == json.loads(
            (args.output / name).read_text()
        )
    report["uninstrumented_complete_outputs_equal"] = True
    report["counter"] = {
        "tool": "DynamoRIO software dynamic instruction instrumentation",
        "client_sha256": digest(args.counter),
        "calibration_library_sha256": digest(args.counter.parent / "icount_calibration.dll"),
        "runtime_sha256": digest(args.dynamorio / "lib64/release/dynamorio.dll"),
        "tool_sources_sha256": source_hashes,
        "calibration": calibration,
        "hardware_pmu": False,
    }
    assert source_hashes == {name: digest(Path(__file__).parent / name) for name in tool_sources}
    for case in report["cases"]:
        for pair in case["pairs"]:
            pair["instructions"] = {
                mode: counts[region]["instructions"] for mode, region in pair["regions"].items()
            }
            pair["threads"] = {
                mode: counts[region]["threads"] for mode, region in pair["regions"].items()
            }
            assert all(n > 0 for n in pair["instructions"].values())
        medians = {
            mode: statistics.median(pair["instructions"][mode] for pair in case["pairs"])
            for mode in ("separate", "shared")
        }
        case["median_instructions"] = medians
        case["reduction_percent"] = 100 * (1 - medians["shared"] / medians["separate"])
        print(
            json.dumps(
                {"input": case["input"], **medians, "reduction_percent": case["reduction_percent"]}
            ),
            flush=True,
        )
    write_json(args.output / "results.json", report)


if __name__ == "__main__":
    main()
