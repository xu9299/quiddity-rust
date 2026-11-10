# Shared-analysis instruction benchmark on Windows

This optional development harness compares separate recognition/evidence calls with
`recognise_with_evidence` using software dynamic instruction instrumentation. It uses
[DynamoRIO](https://dynamorio.org/) and an x64 MSVC developer environment; neither the
client nor DynamoRIO is linked into or required by quiddity.

## What is counted

The client counts executed x86-64 application instructions whose original addresses
are inside quiddity's `_native.pyd`, summed across threads between explicit call
markers. This includes native parsing, analysis and JSON serialization. Instructions
in imported libraries, the Python interpreter, the kernel, and the instrumentation
itself are outside this scope. This is a software dynamic instruction count, not a
hardware PMU reading. No CPU-cycle or elapsed-time conversion is involved.

Each basic block contributes its application-instruction count when it executes.
The client requires DynamoRIO's `-thread_private` mode and embeds a separate counter
for each thread. Marker snapshots subtract all earlier work, including imports and
warm-up calls. The chosen DLL is matched by its filename; keep it loaded throughout
the measurement and use a process dedicated to the benchmark.

Counter updates and boundary reads use MSVC's atomic 64-bit intrinsics, so sampling
a worker's final native instructions does not race a plain C counter access.

The calibration DLL also exports a no-op boundary function. The client intercepts
that function, validates its operation and label, and excludes the marker itself
from calibration counts. Boundary reads fail explicitly instead of silently dropping
a marker. No Windows debug-output or hardware-profiling API is used.

This block-level software measurement does not classify hardware-retired events or
count individual iterations of a `rep` instruction separately. The assembly calibration
checks the counting mechanism; complete-output comparisons check that instrumentation
has not changed the result. Allocation and worker scheduling can still cause small
differences between repeated counts, which are retained in the report.

## Build and run

Use the uninstrumented release Python development binding from the same checkout.
For a Windows Cargo build, copy `target/release/_native.dll` to
`crates/quiddity-py/python/quiddity/_native.pyd` as an independent file after the build.
Use a copy rather than a hard link so subsequent Cargo builds cannot replace the
measured binding. The harness records and rechecks its SHA-256.

Download an official Windows x64 DynamoRIO release and verify its published archive
digest. This harness was developed with DynamoRIO 11.91.20763. In an x64 MSVC developer
PowerShell, with the Windows SDK headers available:

```powershell
cargo build -p quiddity-py --release --locked
Copy-Item target/release/_native.dll crates/quiddity-py/python/quiddity/_native.pyd
./tools/shared_analysis/build_counter.ps1 -DynamoRIO C:/tools/DynamoRIO -OutputDirectory outputs/icount-client
uv run --no-project --python 3.14 tools/shared_analysis/benchmark_instructions.py --dynamorio C:/tools/DynamoRIO --counter outputs/icount-client/instruction_counter.dll --corpus $env:QUIDDITY_CORPUS --output outputs/shared-instructions
```

Use a new output directory for each run. All Python dependencies of this harness are
in the standard library. The controller uses the base interpreter selected by `uv`
when launching under DynamoRIO, avoiding Windows virtual-environment launcher child
processes. It sets the source checkout's Python package path explicitly. The client
build produces only two development DLLs: the counter and its calibration fixture.
No administrator privileges or system-wide profiler registration are needed. Profiler
configuration and temporary files stay under the selected output directory. Git metadata
is collected before instrumentation, so the measured process does not launch Git.

## Counter validation and measurement

Before measuring quiddity, the harness verifies:

- An empty region containing Python work counts zero target-module instructions.
- The assembly loop executes exactly `3*n + 4` instructions for `n = 0, 1, 10, 1000`.
  An earlier warm-up invocation is excluded.
- Calls from two threads contribute exactly 62 instructions in total, with two
  contributing counters. This checks the private-cache and thread-summing behavior.
- Another 128 alternating begin/end regions have exact known loop counts, exercising
  more boundary transitions than the complete five-case benchmark.

The five inputs are the same public cases used in the elapsed-time comparison:
the through-hole, bolt-circle/grid, blind-pocket-pattern and turned-disc fixtures,
plus corpus `inventory_refusal/13975.step.gz`. Both paths warm once before each case;
eight paired measurements alternate order, four pairs in each order. Frame selection
is read before measurement. Complete recognition and evidence documents are compared
after every pair, outside the measured region, with no floating-point tolerance.
The instrumented outputs must also equal an uninstrumented reference run.

`results.json` contains the input, tool-source and binary hashes, calibration results, all per-pair
counts, contributing thread counts, medians and reductions. Raw counter records,
logs and complete outputs remain alongside it. The script fails on an incomplete
counter trace, a calibration mismatch, changed input/binding bytes or unequal output.
These instruction counts describe the stated native-module scope and do not predict
end-to-end GUI latency or every input's performance.
