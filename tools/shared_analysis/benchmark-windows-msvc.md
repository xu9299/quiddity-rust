# Shared-analysis native instruction benchmark

This benchmark compares `recognise` plus `evidence` with `recognise_with_evidence`
on the same release build. Its metric is a software count of executed x86-64
application instructions inside quiddity's `_native.pyd`, summed across threads.
It is not a hardware-PMU retired-instruction measurement. Imported libraries,
the Python interpreter, kernel execution and instrumentation are outside the count.

The [optional harness](README.md) uses
[DynamoRIO 11.91.20763](https://github.com/DynamoRIO/dynamorio/releases/tag/cronbuild-11.91.20763)
without administrator privileges. The downloaded Windows archive has SHA-256
`ea3dad4ae6e39228fe87b5d4d7bc3a53fafdb6161ee5ba6a2e98fb64f8a7e0d4`.
The counter follows DynamoRIO's [basic-block instrumentation API](https://dynamorio.org/API_samples.html)
and [function-wrapping API](https://dynamorio.org/group__drwrap.html).

## Build and sampling

- Engine source: `54d2b7a3e4e97e91f4705453ae6f0fef4f6ecc00`, based on
  `92ad8c6b3a630c77293d8decb9a976741f13063c`.
- Windows 10 x64, Intel Xeon E5-1630 v3, Rust 1.99.1/MSVC, Python 3.14.6.
- Native binding SHA-256:
  `39994ca03abc347857c79689b2e3f5a793abd4fd59e32f9ddaf1802a724b7c6e`.
  The binding is an independent copy of the release DLL, protected from later Cargo relinks.
- Corpus: `pzfreo/quiddity-corpus@3ef06590a7a27ed65ca1f1a6e7f7ea4d5c328f1f`.
- Each path warms once per input. Eight paired measurements alternate call order,
  with four pairs in each order. Frame selection is determined before the regions.
- Regions contain only the native API calls. JSON decoding, comparison and output-file
  writes occur outside them. Complete recognition and evidence documents must match
  after every pair and also match an uninstrumented reference run, without tolerance.

The client uses thread-private code caches and per-thread counters. A no-op function
in the calibration DLL marks each region and returns normally. Its own instructions
are excluded. Before profiling, 134 regions check empty-region exclusion, exact
assembly-loop counts, two-thread aggregation and repeated boundary transitions.
Any mismatched boundary, incomplete trace, changed input/binding digest or unequal
output causes the run to fail.

## Interpretation

The metric measures the native work removed by sharing analysis, while retaining
both complete outputs. Small per-trial count changes can arise from allocation or
worker scheduling; the report retains every trial. These counts do not predict GUI
latency or establish a fixed saving for other inputs. Instrumented wall-clock times
are not used as performance results.

## Results

All 40 pairs returned equal complete recognition and evidence documents, with
`evidence_error: null`, and matched the uninstrumented reference. The table uses the median of eight
measurements per mode. Counts are in millions of native-module instructions.

| Input | Separate calls | Joint call | Reduction |
| --- | ---: | ---: | ---: |
| `golden_simple_through_hole.step` | 824.343 | 471.306 | 42.8% |
| `golden_bolt_circle_and_rectangular_grid.step` | 2930.210 | 1810.925 | 38.2% |
| `golden_blind_pockets_and_pocket_patterns.step` | 1052.142 | 633.557 | 39.8% |
| `turned_body_disc.step` | 864.921 | 503.680 | 41.8% |
| `@corpus/inventory_refusal/13975.step.gz` | 937.900 | 526.663 | 43.8% |

40/40 paired trials reduced the instruction count. [Measurement metadata](results/windows-msvc.json), [all 80 raw call regions](results/windows-msvc-counts.jsonl) and [all 134 calibration regions](results/windows-msvc-calibration.jsonl) record the source/input/binary identities and every count. The harness also saves complete outputs locally for inspection.
