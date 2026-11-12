# Windows integral comparison

I measured direct mass and moment calculations on three public upstream fixtures:
`held_side_knot.step`, `late_turn.step`, and `rotated_seam_band.step`. The raw JSON
retains the returned term vectors as integer bit patterns and per-fixture CPU/wall
timings. The probe checks every repeated result against its expected vector. Each
of the four variants performed five rounds and 3,840 evaluations; all six successful
term vectors were identical across variants. These are measured results, not new
snapshot expectations. No customer geometry or local native binary is included.

Use Windows x64, Rust 1.99.1 and a configured MSVC developer shell. The four
`*-identity.json` files identify the measured commits and compiler. At each exact
commit, insert `integral_benchmark.rs` immediately before the final closing brace
of the existing `tests` module in `crates/haecceity/src/mass.rs`. Keep the test local.
Set `CARGO_BUILD_JOBS=2`, `RUST_TEST_THREADS=1`, `RUST_MIN_STACK=67108864`, and
`QUIDDITY_SPLIT_BENCH_OUTPUT` to an absolute output JSON path. Then run:

```text
cargo test -p haecceity --lib --release --locked compare_upstream_integral_cost -- --ignored --nocapture
```

Restore the test insertion before checking out the next commit. Run main-before,
operands, samples, then main-after. Compare the complete `expected` objects for
exact equality using an integer-preserving JSON parser: their entries are u64
floating-point bit patterns. Sum each round's thread CPU times separately from
its wall times.
The probe warms geometry caches but bypasses mass/moment result caches. STEP reads,
compilation and full recognition are excluded. The current-thread CPU counter is
Windows-specific, so this is not an ARM timing result.

The two control means differ by 9.54%. Operand reuse reduced mean thread CPU time
by 14.67% on these fixtures; sample reuse increased it by 2.19%, and its timing range
overlaps the controls. I am holding sample reuse for better performance evidence.
The new review base `4a284e08a2faee1041e14726d0d621ffc2427bb2` and its operand candidate have the same haecceity
source trees as the original measured versions. The measurements keep their old
identities; validation of the newly based commits is recorded separately.
`review-source-migration.json` records the source-tree and file-blob comparisons.

The final review base is the 0.5.0 release preparation `18eefdc10c07dca3d1d0a3de19fe59e93713b9c6`. All Rust files are unchanged from the preceding #264 base. Its two updated snapshots differ only at 15 package/version fields. See `release-transition.json`. The final candidate heads receive fresh Windows directed checks and ARM full checks; the 56 earlier local Windows test executions keep their original source identities.
