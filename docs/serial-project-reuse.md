# Computation reuse without additional compute workers

The project overlay retains immutable geometry reuse, integral operand and prefix
reuse, section-recess discovery reuse, document-only recognition, one-time evidence
serialization, early ineligible wall rejection and original-slab rejection witnesses.

Freeform recovery, face moments, planar floor readings and volume integration run
on the caller's thread. Recovery follows the original lazy query order. There is no
project face-query pool, heavy-face queue or volume-slab pool. Rejection witnesses
still check up to eight original slabs, sequentially, and cannot establish acceptance.
The optional timing observer reports progress only and performs no geometry work.

Independent parts may still be scheduled by the host application in separate
processes. The overlay does not change its batch concurrency, memory limits or
result format. The official release's unrelated HLR concurrency is unchanged.

Validation must distinguish correctness from throughput. Previous eight-worker
P0426121 timings do not establish this revision's timing. Required checks include
serial slab refusal/rounding, exact integral reuse, shared ray geometry, cached
ownership modes, complete native API bytes, frozen evidence and GUI save/reopen.
Whole-workspace and strict snapshot validation runs on GitHub Actions.
