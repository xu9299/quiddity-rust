# Computation reuse without additional compute workers

Official 0.5.3 supplies immutable ray geometry, NURBS operand reuse, section-recess
discovery reuse, one-time evidence serialization, settled volume-query exits and
wire-owner queries, as well as the released unit metadata and geometry fixes.
The project overlay retains integral sample and ordered-prefix reuse, document-only
recognition, thread filtering, early ineligible wall rejection, floor-proof ordering,
original-slab rejection witnesses and optional native timing. Released implementations
and their additional regression tests and source documentation are retained.

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
