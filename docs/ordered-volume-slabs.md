# Ordered parallel volume slabs

Large through-cut air proofs can spend nearly all their time integrating ray
intervals across slices. Classifier construction and arrangement break creation
can be inexpensive even when the complete integral takes minutes.

The height breaks already partition the outer integral into independent slabs.
Each slab retains the same Gauss points, adaptive refinement, tolerance share,
inner break generation and ray queries. Its unanswered-line flag is local:
once a required line is unanswered the slab has no value. Later local samples
still skip geometry, as they did in the serial pass.

At most eight slabs are evaluated in a batch, subject to `HAECCEITY_THREADS` and
available parallelism. The caller consumes results in ascending slab order and
performs the original additions in that order. A required unanswered result
refuses the pass. A speculative unanswered result beyond an already sufficient
prefix does not invalidate that prefix, because the serial pass would not have
visited it. Acceptance still finishes every pass required by `fills_at_most`;
speculation never turns an unanswered volume into air.

Small arrangements with fewer than 64 participating faces remain serial. Calls
already bound to a face or volume worker remain serial too. Worker binding is
independent of whether a timing observer exists. A failed thread spawn falls
back to calculation on the caller, without changing result order. A batch is
joined before its results are consumed, so early termination can leave at most
that bounded batch of extra work.

Directed checks compare the actual B-Rep volumes bit for bit and at-most answers
between serial and parallel execution, including curved slivers, empty and full
probes, and unanswerable source faces. Separate tests retain ordered rounding,
required refusals, speculative refusals beyond a sufficient prefix, and the
absence of nested worker pools. Complete native-document comparisons, standard
ARM snapshot/workspace regressions and bounded whole-part GUI acceptance remain
separate validation gates. No geometry tolerance or output format changes.
