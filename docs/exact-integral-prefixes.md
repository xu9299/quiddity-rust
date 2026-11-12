# Exact quadrature prefixes

The inner surface integral often revisits a long sequence of complete knot panels
while only its final endpoint changes. Exact node reuse avoids evaluating a sample
again, but still repeats its weight multiplication and ordered addition.

`IntegralPrefixes` stores the running accumulator after each complete panel. A
query may resume only after the same ordered endpoint bit patterns, including the
same starting endpoint. The other surface parameter must also have the identical
finite bit pattern. Changing it clears both sample and prefix caches. Signed zero
and last-bit differences therefore cannot share an accumulator accidentally.

Each stored accumulator is exactly the state reached by the original node loop.
The remaining panels execute that same loop. This uses neither a sum of separate
panel totals nor subtraction of cumulative sums: both would regroup floating-point
operations. A differing endpoint truncates the reusable prefix. A shorter query
can use a matching earlier state, but a reversed walk has a different prefix.

Only finite accumulators and endpoints are retained. Storage is bounded to 1024
panels per active fixed parameter and lives within one face-integral call. The
first visit to a parameter still follows the existing uncached path. Public data,
quadrature nodes, weights, tolerances, refusal rules and output order are unchanged.

Directed tests compare floating-point bits across cancellation, longer and shorter
queries, reversed traversal, changed endpoints, signed zero, nonfinite results and
fixed-parameter changes. They also count actual skipped node callbacks and enforce
the storage bound. The separate manual face probe compares prefix reuse off/on/on/off
with node reuse enabled throughout, using only explicitly selected faces of an
unchanged caller-provided STEP. It is not a substitute for complete API and strict
snapshot validation.
