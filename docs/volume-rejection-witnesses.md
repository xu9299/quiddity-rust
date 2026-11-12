# Early rejection from original volume slabs

The at-most material query with a limit in `[0, 1)` may reject from an individual
original integration slab. It first checks at most eight spread-out slab indices,
using the existing bounded workers. These are complete B-Rep integration slabs
with the original world axes, arrangement breaks, quadrature, tolerances and ray
queries. There are no new geometric approximations or added sampling points.

Every integrated length and quadrature weight is nonnegative. Consequently any
slab with material rules out an empty inset integral, and any raw slab already
beyond the fraction limit also puts the whole raw integral beyond that limit.
A required unanswered slab prevents acceptance. These observations allow early
rejection, including when an obstruction occurs late in the original order.

A witness value is marked partial. It cannot establish an at-most answer. When
the raw integral might be small enough, the full required inset pass is completed
with its original empty/full conventions and refusal checks. Completed values
are cached and consumed in the original order, retaining the original floating
point sums. Other volume APIs and their ordered-prefix behavior are unchanged.

Directed tests cover late material, unanswered required slabs, exact ordered
rounding, curved partial volumes, full-volume conventions and generic rigid
motions of many-sided through-cuts, blind cuts and thin obstructions. Native API
byte comparisons, strict snapshots and whole-part GUI acceptance remain separate
validation gates. No result format changes are involved.
