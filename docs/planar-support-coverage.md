# Bounded straight planar support coverage

A straight finite edge has no sampling chord error. Widening only the scan-line
intervals can erase a physical strip in one parameter direction while leaving it
in the other. A 1 mm square with a 0.0001 mm missing strip reproduces this merely
by changing the first corner of its loop.

The additional path uses the same relative area budget, `1e-9`. It returns a
verdict only when a finite area enclosure is decisive. The existing sampled path
remains responsible for every other query; this change does not claim to resolve
all direction dependence in curved or uncertain geometry.

## Eligibility and finite representations

The target and every same-surface support must have closed, coherent loops. Each
edge must be a native line or a positive degree-one two-pole spline. Available file
pcurves on the queried face and explicitly preferred adjacent pcurves must also
have a bounded straight planar image. Positive two-pole pcurves keep their finite
domain. An unsupported representation cannot be omitted to establish a rejection.
The file's length unit is applied before projecting an affine pcurve construction.
A declared modelling precision is only an endpoint coherence check, never an extra
coverage tolerance.

The reference polygons retain literal topological vertices. Analytic finite feet
and available pcurve images are enclosed with outward arithmetic before rounding.
For a straight edge homotopy, the Jacobian determinant is affine; the mean absolute
determinant at its four corners bounds the complete swept image. Summing these
bounds encloses the change of the closed parity region, including its joins.
Only support error is clipped to the nominal target box. Target error retains
all geometry outside that box.

## Integration and decisions

Exact orientation signs identify crossing events. Outward intervals enclose their
heights, and all vertex heights remain cuts. Between cuts the active edge set is
fixed. A change of crossing order modifies a Boolean/parity scan length by at most
twice the pair separation; its maximum over the band occurs at an endpoint. These
complete pair strips retain uncertain events without paying the whole patch width.
Midpoint scan arithmetic is also enclosed, including overlapping crossing boxes.

Let `G` enclose nominal uncovered area, `A` enclose nominal target area, `E` bound
the complete target/support representation change, and `Et` bound target change.
Coverage requires `G.upper + E <= 1e-9 * (A.lower - Et)`, with a positive denominator.
Rejection requires `G.lower - E > 1e-9 * (A.upper + Et)`. All arithmetic in these
comparisons is rounded outward. Otherwise the new path returns no verdict.

The proof has work limits: at most 512 participating edges, bounded coordinate
magnitudes, and a cap on crossing-band work. These are eligibility limits, not
geometric epsilons. Curved and cylindrical geometry, ambiguous endpoints, invalid
numeric bounds and exhausted work budgets retain the existing query path.

## Validation scope

Synthetic analytic tests retain true strip gaps under loop restarts, rotations
and translations, distinguish a real hole from a filled union, and retain finite
pcurve-domain failures. Exact rational polygon clipping supplies an independent
area reference. Full workspace, corpus, snapshot, API and feature-off checks use
the repository's existing validation requirements; no reference is weakened.

A selected public narrow-channel document is also compared against its unchanged
complete Rust snapshot, retaining the caller-space integration path.

Representation kinds are preflighted before interval construction. Literal feet
are reused only for finite cardinal lines or an error-free endpoint sum; checking
just the rounded endpoint is insufficient. Zero and signed-unit interval identities
are exact, while invalid or unbounded arithmetic still cannot become a finite proof.

Nonfinite or reversed arithmetic inputs propagate an unbounded interval through
all operators, including multiplication by zero. They cannot provide finite proof
error by collapsing a NaN product or an invalid denominator.
