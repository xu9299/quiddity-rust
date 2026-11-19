# Face integration and boundary selection

This change retains the integration densities, quadrature nodes, operation order,
initial cuts, error criteria, depth limit and kink budget. It does not infer failure
from elapsed time, face size, spline degree or a large parent-panel error.

## Reusing work

`face_moments` stores its complete result in the immutable part's face cache, including
`None`; all moment components are retained. A transformed part starts with fresh
geometry caches, and sense repairs also replace the face caches. The mass
cache introduced by the existing shared moment/mass integration retains its own rule.

Within one face integral, the NURBS evaluator borrows one immutable surface and
retains each axis's last basis and derivative by the exact bits of its clamped
parameter. A changed span pair rebuilds the homogeneous control block; non-finite
parameters use the uncached path. Reuse preserves u/v/component accumulation order,
quadrature nodes, convergence criteria and refusal conditions. These operands are
discarded with the integral.

An adaptive panel requests kink locations only if there is a refinement floor, its
depth is below `MAX_DEPTH`, and its interval still has kink budget. Otherwise the
location cannot affect the next step. All quadrature nodes and the walking state are
still evaluated. The discarded corner search has local state, cannot change the walk,
and does not propagate its projection failures. The initial boundary cuts are untouched.

## A certificate for an inevitable refusal

This is a proof about the existing finite algorithm, not a proof that the geometric
area or moments do not exist, or that a different integrator could not compute them.
Gauss--Kronrod differences are numerical error estimates, not rigorous bounds on the
true integration error.

For one requested density, the coarse pass fixes
`F = 1e-11 * sum(coarse term magnitudes)`. It is not recomputed during refinement.
For an accepted leaf with estimated error `e` and magnitude `s`, the original
adaptive walk updates its unresolved excess as follows:

```
if !(e <= max(1e-9 * s, 0)):
    E = max(E, e)
```

`E` starts at zero. Boundary aggregation also takes a maximum, and the caller rejects
the refined loop when `E > F`. The face does not proceed to loop orientation or publish
an integral in that case. The mass and moment refinements each use their own fixed `F`.

Suppose an accepted leaf makes `E > F`. Every subsequent maximum is at least this `E`.
Thus either a later evaluation fails (the original returns `None`), or refinement
finishes with excess still greater than `F` (the original also returns `None`).
No continuation can produce a successful face result. The optimized walk may return
that rejected prefix; the boundary propagates it immediately, and the unchanged
caller applies the same rejection. In particular, projected mass refinement
keeps a known refusal when the master trim was unavailable. No rejected prefix
is published.

This does not reject a parent that can still be subdivided. With ordinary finite
values, an accepted leaf that fails both the relative criterion and `F` was accepted
only because it reached `MAX_DEPTH`. A large error accepted by the relative criterion
does not enter `E`. Equality with `F` is not refusal. With an infinite or NaN floor the
ordered `E > F` guard cannot trigger. Once it does trigger, Rust's `f64::max` cannot
reduce the non-NaN `E`, even if a later estimate is NaN. Evaluation failures remain
failures; no new numerical answer is invented.

## STEP face boundaries

A STEP surface curve has a spatial curve and up to two associated geometric
representations. Its master flag selects the 3D curve or one of those slots.
The pcurve's two coordinates belong to its named carrier. For a different
incident face, the evaluator lifts it on its carrier and locates that spatial
point on the receiving face. A pcurve outside its carrier's finite knot domain
has no valid lift. Clamping it to the edge of that patch changes both the
geometry and the tangent. A retained analytic carrier can also stand alone
without a topological face. A spatial master goes through the existing
projection. Both choices may appear in the same loop. An explicit
`FACE_OUTER_BOUND` supplies the outer contour. Where the file gives only
ordinary face bounds, loop placement infers which one encloses the others.

Analytic elementary surfaces and native B-spline and Bézier patches retain
their own parameter units. For a converted swept patch the mapping accounts
for the extrusion vector's magnitude or the rational quadratic parametrisation
of each quarter revolution. Edge vertices delimit each pcurve in its own
parameter. Periodic seams can shift a parameter by a whole turn, and a
collapsed side can have several parameter values for one physical point.

The kernel first certifies the declared region. The specialised integrator for
one open spline patch remains available. Other evaluable regions use the
general boundary integral, including holes and loops that wind around a
periodic direction. Its density and orientation are those of the receiving
surface. The same accepted region supplies area, flux, moments, containment,
trim bounds and edge queries. Solid volume consumes the face flux. A whole
sphere with a vertex loop uses its analytic formula.

Endpoint matches and closing paths are checked against the file's length
uncertainty plus the coordinate floor. Interior checks sample the union of the
two representations' knot spans in physical coordinates. A sampled departure
over one thousandth of the boundary's reference diagonal, or a departure times
perimeter greater than one thousandth of area, rejects the declared region.
These are sampled representation checks. Independent curve parameters can
advance at different rates on the same spatial arc. The spatial curve remains
available if a declaration cannot be evaluated or the integral fails its error
checks. The mass cache publishes the chosen method along with its result so
the first query for a centroid cannot choose another boundary.

The bowed edge in `master-trim-square.step` has height
`y = 2 h x (1 − x)` for `0 ≤ x ≤ 1`. Its 3D approximation is straight.
The two areas are 1 and `1 − h/3` mm². For `h = 0.0004` the file selects
0.9998666667 mm². Regression cases use the same construction on a plane,
a trimmed patch with a hole, a seam-crossing cylinder and swept surfaces.
STEP export retains selected raw pcurves and the explicit outer bound without
requiring an area calculation.

## Settled volume queries

A volume pass refuses as soon as a sampled line is unresolved. Later samples
cannot restore an answer; placeholder zeros only finish discarded integration
work and must never be published as material or air. Successful passes retain
their sampling and accumulation order.

For `fills_at_most`, a material prefix can stop the inset pass when `limit < 1`,
even before air is found; otherwise both material and air are needed. The exact
pass can reject once `material / whole > limit`: lengths are nonnegative and
quadrature weights positive, so its running total cannot decrease. Equality does
not reject. Acceptance requires finishing any shortened inset pass, preserving
the empty/full conventions and refusing unanswered lines. An unknown or
nonpositive whole volume also refuses the at-most query.

## Verification obligations

The tests compare exhaustive and early walks, including accepted sums and end states,
ascending/descending intervals, a bad parent resolved by subdivision, acceptance by
the relative criterion, adjacent floating-point values around the floor, non-finite
floors and exhaustion of the depth budget. A rejected prefix is explicitly shown to
have an incomplete sum, so callers must check its excess before using that sum.
The existing quadrature and kink-budget tests retain exhaustive mode.

Independent face-area/moment references, placements, existing recognition snapshots,
inventory floors and tolerances remain unchanged. Repeated-moment speedups measure
cache reuse; they are not claims about cold integrations or whole-part recognition.
