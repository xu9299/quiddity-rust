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
geometry and the tangent. An elementary surface or a native spline carrier can stand alone
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
The perimeter is the physical length of the master on its STEP carrier. For a
foreign pcurve the integrand uses the carrier's partial derivatives directly.
Transferring the trim to the receiving surface is necessary for the area
integral. The length check uses the same adaptive error budget on its carrier.
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

A spatial query on a topological edge first resolves the selected carrier. A
plane lifts its parameter curves to lines, conics or NURBS in space, and
straight isoparameters of elementary quadrics can lift to lines or circles.
Section outlines, hidden lines, closest points, overlap sampling and feature
measurements use this geometry. A general surface curve is evaluated on its
carrier; a chord only locates the interval for a closest-point search. Mesh
coedges share the STEP vertex at their regular endpoints. The mass integral
still uses the actual trim and its checked closing segment.

Turned chamfers and countersinks use a cone's finite rims. A linear pcurve
with a small generator component on a cylinder or cone stays within half
its generator excursion of its middle parallel. For a planar rational
quadratic, the cone supplies a candidate circle at the edge's axial station.
On each quadratic Bézier span the degree-four Bernstein coefficients bound
the squared radial residual, and the control points bound the axial residual.
A resolved arc that stays within coordinate precision can supply a circular
rim for recognition. Its spline is still used for measurement.

On a plane, position is affine in both parameters and the normal is constant.
Consequently every area, flux and first- or second-moment density is at most
quadratic along the inner integration line. The two-point Gauss rule is exact
there. This avoids a spurious refusal on a short, fast-parameterised curve
when a translated world coordinate loses precision in the embedded error
estimate. The outer boundary integral retains its error test.

The STEP fixture `parabola_cut_cone.step` gives two independently approximated
degree-seven Bézier pcurves and selects `PCURVE_S1` on its six surface curves.
Integrating the two Bernstein curves on the cone and the cutting plane, and
adding the analytic caps, gives 723.863727303646 mm². The same file with its
six selections set to `CURVE_3D` has the analytic spatial-parabola area
723.863727338942 mm². The regression exercises both representations under
translation and rotation at the original absolute tolerance.

Feature proofs that require a straight or circular boundary query the selected
primitive and its finite parameter interval. Length, points at a fraction of arc
length and nearest-point queries also evaluate a general curve on its carrier.
A cylinder or cone incident on the selected edge supplies a candidate centre
and radius for a circular rim. For positive-weight rational quadratic Bézier
spans, Bernstein coefficients bound radial and axial departure over every span.
The angular derivative must keep its sign through at most one turn. This lets
a double-D bore, a turned chamfer or a recess use a rim written as a spline.
Points, tangents and arc length still come from the original parameterisation.

A planar extrusion copies the selected spatial curve into both caps. A selected
planar rational spline is first restricted to the interval between its STEP
vertices, including a segment that crosses the knot seam. Homogeneous de Boor
evaluation gives the Bézier poles and weights on each retained knot span.
Translating those poles produces an exact tensor-product side wall. For a
quarter-circle cut, the other three quarters of its generating circle must
not seed the wall's ray intersections. For a planar
recess wall adjacent to a cylinder, a straight UV line with a small angular
component defines a helix. Its maximum distance from the segment joining the
STEP vertices is bounded by the larger endpoint mismatch plus `r Δu²/8`.
The wall classifier and the adjoining slot end-wall proof accept it as straight at the kernel coordinate floor
when this bound fits. Area, length and proximity still evaluate the helix.

The two rims of a swept closed spline can start and end at the same knot
value despite traversing a full period. Summing successive parameter changes
on the nearest branch gives the winding used by the boundary integral. A
short trim that crosses the seam and returns has zero net winding. On
`mfcadpp/10155.step` the hole in face #2387 sweeps to a full band:
the side area is its rim length times the sweep height, and independent
OpenCascade extrusion gives a volume of 1987.3563215954 mm³ with zero
intersection with the original solid.

In a periodic face domain, regular ends of selected coedges use their shared
STEP vertex in the parameter branch of each coedge. A winding ray exactly on
a cylindrical seam must see the same finite cut as the two neighbouring
generators. For `mfcadpp/10096.step` face #1989 the file's two planes put
the lower end at −14.83003931168 mm at zero longitude and −19.24211871677 mm
at half a turn. The wall contains −17 mm only on the latter generator.

In `cadgenbench_inputs/cgb243.step.gz`, spherical face `#55337` has radius
2.5 mm. Its three surface curves select `PCURVE_S1`. One is on the sphere and
two are on rational B-spline patches. An independent calculation evaluates
their two-dimensional splines and the two rational tensor products, maps the
lifted points to sphere longitude and latitude, and integrates the Cartesian
moment densities with the spherical Jacobian. Gauss integration with 24 and 48
nodes per span gives 3.541993373916 mm² to twelve decimal places. Selecting
`CURVE_3D` on just those three edges reproduces the independently pinned
spatial area of 3.541992339312 mm². The regression checks both centroids and
all six independent central second moments at their original tolerances.

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
