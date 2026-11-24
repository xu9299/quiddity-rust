# Reducing work without changing face-integral decisions

The original work-reuse optimizations retained the integration densities, quadrature nodes, operation order,
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

## Faces with one very large B-spline net

A face whose B-spline net has many knot spans in u and few in v (a threaded part's
helical flank, a sweep along a spline: 300 to 800 spans against one to four), or a
cylinder carrying a helical thread edge, integrates along u for its inner integral at
every boundary node, so the work grew with the spans (or radians) times the nodes (#302).
These changes remove that, none touching the integrands, nodes or the boundary quadrature:

- An open, unpinched face with at least 32 u spans and eight times v's takes the other
  form of Green's theorem, the inner integral along v from the loop's first point, whose
  few spans each cost one 16-point panel. Each panel's rule is compared with the same
  rule on its halves, to 1e-9 of its terms' magnitude in each component (and 1e-12 of the
  largest component's); a face that fails the comparison anywhere, or that this form
  refuses or finds invalid, is integrated again by the stock form, so only a completed,
  checked result is used. A loop's direction still reads the swept parameter area.
  xu9299 independently found the same approach (the other Green's-theorem form) and the
  narrow-dip counterexample in #302.
  xu9299 also found that the stationary-node test must judge motions by the active
  parameter ranges (`domain()`), not the knot vector's end values, which a writer may pad.
- A cylinder or cone whose loops (closed in the lifted parameters, not winding round u) run
  through at least 32 radians of u, as a thread's helical edge does, takes the same form:
  along u the inner integral costs a panel per radian at every boundary node, along v one
  panel. The same check and fallback apply.
- On a face whose inner integral's parameter has many spans (the same threshold: a face
  closed in v, whose loops wind round it, must keep integrating along u), a node that
  moves along its surface's own parameter line, within 1e-12 of its motion across it
  (each motion taken as a share of the surface's active parameter range, `domain()`, so u
  and v of different scales are judged alike), weighs nothing and is not evaluated. This applies to the existing form as
  well as the new one: a closed B-spline face whose loops wind round u and whose v has many
  spans skips such nodes too. Locating the polynomial piece of a node for its
  velocity now searches the knots instead of rebuilding their distinct list (the same
  piece), which had cost one pass over the knots per node.

Results agree with the stock form to round-off (about 1e-10 mm in a recognised location).
`tests/face_moments.rs` and the face-area checks against OpenCascade are unchanged.

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
caller applies the same rejection. In particular, the canonical mass-refinement
caller records the canonical legacy refusal before considering an explicitly
selected master trim. No rejected prefix is published.

This does not reject a parent that can still be subdivided. With ordinary finite
values, an accepted leaf that fails both the relative criterion and `F` was accepted
only because it reached `MAX_DEPTH`. A large error accepted by the relative criterion
does not enter `E`. Equality with `F` is not refusal. With an infinite or NaN floor the
ordered `E > F` guard cannot trigger. Once it does trigger, Rust's `f64::max` cannot
reduce the non-NaN `E`, even if a later estimate is NaN. Evaluation failures remain
failures; no new numerical answer is invented.

## Choosing the STEP boundary for measurement

A STEP surface curve can provide several descriptions of one edge. Its master
representation selects the spatial curve or a numbered pcurve. Public haecceity mass
properties use those selections for face areas, fluxes and moments and their shell and
solid sums. Recognition uses the master boundary only when a selected pcurve deviates
from its 3D edge beyond the face's STEP uncertainty. Otherwise it integrates the
projected 3D boundary in caller coordinates, including under an integration budget.

A pcurve is evaluated in the units and parameter interval of its named carrier.
If that carrier differs from the face being measured, the point and tangent
are transferred through 3D onto the receiving surface. Endpoint agreement,
finite carrier domains, intersection and integration checks govern acceptance.
If the assembled loop fails, measurement retries the whole face using the 3D
edges. Mass and moments share the accepted choice across query orders.

The selected boundary goes through the same integration as the projected one: the
cheaper along-v form and its retry for a face with one very large net, the per-face
integration budget (integrand, deviation sample and selected-boundary gap evaluations spend work),
and the shared moment integration. A remembered refused along-v trial is keyed by the
integrand's width, whether mass is observed and whether the loops are the selected trim,
since the two boundaries are different loops. The choice is parameter-space data (each
selection's pcurve and interval, accepted by distances in millimetres), so a rigid
motion leaves it unchanged: carrying a face's integrals through a motion carries the
boundary choice with its mass, and the moved part's refusals keep their boundary bit.
Where the moment group cannot be integrated against world coordinates (a small selected
trim far from the origin), it is integrated once more against the B-spline net
translated to the face's moment origin, over the same boundary and error test.

## Following the selected boundary in the mesh

The mesher samples each topological edge once for all its incident faces. It
prefers a geometrically qualified selected trim on the carrier face,
then checks the other incident faces for that same selection. Every receiving
face must be able to resolve the selected trim on its own surface. An edge
without compatible receiving trims or a sampleable selection uses
its existing spatial samples. The faces map those common points to their own
parameter charts. On a closed spline, an inverse stopped at the seam also
tries the other end of the chart.
If a selected loop still refuses triangulation after the usual thinning
retries, the mesher retries with spatial edges for that loop's face. Every
incident face uses the same replacement samples on those edges.
This is one compatibility retry, not proof that the selected geometry was
invalid: a qualified selected loop may also exceed the mesher's capabilities.
Its replacement can differ from the boundary used by the face integral. The
mesher does not recompute or change the accepted mass boundary. If the retry
also fails, the existing `MeshRefusal` is returned.

Ordinary boundary selection reads the existing geometric boundary cache without
requesting a face integral or changing the accepted mass-boundary decision.
Each mesh call constructs its own shared edge plans and parameter charts.
Boundary retries rebuild only affected charts and plans; interior refinement
waits until boundary choices settle.

The mesh keeps the shared STEP vertices at the ends of an edge. An unresolved
run of samples next to both the STEP vertex and the selected curve's true end
does not add a second boundary key there.

A selected spherical cap can also be triangulated in a central projection.
All its boundary points must lie well inside one hemisphere and its measured
area must fit in a hemisphere. In that chart, the signed planar area of a
triangle has the same sign as the determinant of the three spatial directions
from the sphere centre. Refinement places points back on the same sphere.
The mesh still reports the surface's usual parameters and original shared
boundary positions.

## Rebuilding affected mesh preparation

A boundary retry invalidates only the faces and shared-edge plans that depend
on that boundary. Interior refinement waits until the shared boundary plans
settle, so a refused selected chart does not discard completed refinement of
every face. A spherical projection is attempted only after the ordinary
parameter chart folds, with the same area and orientation guards.

The dependency state belongs to one meshing call. Retrying a boundary can still
change the triangulation and its refinement cost. This does not change the
general geometry, NURBS evaluators, integration rules or recognition policy.

## Reusing selected-boundary preparation

Mesh qualification and UV loops have separate lazy caches from integration.
Each face retains one selected-edge key and its UV loops. A different key
computes a temporary result instead of growing a map. These caches and the
immutable curve pool belong to the Part.
Loop keys retain selected topological-edge identities. Identical immutable curve
operands share ordinary inversion tables within each Part; parameter units and
every stored scalar bit participate in the key. Intervals, receiving surfaces
and transfers retain their own checks. Curve closure is refreshed after control
point transformations.

Receiving-surface inversions reuse basis operands and nearby parameter hints,
with the existing residual and finite-domain checks. Analytic midpoint bounds
and grouped crossing exclusions keep the original checks for unresolved cases.
The existing recognition boundary policy, metered deviation checks, combined
point/tangent evaluation and carried integral boundaries are preserved. Mesh
and integration interval selection share one resolver, with mesh hints optional.

There is no cross-curve coefficient pool or separate native finite-value
certificate. The integration rule and geometry tolerances are unchanged.

## The native spline integration rule introduced in #281

The original #281 fix tried this specialised integral after a refused projection.
The current selector tries it first when a native, nonperiodic, noncollapsed
B-spline or Bézier face has one complete, explicit master-pcurve loop. The
general boundary integrator handles mixed representations and other face types.
Near parallel surface tangents
can make Gauss–Newton oscillate even with small position residuals; using the file's
selected trim supplies positions and derivatives without that ill-conditioned solve.
Recovered faces measure the region bounded by the **file's master pcurves**. This
can differ from a 3D-edge answer by the file's own representation disagreement.

Each edge needs exactly one master pcurve, finite controls within the surface
domain, and matching endpoints within the file's stated precision plus the kernel
coordinate floor. Recovery additionally samples `|S(pc(t)) - C3d(t)|` at 33 equally
spaced parameters (including endpoints) on **every span in the union of the 2D and
3D curves' knot breaks**. Nonfinite distances or a sampled maximum greater than
`0.001 * L` refuse recovery. Here `L` is the diagonal of this trimmed face's
`Part::face_bounds` axis-aligned box, in millimetres, with no unit-size floor. This
is the bounding-box length normalization used for GVP centroid evidence; tensor
comparisons separately normalize by the reference inertia trace. The interior
check bounds sampled disagreement, not an unsampled mathematical supremum. It
uses the existing GVP band, not STEP's usually much tighter declared precision.
Endpoint and small-gap connectivity checks retain their original precision band.

Recovery also requires `maximum_deviation * P / A <= 0.001`, using that **same
sampled maximum**, the physical length `P` of the master trim on the surface
(including the short UV closing segments), and its integrated area `A`. This
rejects long thin strips whose disagreement is small against the diagonal but
large against their width. Length uses bounded adaptive Gauss–Kronrod with a
face-relative error budget; its summed estimated length error is added to `P`.
The area already passes the recovery integrator's numerical criteria before
this guard runs, and neither mass nor moments are published before it passes.

A known recall limit is that a NURBS pcurve's **whole knot domain must equal the
edge interval**: a valid subinterval trim is still ineligible. Unsupported,
missing, ambiguous or unresolved master trims remain unavailable.

The outer rule retains its coarse `F = 1e-11 * magnitude`, `1e-9` local criterion,
`MAX_DEPTH = 12`, and kink budget. Inner adaptive primitives sum their error
estimates, use the same depth bound, and receive half of F through a conservative
pcurve variation bound. Their explicitly node-paired weighted errors enter each
outer panel's estimate. Unresolved integrals still refuse. These are numerical
error estimates, as in the legacy walk, rather than rigorous true-error bounds.

Canonical mass and boundary choice are published under one once-cell before any
moments escape, including concurrent first queries. An unresolved canonical legacy
quadrature retries only the eligible master boundary, without repeating legacy mass.
A completed but invalid canonical sum stays refused, as in a standalone mass query.
A failure evaluating additional moment fields before canonical mass was determined
still needs an independent canonical query; it cannot poison a valid mass cache.
### Export of source master flags

STEP export preserves the source `PCURVE_S1`/`PCURVE_S2` flag on raw master slots
of native B-spline and Bézier faces, even when no mass query or recovery has taken
place. A source `CURVE_3D` designation stays `CURVE_3D`. This deliberately changes
1,126 faces across the 11 flagged corpus files relative to the previous writer
on main, which always emitted `CURVE_3D` for them. Export copies retained source
metadata without integrating or consulting query-populated caches, preserving
identical bytes with and without prior mass queries.

The synthetic regression uses polynomial normal-offset edges on
`S(u,v)=(w*u+v, epsilon*v, v²)`, with independent one-dimensional references for
area, flux, centroid and central moments. The original narrow cases are rejected by the representation guards; a broad
patch retains a localized ill-conditioned projection near its boundary while
satisfying both bands. Endpoint-matching interior bulges also refuse. Neither reference bands
nor quadrature tolerances are enlarged (ADR 0008).

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

The work-reuse tests retain their independent face-area and moment references,
placement checks, inventory floors and tolerances. Repeated-moment speedups measure
cache reuse; they are not claims about cold integrations or whole-part recognition.


## Selected-mesh resource measurements

Inputs are from `pzfreo/quiddity-corpus` at `3ef06590a7a27ed65ca1f1a6e7f7ea4d5c328f1f`. Each Part is read with `haecceity::read_step_file`. Both timed calls use `part.mesh(deflection, 0.3)`, with `deflection = 0.001 * diagonal` for a finite positive `part.bounds().diagonal()` (otherwise `0.001`). Bounds and deflection setup occur before timing. The first mesh remains alive during the second timed call for exact output comparison; both are checked and released outside timing. The separate memory probe drops each mesh before the next call.

| Input | Corpus-relative file | Paired timing rounds |
| --- | --- | --- |
| cgb245 | `cadgenbench_inputs/cgb245.step` | 30 |
| string_post | `gramel/string_post.step` | 30 |
| cgb243 | `cadgenbench_inputs/cgb243.step.gz` | 30 |
| cgb202 | `cadgenbench_inputs/cgb202.step.gz` | 10 |
| cgb203 | `cadgenbench_inputs/cgb203.step` | 10 |
| cgb207 | `cadgenbench_inputs/cgb207.step` | 10 |
| cgb217 | `cadgenbench_inputs/cgb217.step.gz` | 10 |
| cgb242 | `cadgenbench_inputs/cgb242.step.gz` | 10 |

### Resource comparison: Use STEP-selected boundaries for face meshing

Windows/i7-11800H, Rust 1.99.1/MSVC release, x86-64-v3, one pinned CPU. A single balanced five-version comparison used 30 paired rounds for cgb245, string_post and cgb243, and ten for each other input. Every version used a fresh Part and two mesh calls; import and deflection setup were outside mesh timing, with no warm-up or discarded rounds. Ratios below 1 are faster; intervals are paired 95% bootstrap intervals (5,000 samples, seed 317). The comparison is against merged #317.

| Input | First-call ratio [95% CI] | Repeated-call ratio [95% CI] |
| --- | --- | --- |
| cgb245 | 2.750 [2.684, 2.813] | 1.390 [1.358, 1.417] |
| string_post | 2.167 [2.108, 2.218] | 1.296 [1.268, 1.328] |
| cgb243 | 7.878 [7.793, 7.960] | 4.498 [4.452, 4.539] |
| cgb202 | 2.818 [2.782, 2.853] | 1.071 [1.056, 1.083] |
| cgb203 | 1.077 [1.055, 1.098] | 1.025 [1.009, 1.038] |
| cgb207 | 1.839 [1.726, 1.924] | 1.119 [1.091, 1.144] |
| cgb217 | 1.340 [1.296, 1.383] | 1.032 [0.995, 1.064] |
| cgb242 | 1.929 [1.874, 1.978] | 1.033 [1.014, 1.056] |

| Input | Import ratio [95% CI] | Import plus first mesh ratio [95% CI] |
| --- | --- | --- |
| cgb245 | 1.031 [1.012, 1.049] | 1.481 [1.456, 1.505] |
| string_post | 1.007 [0.981, 1.031] | 1.563 [1.534, 1.592] |
| cgb243 | 1.019 [1.002, 1.038] | 5.350 [5.300, 5.398] |
| cgb202 | 1.019 [1.004, 1.033] | 2.026 [2.006, 2.046] |
| cgb203 | 0.948 [0.840, 1.047] | 1.046 [1.019, 1.072] |
| cgb207 | 0.879 [0.682, 1.024] | 1.493 [1.297, 1.626] |
| cgb217 | 1.022 [0.976, 1.073] | 1.206 [1.167, 1.246] |
| cgb242 | 1.001 [0.944, 1.054] | 1.583 [1.543, 1.622] |

Import plus first mesh excludes bounds/deflection setup, validation and teardown. It checks for work moved into preparation without claiming end-to-end application latency.

Eight-input memory measurements used five fresh processes per version/input and ten calls per process. Median peak private commit changed by +9.0% to +92.4%; cgb243 changed from 140.11 to 221.57 MiB. The correctness fix increases resource use; this is an observed implementation cost, not a lower bound on the cost of correct boundaries. Import plus first mesh is reported separately in the source documentation, including intervals crossing 1. These measurements are not application timings or corpus-wide resource bounds.

| Input | Peak private commit: preceding / current (MiB) | Peak working set: preceding / current (MiB) | Retained Rust increase: preceding / current (MiB) |
| --- | --- | --- | --- |
| cgb245 | 31.5000 / 60.6133 | 36.9414 / 64.6094 | 0.0000 / 17.4228 |
| string_post | 8.4766 / 12.8711 | 14.7578 / 18.9609 | 0.0000 / 2.2241 |
| cgb243 | 140.1094 / 221.5703 | 140.7148 / 219.9688 | 0.0000 / 48.6644 |
| cgb202 | 314.7812 / 427.1523 | 312.0820 / 421.0781 | 0.0000 / 80.3429 |
| cgb203 | 74.1289 / 80.7773 | 78.2891 / 84.7109 | 0.0000 / 5.3766 |
| cgb207 | 48.9102 / 67.7812 | 54.1133 / 72.4453 | 0.0000 / 12.9470 |
| cgb217 | 71.1250 / 103.0117 | 75.1641 / 105.4336 | 0.0000 / 20.9102 |
| cgb242 | 119.1289 / 179.7969 | 122.1211 / 181.6172 | 0.0000 / 43.0728 |

Process peaks cover import and meshing. Private commit is not resident physical memory; working set is reported separately. Rust counts measure requested allocation bytes. Retained allocation is the increase over the post-import, post-deflection-setup state, with the mesh dropped and the Part alive. Every process releases each mesh before the next call and drops the Part at the end. Retained Rust bytes did not grow over ten calls. Valid triangulations need not be bit-identical between stages; geometric regression tests retain their reference bands and fallbacks.

### Resource comparison: Rebuild only affected mesh preparation after boundary fallback

Windows/i7-11800H, Rust 1.99.1/MSVC release, x86-64-v3, one pinned CPU. A single balanced five-version comparison used 30 paired rounds for cgb245, string_post and cgb243, and ten for each other input. Every version used a fresh Part and two mesh calls; import and deflection setup were outside mesh timing, with no warm-up or discarded rounds. Ratios below 1 are faster; intervals are paired 95% bootstrap intervals (5,000 samples, seed 317). The comparison is against the selected-boundary mesh fix.

| Input | First-call ratio [95% CI] | Repeated-call ratio [95% CI] |
| --- | --- | --- |
| cgb245 | 0.924 [0.906, 0.942] | 0.999 [0.976, 1.021] |
| string_post | 1.002 [0.973, 1.034] | 1.001 [0.973, 1.031] |
| cgb243 | 0.662 [0.659, 0.666] | 0.485 [0.480, 0.490] |
| cgb202 | 0.589 [0.576, 0.608] | 1.097 [1.083, 1.112] |
| cgb203 | 0.991 [0.976, 1.005] | 1.009 [0.986, 1.034] |
| cgb207 | 0.941 [0.922, 0.963] | 1.015 [0.990, 1.042] |
| cgb217 | 0.972 [0.950, 0.994] | 0.998 [0.963, 1.035] |
| cgb242 | 0.959 [0.947, 0.973] | 0.995 [0.977, 1.016] |

| Input | Import ratio [95% CI] | Import plus first mesh ratio [95% CI] |
| --- | --- | --- |
| cgb245 | 0.991 [0.975, 1.006] | 0.958 [0.944, 0.972] |
| string_post | 1.011 [0.988, 1.036] | 1.005 [0.982, 1.031] |
| cgb243 | 0.992 [0.973, 1.011] | 0.685 [0.682, 0.689] |
| cgb202 | 0.997 [0.976, 1.015] | 0.679 [0.669, 0.693] |
| cgb203 | 0.985 [0.957, 1.014] | 0.989 [0.978, 1.002] |
| cgb207 | 0.990 [0.977, 1.002] | 0.951 [0.935, 0.970] |
| cgb217 | 0.965 [0.923, 1.008] | 0.969 [0.943, 0.996] |
| cgb242 | 0.993 [0.960, 1.028] | 0.967 [0.951, 0.983] |

Import plus first mesh excludes bounds/deflection setup, validation and teardown. It checks for work moved into preparation without claiming end-to-end application latency.

Eight-input memory measurements used five fresh processes per version/input and ten calls per process. Median peak private commit changed by -1.6% to +8.3%; cgb243 changed from 221.57 to 218.07 MiB. The largest confirmed mesh-time regression was cgb202 (repeated call): 1.097 [1.083, 1.112]. Import plus first mesh is reported separately in the source documentation, including intervals crossing 1. These measurements are not application timings or corpus-wide resource bounds.

| Input | Peak private commit: preceding / current (MiB) | Peak working set: preceding / current (MiB) | Retained Rust increase: preceding / current (MiB) |
| --- | --- | --- | --- |
| cgb245 | 60.6133 / 62.5742 | 64.6094 / 66.9844 | 17.4228 / 17.3707 |
| string_post | 12.8711 / 13.9414 | 18.9609 / 19.8984 | 2.2241 / 2.2241 |
| cgb243 | 221.5703 / 218.0664 | 219.9688 / 214.9180 | 48.6644 / 47.5904 |
| cgb202 | 427.1523 / 449.5352 | 421.0781 / 441.7812 | 80.3429 / 73.6665 |
| cgb203 | 80.7773 / 83.6055 | 84.7109 / 88.3320 | 5.3766 / 5.3766 |
| cgb207 | 67.7812 / 71.1250 | 72.4453 / 75.7461 | 12.9470 / 12.8347 |
| cgb217 | 103.0117 / 108.1914 | 105.4336 / 110.4727 | 20.9102 / 20.9102 |
| cgb242 | 179.7969 / 188.9805 | 181.6172 / 190.0469 | 43.0728 / 42.8993 |

Process peaks cover import and meshing. Private commit is not resident physical memory; working set is reported separately. Rust counts measure requested allocation bytes. Retained allocation is the increase over the post-import, post-deflection-setup state, with the mesh dropped and the Part alive. Every process releases each mesh before the next call and drops the Part at the end. Retained Rust bytes did not grow over ten calls. Valid triangulations need not be bit-identical between stages; geometric regression tests retain their reference bands and fallbacks.

### Resource comparison: Reuse STEP-selected boundary preparation

Windows/i7-11800H, Rust 1.99.1/MSVC release, x86-64-v3, one pinned CPU. A single balanced five-version comparison used 30 paired rounds for cgb245, string_post and cgb243, and ten for each other input. Every version used a fresh Part and two mesh calls; import and deflection setup were outside mesh timing, with no warm-up or discarded rounds. Ratios below 1 are faster; intervals are paired 95% bootstrap intervals (5,000 samples, seed 317). The comparison is against mesh-flow reuse.

| Input | First-call ratio [95% CI] | Repeated-call ratio [95% CI] |
| --- | --- | --- |
| cgb245 | 0.553 [0.540, 0.566] | 0.705 [0.691, 0.718] |
| string_post | 0.576 [0.563, 0.587] | 0.763 [0.744, 0.778] |
| cgb243 | 0.268 [0.266, 0.271] | 0.410 [0.406, 0.414] |
| cgb202 | 0.786 [0.758, 0.831] | 0.874 [0.863, 0.887] |
| cgb203 | 0.885 [0.873, 0.897] | 0.896 [0.879, 0.914] |
| cgb207 | 0.674 [0.657, 0.692] | 0.903 [0.878, 0.924] |
| cgb217 | 0.839 [0.822, 0.859] | 0.885 [0.849, 0.925] |
| cgb242 | 0.606 [0.598, 0.616] | 0.848 [0.814, 0.897] |

| Input | Import ratio [95% CI] | Import plus first mesh ratio [95% CI] |
| --- | --- | --- |
| cgb245 | 1.001 [0.983, 1.020] | 0.791 [0.778, 0.804] |
| string_post | 0.967 [0.945, 0.987] | 0.707 [0.693, 0.720] |
| cgb243 | 0.966 [0.951, 0.979] | 0.339 [0.337, 0.341] |
| cgb202 | 0.996 [0.970, 1.026] | 0.856 [0.834, 0.887] |
| cgb203 | 0.951 [0.921, 0.981] | 0.898 [0.887, 0.910] |
| cgb207 | 1.008 [0.988, 1.027] | 0.744 [0.727, 0.761] |
| cgb217 | 0.975 [0.937, 1.014] | 0.884 [0.870, 0.904] |
| cgb242 | 0.973 [0.951, 0.996] | 0.695 [0.687, 0.703] |

Import plus first mesh excludes bounds/deflection setup, validation and teardown. It checks for work moved into preparation without claiming end-to-end application latency.

Eight-input memory measurements used five fresh processes per version/input and ten calls per process. Median peak private commit changed by -26.1% to +2.5%; cgb243 changed from 218.07 to 192.55 MiB. No mesh-time regression had its entire 95% interval above 1 in these eight inputs; that is not a universal no-regression guarantee. Import plus first mesh is reported separately in the source documentation, including intervals crossing 1. These measurements are not application timings or corpus-wide resource bounds.

| Input | Peak private commit: preceding / current (MiB) | Peak working set: preceding / current (MiB) | Retained Rust increase: preceding / current (MiB) |
| --- | --- | --- | --- |
| cgb245 | 62.5742 / 46.2148 | 66.9844 / 50.9297 | 17.3707 / 9.1376 |
| string_post | 13.9414 / 11.6914 | 19.8984 / 18.0430 | 2.2241 / 1.3409 |
| cgb243 | 218.0664 / 192.5508 | 214.9180 / 189.9453 | 47.5904 / 36.8331 |
| cgb202 | 449.5352 / 400.9414 | 441.7812 / 394.1250 | 73.6665 / 54.5745 |
| cgb203 | 83.6055 / 85.7148 | 88.3320 / 89.1328 | 5.3766 / 6.6449 |
| cgb207 | 71.1250 / 63.1992 | 75.7461 / 68.0391 | 12.8347 / 9.3043 |
| cgb217 | 108.1914 / 91.9062 | 110.4727 / 95.1875 | 20.9102 / 13.5787 |
| cgb242 | 188.9805 / 161.9883 | 190.0469 / 163.5859 | 42.8993 / 30.0853 |

Process peaks cover import and meshing. Private commit is not resident physical memory; working set is reported separately. Rust counts measure requested allocation bytes. Retained allocation is the increase over the post-import, post-deflection-setup state, with the mesh dropped and the Part alive. Every process releases each mesh before the next call and drops the Part at the end. Retained Rust bytes did not grow over ten calls. Valid triangulations need not be bit-identical between stages; geometric regression tests retain their reference bands and fallbacks.



## ARM selected-mesh process measurements

### ARM process resources: pr2

On one standard ARM macOS runner (Rust 1.99.1, release, default features), I measured eight inputs and five fixed versions with `/usr/bin/time -l`: ten balanced rounds per input/call count, 800 fresh processes and 4,400 mesh calls. Each process imports one Part and makes one or ten mesh calls, releasing each mesh before the next call. The counters include startup, import/decompression, bounds setup, meshing, count checks and teardown. Ratios are candidate/comparator; peak resident sizes are medians in MiB.

| Input | Calls | Process time / 6e56fd77 [95% CI] | Process time / preceding [95% CI] | Peak RSS: 6e56fd77 / preceding / current (MiB) |
| --- | --- | --- | --- | --- |
| cgb245 | 1 | 1.396 [1.297, 1.525] | 0.995 [0.937, 1.051] | 41.48 / 68.34 / 70.78 |
| cgb245 | 10 | 1.370 [1.255, 1.549] | 0.986 [0.871, 1.119] | 41.80 / 69.09 / 71.21 |
| string_post | 1 | 1.338 [1.207, 1.457] | 0.948 [0.875, 1.010] | 14.80 / 18.25 / 19.06 |
| string_post | 10 | 1.314 [1.211, 1.447] | 1.049 [0.970, 1.139] | 15.05 / 18.70 / 19.35 |
| cgb243 | 1 | 3.036 [2.929, 3.149] | 0.675 [0.643, 0.703] | 154.05 / 230.62 / 235.91 |
| cgb243 | 10 | 2.305 [2.243, 2.373] | 0.511 [0.498, 0.522] | 154.16 / 232.03 / 235.98 |
| cgb202 | 1 | 1.276 [1.169, 1.386] | 0.694 [0.660, 0.727] | 347.33 / 452.10 / 474.85 |
| cgb202 | 10 | 1.157 [1.082, 1.228] | 0.944 [0.917, 0.971] | 348.51 / 453.20 / 478.45 |
| cgb203 | 1 | 1.062 [1.009, 1.121] | 1.006 [0.973, 1.040] | 87.36 / 93.25 / 97.38 |
| cgb203 | 10 | 1.073 [1.038, 1.112] | 1.010 [0.940, 1.077] | 89.60 / 94.75 / 97.85 |
| cgb207 | 1 | 1.520 [1.438, 1.594] | 0.999 [0.971, 1.027] | 58.22 / 76.34 / 80.51 |
| cgb207 | 10 | 1.313 [1.214, 1.458] | 1.062 [0.971, 1.209] | 58.70 / 76.99 / 80.70 |
| cgb217 | 1 | 1.270 [1.242, 1.301] | 0.963 [0.931, 0.999] | 81.86 / 112.70 / 118.83 |
| cgb217 | 10 | 1.030 [0.976, 1.078] | 0.968 [0.920, 1.014] | 82.21 / 114.21 / 120.58 |
| cgb242 | 1 | 1.423 [1.354, 1.483] | 0.999 [0.942, 1.044] | 146.64 / 204.51 / 215.68 |
| cgb242 | 10 | 1.088 [1.007, 1.163] | 0.953 [0.911, 0.984] | 146.48 / 204.45 / 215.42 |

Across eight inputs, 3 one-call and 3 ten-call workloads had process-time intervals wholly below 1. Median peak RSS changed by +1.7% to +5.6% (largest increase: 25.25 MiB). Supplemental mesh-only timings also showed later-call regressions on string_post (+7.1%) and cgb202 (+4.9%).

`time` prints elapsed seconds to 0.01 s; the smallest one-call workloads are about 0.07–0.10 s. Intervals are pointwise paired bootstrap intervals.

The runner identifies itself as an Apple M1 virtual machine. A native 60-million-iteration probe confirmed that `/usr/bin/time -l` exposes neither retired instructions nor cycles there. These are measured time and peak-memory results, not an instruction-count comparison. The [fixed-source harness](https://github.com/xu9299/quiddity-rust/tree/ffaf56e01889d64333aef1da74933fc75ec29946/.mesh-resource) and [execution recipe](https://github.com/xu9299/quiddity-rust/blob/ffaf56e01889d64333aef1da74933fc75ec29946/.github/workflows/selected-mesh-resources-arm.yml) are available for a physical Mac. Peak physical footprint is recorded separately from peak RSS in the complete report.

[Raw measurement run](https://github.com/xu9299/quiddity-rust/actions/runs/45297463459).

| Input | Calls | Peak footprint: 6e56fd77 / preceding / current (MiB) |
| --- | --- | --- |
| cgb245 | 1 | 37.61 / 64.22 / 66.80 |
| cgb245 | 10 | 37.92 / 64.96 / 67.23 |
| string_post | 1 | 10.88 / 14.25 / 15.06 |
| string_post | 10 | 11.13 / 14.70 / 15.35 |
| cgb243 | 1 | 150.24 / 223.84 / 232.06 |
| cgb243 | 10 | 150.36 / 224.13 / 232.15 |
| cgb202 | 1 | 316.38 / 421.01 / 443.87 |
| cgb202 | 10 | 316.58 / 421.03 / 446.74 |
| cgb203 | 1 | 83.52 / 89.29 / 93.46 |
| cgb203 | 10 | 85.77 / 90.78 / 93.95 |
| cgb207 | 1 | 54.38 / 72.24 / 76.57 |
| cgb207 | 10 | 54.87 / 72.88 / 76.77 |
| cgb217 | 1 | 77.99 / 108.74 / 114.92 |
| cgb217 | 10 | 78.33 / 110.24 / 114.99 |
| cgb242 | 1 | 142.82 / 200.48 / 211.74 |
| cgb242 | 10 | 142.67 / 200.41 / 211.49 |

The same ten-call processes also recorded `Instant` timings inside the probe. These explanatory measurements exclude count checks and mesh release; the last nine calls are averaged within each process before comparing ten paired rounds, rather than treating 90 calls as independent. Import includes decompression. These phase timings do not replace whole-process time or hardware instructions.

| Input | Import ratio [95% CI] | First mesh ratio [95% CI] | Later mesh ratio [95% CI] |
| --- | --- | --- | --- |
| cgb245 | 1.031 [0.988, 1.075] | 0.924 [0.733, 1.142] | 0.989 [0.872, 1.127] |
| string_post | 0.993 [0.908, 1.088] | 1.004 [0.868, 1.170] | 1.071 [1.003, 1.169] |
| cgb243 | 0.921 [0.846, 1.002] | 0.667 [0.620, 0.724] | 0.469 [0.458, 0.480] |
| cgb202 | 0.900 [0.814, 0.983] | 0.603 [0.558, 0.651] | 1.049 [1.009, 1.089] |
| cgb203 | 0.990 [0.900, 1.097] | 0.967 [0.881, 1.055] | 1.016 [0.947, 1.082] |
| cgb207 | 1.079 [0.995, 1.177] | 0.988 [0.939, 1.045] | 1.073 [0.963, 1.257] |
| cgb217 | 0.943 [0.877, 1.020] | 0.925 [0.853, 1.000] | 0.975 [0.920, 1.030] |
| cgb242 | 1.010 [0.879, 1.153] | 0.926 [0.860, 0.986] | 0.955 [0.919, 0.981] |

The common comparison starts at merged #317 (`6e56fd775b43c748430a2bf9dc3f7b561a081822`), whose mesher still uses spatial boundaries. The incremental comparison isolates this change on its direct predecessor. These are bounded input measurements, not corpus-wide resource limits.

### ARM process resources: pr3

On one standard ARM macOS runner (Rust 1.99.1, release, default features), I measured eight inputs and five fixed versions with `/usr/bin/time -l`: ten balanced rounds per input/call count, 800 fresh processes and 4,400 mesh calls. Each process imports one Part and makes one or ten mesh calls, releasing each mesh before the next call. The counters include startup, import/decompression, bounds setup, meshing, count checks and teardown. Ratios are candidate/comparator; peak resident sizes are medians in MiB.

| Input | Calls | Process time / 6e56fd77 [95% CI] | Process time / preceding [95% CI] | Peak RSS: 6e56fd77 / preceding / current (MiB) |
| --- | --- | --- | --- | --- |
| cgb245 | 1 | 1.109 [1.017, 1.199] | 0.794 [0.736, 0.858] | 41.48 / 70.78 / 55.24 |
| cgb245 | 10 | 1.005 [0.944, 1.065] | 0.733 [0.675, 0.791] | 41.80 / 71.21 / 55.54 |
| string_post | 1 | 1.116 [0.987, 1.262] | 0.834 [0.766, 0.937] | 14.80 / 19.06 / 17.14 |
| string_post | 10 | 0.987 [0.944, 1.024] | 0.751 [0.686, 0.815] | 15.05 / 19.35 / 17.45 |
| cgb243 | 1 | 1.155 [1.091, 1.236] | 0.381 [0.357, 0.407] | 154.05 / 235.91 / 214.30 |
| cgb243 | 10 | 0.986 [0.957, 1.018] | 0.428 [0.415, 0.442] | 154.16 / 235.98 / 213.98 |
| cgb202 | 1 | 1.140 [1.041, 1.255] | 0.894 [0.860, 0.927] | 347.33 / 474.85 / 437.59 |
| cgb202 | 10 | 1.054 [0.999, 1.109] | 0.911 [0.878, 0.945] | 348.51 / 478.45 / 440.92 |
| cgb203 | 1 | 0.998 [0.941, 1.052] | 0.940 [0.883, 1.004] | 87.36 / 97.38 / 99.97 |
| cgb203 | 10 | 1.049 [0.986, 1.110] | 0.978 [0.924, 1.037] | 89.60 / 97.85 / 103.33 |
| cgb207 | 1 | 1.119 [1.042, 1.190] | 0.736 [0.710, 0.765] | 58.22 / 80.51 / 73.70 |
| cgb207 | 10 | 1.071 [1.015, 1.127] | 0.816 [0.705, 0.911] | 58.70 / 80.70 / 74.50 |
| cgb217 | 1 | 1.129 [1.085, 1.173] | 0.889 [0.857, 0.920] | 81.86 / 118.83 / 104.45 |
| cgb217 | 10 | 0.971 [0.885, 1.063] | 0.942 [0.877, 1.014] | 82.21 / 120.58 / 104.14 |
| cgb242 | 1 | 1.036 [1.010, 1.067] | 0.728 [0.702, 0.760] | 146.64 / 215.68 / 190.28 |
| cgb242 | 10 | 0.870 [0.801, 0.930] | 0.800 [0.753, 0.836] | 146.48 / 215.42 / 190.33 |

Across eight inputs, 7 one-call and 6 ten-call workloads had process-time intervals wholly below 1. Median peak RSS changed by -22.0% to +5.6% (largest increase: 5.48 MiB). cgb203 accounts for the RSS increase; the other seven inputs use less peak RSS.

`time` prints elapsed seconds to 0.01 s; the smallest one-call workloads are about 0.07–0.10 s. Intervals are pointwise paired bootstrap intervals.

The runner identifies itself as an Apple M1 virtual machine. A native 60-million-iteration probe confirmed that `/usr/bin/time -l` exposes neither retired instructions nor cycles there. These are measured time and peak-memory results, not an instruction-count comparison. The [fixed-source harness](https://github.com/xu9299/quiddity-rust/tree/ffaf56e01889d64333aef1da74933fc75ec29946/.mesh-resource) and [execution recipe](https://github.com/xu9299/quiddity-rust/blob/ffaf56e01889d64333aef1da74933fc75ec29946/.github/workflows/selected-mesh-resources-arm.yml) are available for a physical Mac. Peak physical footprint is recorded separately from peak RSS in the complete report.

[Raw measurement run](https://github.com/xu9299/quiddity-rust/actions/runs/45297463459).

| Input | Calls | Peak footprint: 6e56fd77 / preceding / current (MiB) |
| --- | --- | --- |
| cgb245 | 1 | 37.61 / 66.80 / 51.24 |
| cgb245 | 10 | 37.92 / 67.23 / 51.52 |
| string_post | 1 | 10.88 / 15.06 / 13.11 |
| string_post | 10 | 11.13 / 15.35 / 13.42 |
| cgb243 | 1 | 150.24 / 232.06 / 210.43 |
| cgb243 | 10 | 150.36 / 232.15 / 210.10 |
| cgb202 | 1 | 316.38 / 443.87 / 406.45 |
| cgb202 | 10 | 316.58 / 446.74 / 408.03 |
| cgb203 | 1 | 83.52 / 93.46 / 96.02 |
| cgb203 | 10 | 85.77 / 93.95 / 99.38 |
| cgb207 | 1 | 54.38 / 76.57 / 69.73 |
| cgb207 | 10 | 54.87 / 76.77 / 70.53 |
| cgb217 | 1 | 77.99 / 114.92 / 100.49 |
| cgb217 | 10 | 78.33 / 114.99 / 100.17 |
| cgb242 | 1 | 142.82 / 211.74 / 186.26 |
| cgb242 | 10 | 142.67 / 211.49 / 186.31 |

The same ten-call processes also recorded `Instant` timings inside the probe. These explanatory measurements exclude count checks and mesh release; the last nine calls are averaged within each process before comparing ten paired rounds, rather than treating 90 calls as independent. Import includes decompression. These phase timings do not replace whole-process time or hardware instructions.

| Input | Import ratio [95% CI] | First mesh ratio [95% CI] | Later mesh ratio [95% CI] |
| --- | --- | --- | --- |
| cgb245 | 1.008 [0.931, 1.088] | 0.522 [0.453, 0.602] | 0.712 [0.651, 0.770] |
| string_post | 0.970 [0.904, 1.037] | 0.604 [0.546, 0.665] | 0.747 [0.680, 0.808] |
| cgb243 | 0.947 [0.862, 1.045] | 0.276 [0.255, 0.298] | 0.444 [0.434, 0.456] |
| cgb202 | 1.024 [0.940, 1.115] | 0.725 [0.687, 0.769] | 0.941 [0.907, 0.973] |
| cgb203 | 0.960 [0.908, 1.009] | 1.013 [0.944, 1.088] | 0.977 [0.920, 1.039] |
| cgb207 | 0.839 [0.748, 0.937] | 0.629 [0.579, 0.677] | 0.850 [0.716, 0.964] |
| cgb217 | 1.044 [0.954, 1.145] | 0.862 [0.797, 0.936] | 0.949 [0.878, 1.024] |
| cgb242 | 0.940 [0.852, 1.033] | 0.609 [0.571, 0.645] | 0.831 [0.781, 0.868] |

The common comparison starts at merged #317 (`6e56fd775b43c748430a2bf9dc3f7b561a081822`), whose mesher still uses spatial boundaries. The incremental comparison isolates this change on its direct predecessor. These are bounded input measurements, not corpus-wide resource limits.
