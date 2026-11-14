# Face quadrature and parameter directions

For a density `f(u,v)`, Green's theorem permits either `∮ G dv`, with
`G(u,v) = ∫ f(s,v) ds`, or `−∮ H du`, with `H(u,v) = ∫ f(u,s) ds`.
On an eligible nonperiodic, noncollapsed freeform patch with contractible trims,
the shorter spline axis reduces the work spent constructing each antiderivative.
Winding and loop orientation are determined separately from this numerical choice.
Periodic seams, collapsed sides and open parameter walks retain their existing path.
Before changing the primitive, a separate boundary pass checks the necessary closed
differentials `du`, `dv` and `d(uv)`. Each parameter is normalized to its domain span.
A projection-velocity defect may survive a small local embedded-rule error; failing
this consistency check keeps the established integration axis. This check does not
prove uniqueness of a projected foot point or a global quadrature error bound.

Changing the axis changes quadrature nodes and refinement. The error coefficients
remain unchanged; computed magnitudes and absolute floors need not be identical.
Neither direction supplies an independent exact answer for the other.

Quadrature-specific rational evaluation accumulates local homogeneous coordinates relative to an active
control point. The identity `S = O + Σ N w (P−O) / Σ N w` preserves the curve or
surface while reducing cancellation in quotient-rule derivatives. Active controls
that cannot form finite homogeneous coordinates retain the original arithmetic
path instead of masking a nonfinite result. Ordinary curve values, surface normals,
inversion queries and shared UV topology checks retain their established arithmetic.
Local evaluation is used only inside the density and projected-boundary integration
of open, noncollapsed parameter patches. Closed and collapsed frames retain their
ordinary point, derivative, inverse and curve paths consistently. Mixing a local
inverse or tangent frame with their established sheet and winding choices can count
the wrong parameter turn, even when the surface values differ only at round-off.
The eligibility cache contains no source geometry and is excluded from Debug identity.
Control-point transforms clear this cache before the next eligibility query.
The integral-local operand cache retains the same frame for every span pair, including
its finite-control fallback. Cached and direct evaluation must agree bit for bit;
ordinary evaluation retains the upstream operand order and cache behavior.

The projected three-dimensional boundary can enter a held parameter side extremely
close to an edge endpoint. Such a short interval can lie outside both embedded
quadrature rules' outermost nodes. A small embedded-rule difference therefore does
not exclude a missed endpoint transition. Once located, every representable interior
transition on an open, noncollapsed patch is retained; its distance is not compared
with a fixed fraction of the neighboring parameter interval. Closed and collapsed
frames retain their existing near-cut treatment. This does not replace the bounded refinement or
nonfinite rejection rules.

The synthetic unit-square control uses a three-dimensional trim extending slightly
beyond the patch. Its constrained projection still bounds an area of exactly one.
Both loop directions test this analytic result without using confidential geometry.
Additional translated polynomial and rational controls check evaluation derivatives.
