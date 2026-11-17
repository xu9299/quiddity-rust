# Finite support coverage

`haecceity::cover::assess` reports covered, uncovered, or indeterminate. A failed proof is
indeterminate, rather than evidence of a physical gap. The boolean convenience accepts only
established coverage. Neither path returns to the former directional sampling margin.

## Planar straight arrangements

Finite feet on the analytic lines form the nominal boundary, joined through each shared
topological vertex. Local endpoint boxes charge the uncertain connection. Independent file
pcurves contribute the area swept between their finite segments and the nominal segments.
No infinite line intersection is substituted for an edge endpoint. Crossings partition the
polygon arrangement into bands where uncovered width is affine. Its integral and the trim
uncertainty bracket uncovered area against the unchanged relative budget, `1e-9 * area`.

A curved edge replaced by its vertex chord can change parity only inside a conservative box
containing the entire finite curve and that chord. Analytic extrema bound conics; positive
rational weights allow the control hull to bound a spline. Independent pcurves must also be
bounded. An unsupported representation makes this path indeterminate. File uncertainty only
guards endpoint coherence; it is not an acceptance margin.
For two finite straight representations, a straight homotopy has an affine Jacobian
determinant on the parameter square. Convexity bounds its absolute integral by the mean
of the four corner determinants, which bounds swept area even for folds. Roundoff is
added, and the minimum with the independently enclosing capsule and clipped-box bounds
is retained. This preserves transverse displacement instead of charging a longitudinal
slide as a full-width strip. It changes only an error bound, not the relative area budget.
Constructed parts without a declared modelling precision retain the vertex path and analytic
segment together in the area bracket. Missing metadata is not interpreted as zero precision
and does not permit silently replacing an inconsistent boundary.

## Corresponding complete boundaries

This additional path proves acceptance only. Target and support need a bijection of their
loops, including holes, and a cyclic edge pairing in either orientation. A conic is written
as `o + x*cos(t) + y*sin(t)` over its finite interval. Absorbing its start angle and sense
into the basis bounds equal fractions of two intervals by the sum of the origin and basis
differences plus a length coefficient times their span difference. Full closed conics may
align their seam phase; partial arcs must retain their actual span and sense.

If corresponding points differ by at most `d`, their straight homotopy remains in the
`d`-tube of the reference curve. A rectifiable curve of length at most `L` has tube area at
most `2*d*L + pi*d*d`. The parity regions can differ only in these tubes and the separately
charged endpoint connection boxes. Independent pcurves are included even when they are not
the declared master. Unevaluable pcurves cannot disappear from the proof.

For a cylinder, physical parameter coordinates are `(radius*u, v)`. A lower radial bound
keeps each homotopy away from the axis. At radial distance at least half the cylinder's
radius the parameter map has bounded derivative. Its Jacobian has orthogonal rows with norms
`R/rho` and `1`. If `r` is the reference curve's lower radial bound, displacement uses the
factor `max(1,R/(r-d))`, and reference length uses `max(1,R/r)`, each rounded conservatively.
This uses the actual guarded radial interval rather than fixed worst-case factors. Curves
failing those guards remain indeterminate. Finite coaxial conics and axial lines give the target's unwrapped
parameter polygon. A loop with unresolved winding is refused.

The denominator has a separate conservative lower area bound. Coarse conic chords or ideal
cylindrical parameter segments give a polygon; their swept-area bounds, representation
uncertainty and numerical allowance are subtracted from its area. This polygon never proves
support coverage. Acceptance requires the complete boundary correspondence error to fit
`1e-9` times that lower bound. A nonpositive lower bound refuses the proof.

This path does not establish arbitrary unions of curved support regions, compare unsupported
spline boundaries, or infer a gap when no correspondence is found. It makes no changes to
the input geometry or to saved recognition evidence.

Two-pole degree-one positive rational pcurves retain their finite straight segment, including
endpoint clamping. Axial cylindrical pcurves of this form map to the same finite segment in
3D. Other positive rational splines are split at their knots into homogeneous Bezier nets,
then restricted to the finite interval. Periodic wraps retain both pieces of the domain.
Subdivision tightens the positive rational control hulls. Endpoint observations select which
extrema need more work; they never replace the enclosing hulls. Exhaustion retains every
remaining hull. Their combined box encloses the finite curve and its vertex chord, without
charging unrelated spans. Invalid intervals or unsuccessful extraction retain the full hull.
Temporary knot insertion and subdivision operate on copied control data; the input B-Rep
remains unchanged.

For a finite positive pcurve on a native B-spline or Bezier surface, the complete UV control
hull bounds the parameter domain. Inside the original domain, the surface's positive control
net bounds its entire image. A small excursion past a clamped domain end uses the end
polynomial's Bernstein extension: its basis total variation is bounded by
`(1+2*excursion/span)^degree`. The tensor product gives a bound on negative coefficient mass,
and the positive weight range bounds the rational denominator away from zero. The complete
control box is expanded by that normalized negative mass times each coordinate range.
This is a geometric enclosure of the extension, not a precision or acceptance margin.
An excessive excursion or a nonpositive denominator refuses the bound.
A planar support may charge the intersection of that full projected
box, including the edge's vertices, with the queried region. A proved remote hull has zero
local area charge. An intersecting hull remains charged even if the nominal 3D edge is
remote. This bound does not assume an exact curve composition, use samples, extend an
unbounded parameter line, or reuse a converted swept surface's different parameter units.

## Sufficient containment

A straight convex planar support admits a second acceptance proof. Its vertex polygon is
checked for strict convexity, and every finite straight representation stays within a
bounded displacement of that polygon. Every target edge, including complete conics,
positive rational control hulls and independent pcurves, must lie strictly inside all of
the support's half-planes after that displacement and roundoff allowance. A supporting face
with additional holes cannot use this proof. An interior seam vertex alone is insufficient.

A cylindrical face with two closed winding-one contours and one shared seam bounds a full
band. Opposite contour windings and a closed topological walk are required. Finite curve
extrema bound each contour's axial range; the conservative interior band excludes both
ranges. Independent seam representations contribute a physical parameter tube area.
For a target's additional loops, positive control hulls and all retained pcurves must stay
in an open radial hemisphere and within the outer band. The hemisphere establishes zero
winding. The complete circumference times each hole's axial span bounds the removed area,
including non-linear cylindrical pcurves whose UV control hull is positive and finite.
The resulting lower target area and upper missing band area use the unchanged relative
budget. A shorter support, an unresolved boundary or a hole outside the band cannot be
accepted by this path. Neither containment proof establishes rejection.
