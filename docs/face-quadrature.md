# Face quadrature and parameter directions

For a density `f(u,v)`, Green's theorem permits either `∮ G dv`, with
`G(u,v) = ∫ f(s,v) ds`, or `−∮ H du`, with `H(u,v) = ∫ f(u,s) ds`.
On an eligible nonperiodic, noncollapsed freeform patch with contractible trims,
the shorter spline axis reduces the work spent constructing each antiderivative.
Winding and loop orientation are determined separately from this numerical choice.
Periodic seams, collapsed sides and open parameter walks retain their existing path.

Changing the axis changes quadrature nodes and refinement. The error coefficients
remain unchanged; computed magnitudes and absolute floors need not be identical.
Neither direction supplies an independent exact answer for the other.

Rational evaluation accumulates local homogeneous coordinates relative to an active
control point. The identity `S = O + Σ N w (P−O) / Σ N w` preserves the curve or
surface while reducing cancellation in quotient-rule derivatives. Active controls
that cannot form finite homogeneous coordinates retain the original arithmetic
path instead of masking a nonfinite result.

The projected three-dimensional boundary can enter a held parameter side extremely
close to an edge endpoint. Such a short interval can lie outside both embedded
quadrature rules' outermost nodes. A small embedded-rule difference therefore does
not exclude a missed endpoint transition. Once located, every representable interior
transition is retained; its distance is not compared with a fixed fraction of the
neighboring parameter interval. This does not replace the bounded refinement or
nonfinite rejection rules.

The synthetic unit-square control uses a three-dimensional trim extending slightly
beyond the patch. Its constrained projection still bounds an area of exactly one.
Both loop directions test this analytic result without using confidential geometry.
Additional translated polynomial and rational controls check evaluation derivatives.
