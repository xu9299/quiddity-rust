# NURBS evaluation within face integrals

The face integrator repeatedly evaluates one parameter while the other stays fixed.
`nurbs::integral::IntegralEvaluation` specializes that access pattern. Its state belongs
to one immutable surface, one integral, and one integration direction. Surface inversion,
ray intersection, boundary projection, and public point evaluation retain their existing
evaluators. No additional computation thread or shared cache is introduced.

## Fixed-parameter contraction

Write the homogeneous surface as `H(u,v) = sum_i sum_j N_i(u) M_j(v) H_ij`, where
`H_ij = (w_ij P_ij, w_ij)`. For a fixed `v`, prepare

```
Q_i = sum_j M_j(v)  H_ij
R_i = sum_j M'_j(v) H_ij
H   = sum_i N_i(u)  Q_i
H_u = sum_i N'_i(u) Q_i
H_v = sum_i N_i(u)  R_i
```

The roles of `u` and `v` reverse for the other integration direction. Retaining `R_i`
is necessary: an isocurve tangent alone does not supply both partials of the surface.
Cartesian position and both partials follow the usual rational quotient rule.
The fixed parameter's exact clamped bits and both knot spans key the contraction;
moving either span or the fixed parameter rebuilds it.

## Centered span polynomials

On each nonzero knot interval, the B-spline basis is polynomial. The basis recurrence
builds coefficients in the local coordinate `x = (t-centre)/half`; Horner evaluates
the basis and its derivative, scaled by `1/half`. The original knot sequence and
span-selection convention remain authoritative, including repeated knots and endpoints.
This is an algebraic change of basis, not fitting or simplifying the input geometry.

The specialized path supports degrees zero through seven and uses eight-entry operands.
Each axis has at most 2 MiB of polynomial slots, including empty slots at repeated knots.
An axis exceeding that bound retains basis recurrence. Zero-width, very short, poorly
resolved, or nonfinite polynomial intervals also retain recurrence. Higher degrees use
the original integral evaluator. Nonfinite clamped parameters, weighted controls, or specialized
results fall back to the original evaluator. The parameter direction is immutable for
an evaluator's lifetime.

An active control net with a largest-to-smallest weight ratio above 256 also uses the
original evaluator. Near a span endpoint, a tiny absolute error in a nearly zero basis
value can otherwise be amplified by a large rational weight, even while every value
stays finite. The cutoff is a conservative applicability guard, not a uniform error
bound. A flat rational sweep with a closed-form parameterization tests this failure
mode, including weight contrasts through `1e16`, endpoints, and both directions.

The floating-point association changes. Mathematically equal formulas need not produce
bit-identical measurements or identical adaptive subdivision decisions. Existing integral
rules, error thresholds, and refusal conditions are unchanged; those rules still decide
whether a measurement is available. Numerical checks and snapshot differences must be
reviewed separately from timing results.

## Validation and reproducible timing

Directed checks cover both partials and integration directions, changing spans and fixed
parameters, repeated knots, clamping, degrees through the original limit, bounded memory,
ill-resolved intervals, and nonfinite fallback. An independent composite Simpson integral
checks the analytic area, position flux, and area vector of a rational quarter-cylinder.
These examples do not establish a uniform error bound for arbitrary input geometry.

The ignored `nurbs::integral::tests::benchmark_synthetic_inner_sweeps` check constructs
its own rational spline control grids and compares the original and specialized evaluators
in the same executable. It brackets candidate timings with original timings, prints raw
seconds and checksums, and makes no timing assertion. Run it serially in a release build;
report the compiler, CPU flags, measured source, and machine contention alongside results.
It measures repeated surface evaluations, not a complete recognition job or batch throughput.

The algebra is described in the [nurbspy isocurve documentation](https://turbo-sim.github.io/nurbspy/theory/surface_isocurves.html).
Local power-basis evaluation is also discussed by
[Gatilov, Computer-Aided Design (2016)](https://doi.org/10.1016/j.cad.2015.10.006).
Speedups reported in that paper are not measurements of this implementation.
