# Prismatic through-cut pair eligibility

A pair can enter a through-cut run group only when both faces have the existing
native-plane normal and their cross product has a unit direction. Previously,
the shared-edge turn was evaluated before checking the second normal or that
direction. A rejected nonplanar neighbour could therefore trigger curve
evaluation and surface-parameter inversion whose result was never used.

The grouping pass now evaluates these existing prerequisites first. It still
visits the same neighbours in the same order, uses the same normal and unit
functions, and requires the exact same concave turn for every eligible pair.
There is no new tolerance or surface classification. Run grouping, ownership,
wall coverage, axial extent and material/air proofs are unchanged. No pair is
accepted solely from its normal, adjacency or a mesh approximation.

Optional timing separates normal collection, pair grouping and group proofs.
It adds no geometry query and does not change published documents. Directed
checks retain the interrupted polygonal cuts, blind refusals, rigid motion,
nearly aligned walls, slot ownership, degenerate normals and curved-trim extent
controls. Complete native output comparisons and bounded whole-part GUI
measurements remain separate from those geometric regressions.

Proof diagnostics further separate ownership, cycle ordering, section creation,
exact bounds, tool construction, wall coverage, intact-core coverage and the air
proof. Component spans retain the first source-face index. These observations do
not alter the proof or its order.

The ignored `profile_original_step_prismatic_components` test reads an explicitly
provided `QUIDDITY_PRISMATIC_PROFILE` STEP path and applies the normal frame policy.
With `QUIDDITY_PRISMATIC_COMPONENTS_ONLY` it reports candidate component sizes and
skips proofs inside test-only code. Otherwise it prepares recovery and runs this
family's proofs, optionally under `QUIDDITY_PRISMATIC_DEGRADED`. It deliberately
does not reconstruct other families' publication ownership and is a diagnostic,
not complete recognition or an acceptance result. These environment hooks are
absent from production builds; the manual test is ignored by normal regressions.
`QUIDDITY_PRISMATIC_SKIP_RECOVERY_WARMUP` leaves the recovery cache cold for a
focused proof diagnosis. Geometry is still read and queried normally; this mode
does not measure the whole-part pipeline or its warmed cache performance.

Volume observations distinguish classifier construction, probe mass, finite-side
preparation, height break generation, source/probe edge piercing, slice breaks
and the two nested integrals. Timers do not add a geometry query, reorder samples,
alter the quadrature rule or change early refusal. A long component air span can
therefore be attributed to setup, ray intersections or integration before any
algorithm is changed.

Periodic observations use a wall-clock deadline independent of completed-span
events. A steady stream of slow slice completions must not hide the active
integral or its accumulated subtotals. The manual family diagnostic includes
volume subtotals in each heartbeat, preserving observations if it is stopped.
