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
