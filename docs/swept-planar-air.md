# Planar air along a known sweep

A many-sided prismatic tool can produce hundreds of height breaks when the
volume integral slices across an unrelated world axis. Its known sweep direction
provides a simpler arrangement without replacing, transforming or resampling any
source surface. Rays still query the original finite B-Rep in world coordinates.

`planar_empty_along` is an additional sufficient proof. It requires a positive
finite probe volume, a valid direction, and planar source/probe faces wherever
their conservative boxes meet the overlap region. Curved or unresolved regions
decline this proof. A boundary-free region may still use the existing exact
point-classification shortcut.

The same arrangement breaks and exact planar quadrature are evaluated in an
orthonormal frame along the sweep. Every product and sum in projecting search
boxes rounds outward; boxes only bound the search and never establish air.
The proof uses **zero inset**. Any measured material stops it unsuccessfully;
only a complete integral with zero material and no unanswered required line
accepts. It does not substitute sampling, adjacency, a mesh or a relaxed volume
threshold for the B-Rep proof. The caller retains its ordinary material-fraction
policy for every declined or unsuccessful attempt.

Directed controls include an 80-sided interrupted through-cut, a blind cut,
an obstruction thinner than the ordinary coordinate-floor inset, and generic
rigid motion. Curved and unanswerable boundaries decline the shortcut. Existing
through-cut, native-document, evidence, snapshot and whole-part GUI checks remain
separate gates; no output format changes.

The through-cut caller validates the prospective public record before the air
proof. Its canonical frame, centred section and serialization displacement
constraints are prerequisites regardless of the material result. Moving these
checks does not publish a candidate or weaken its subsequent air proof; it avoids
that expensive proof when the record cannot be represented at all. A manual,
test-only original-STEP diagnostic records the final air decision per component.
