# Project computation reuse on 0.5.0

This branch ports the project kernel's existing computation-reuse changes onto the
official 0.5.0 release. It is a source preparation branch until its directed tests,
native output comparisons, project checks and strict snapshot checks have passed.

The release already includes the face-integral early-refusal work and open-shell
recess evidence handling. Those older patches are not reapplied. The release's
fixture bytes, snapshot provenance, package versions and dependency lock remain
unchanged.

The port preserves the release's joint recognition/evidence API and document
assembly. The document-only project entry shares the same publication path while
omitting only correspondence measurements. Normalisation refusals and inferred
frames pass through the same publisher for full and document-only recognition.
All native Python entry points retain the release's panic translation.

Moment queries retain 0.5.0's independently verified mass components, refinement
criteria and cache publication. The existing integral-local surface evaluator
supplies the same point and partial derivatives to those callbacks; it does not
replace the mass-validity rules. Closed-surface wrapping and tangent-rank checks
also remain those of the release.

The project patches retain exact operand and prefix reuse, shared ray geometry,
section-discovery and evidence-record reuse, ordered floor-proof refusal, and
sequential original-slab rejection witnesses. Additional compute workers and
recovery scheduling are removed; see [serial reuse](serial-project-reuse.md).
No part-specific identifier or sample input is included. Private real-part timings remain outside this source
repository.
