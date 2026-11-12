# Scheduling independent surface recovery

Full recognition already prepares unqueried freeform faces in a bounded parallel
queue, then consumes the cached results in the caller's original order. Starting
a long face integral late can leave other workers with no remaining queue items.

For queues large enough to use that existing parallel path, recovery now orders
pending faces by a cheap hint from their immutable control-net dimensions and
degrees. The hint is the product of the two control-count-minus-degree values and
the two degree-plus-one values. Saturating integer arithmetic avoids overflow.
It estimates work only: it neither accepts geometry nor removes any query.
Stable ties retain the caller's order. No face identifiers or measured durations
are embedded in the policy.

Only recovery cache preparation is reordered. Face traversal, integrations,
fitting, tolerances, missing-result caching and subsequent result reduction remain
unchanged. Already answered faces and faces outside the requested subset are not
scheduled. The existing four-worker cap, small-queue fallback and thread-stack
policy are unchanged. Other parallel callers retain their existing scheduling.

Directed checks compare the exact recovered plane and all five mass-term bits
against ordered evaluation on differently subdivided trimmed spline faces. They
include a valid saddle that cannot recover to an analytic primitive, a previously
cached refusal, an already recovered face,
an unrequested face, nonfreeform geometry and the small-queue lazy path. These
checks do not establish a whole-part performance gain; complete API comparisons,
strict snapshot validation and bounded real-part measurement remain separate.
