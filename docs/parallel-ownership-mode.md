# Ownership mode in parallel floor queries

I bind each floor reader to the ownership mode active when the reader is created.
The serial inventory's locally degraded retry stores this mode in thread-local
state. A new scoped worker starts in strict mode, while the parallel executor
previously propagated only its timing context. A degraded discovery could thus
mix strict worker readings with degraded readings from the caller thread, and its
mode-specific discovery cache could retain that mixed result.

Each bound query sets the captured mode for its invocation and restores the
worker's previous mode afterwards, including during panic unwinding. The same
bound reader serves prepared and lazy floor readings. Strict and degraded cache
slots remain separate. Geometry predicates, exact volume conditions, tolerances,
ordered replay and candidate numbering are unchanged.

The directed checks send a captured reader through scoped worker execution and
verify both modes, the single-worker fallback, nested opposite modes and unwind
restoration. Existing B-Rep, inventory-projection and byte-comparison checks remain
required. These checks do not replace the pinned full workspace and strict
snapshot regression, or complete original-part GUI acceptance.
