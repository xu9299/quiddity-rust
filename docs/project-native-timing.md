# Project native timing

The project kernel can report existing recognition and evidence stages through an optional
native Python callback. A session retains its callback only for the enclosing call. Nested
stage guards and a periodic heartbeat identify work that is still running; completed stage
and face records supply durations for diagnosis.

Timing does not select faces, change numerical criteria, omit recognisers or supply geometric
evidence. Without a session, timing guards do not allocate a reporting thread. The stable
0.5.0 entry-point wrapper still translates internal Rust panics into `QuiddityError`, and the
STEP reader retains its complex-reference checks before constructing geometry.

This diagnostic interface is a project patch, separate from the public result schemas.
