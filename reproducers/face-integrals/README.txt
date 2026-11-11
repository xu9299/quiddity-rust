Reduced STEP samples for face integration
========================================

I extracted these three trimmed faces from a larger STEP input to reproduce slow
face_mass / face_moments calls. The original part is not included. These samples
retain real geometry, absolute coordinates, complete supporting spline control
nets, and trimming curves. They are not synthetic or geometrically anonymized;
the supporting surfaces can extend beyond the visible trimmed faces.

Files
-----

| File | Original STEP face entity | Bytes | STEP entities |
| --- | --- | ---: | ---: |
| slow_face_a.step | #18820 | 55,134 | 606 |
| slow_face_b.step | #19406 | 55,173 | 606 |
| slow_face_c.step | #19992 | 64,052 | 692 |

Each single-face file contains exactly one ADVANCED_FACE and one supporting
B-spline surface. It is an open surface model, not a solid. Quiddity's face index
is 0 in each single-face file; the original STEP entity identifiers are retained.
Only the three individual STEP files are supplied here, together with this
README, face_reproducer.rs and validation.json. The validation data also retains
metadata from a local combined-file check; that viewing file is not included.

I copied the selected entities and their dependency closure directly from the
source text. All retained geometric entity definitions, their identifiers and
their relative order are unchanged. Only the file header, product label (#7), and
top-level representation item list (#10) were rewritten. Units remain millimetres
and the original declared uncertainty remains 1e-4 mm. No repair, sewing, CAD
re-export, approximation, translation or scaling was applied.

Independent geometry check
--------------------------

OCCT successfully reads the three individual files, reports valid shapes, and
finds one face and no solids in each file. For each
selected finite face, OCCT's area, centroid, and all nine inertia matrix entries
are bit-identical before and after extraction in the same environment.

| Face | OCCT surface area (mm^2) |
| --- | ---: |
| A | 10.759968490101977 |
| B | 10.75950079907782 |
| C | 10.763360456014112 |

These are independent kernel measurements, not rigorous bounds on true integral
error. File hashes, measured timings and diagnostic witnesses are in validation.json.

Reproduce the calls
-------------------

The tested unmodified main baseline is
b590592035733d8ac9fa5b3144bb5d969bda89a0. It includes #232 and #239.
The performance candidate is 29497fc9887761b27365991c72362116234536c3 on that baseline.
The observations use Windows x64, Rust 1.99.1, MSVC and release builds.

Copy face_reproducer.rs into crates/haecceity/tests/face_reproducer.rs in the
upstream source checkout. From its workspace root, using a POSIX shell:

    FACE_REPRO_INPUT=/absolute/path/slow_face_a.step FACE_REPRO_OPERATION=mass \
      cargo test -p haecceity --release --test face_reproducer -- \
        --ignored --exact measure_extracted_faces --nocapture

For PowerShell, set the two environment variables before the same cargo command:

    $env:FACE_REPRO_INPUT = 'C:/absolute/path/slow_face_a.step'
    $env:FACE_REPRO_OPERATION = 'mass'
    cargo test -p haecceity --release --test face_reproducer -- --ignored --exact measure_extracted_faces --nocapture

Repeat with B and C, and use FACE_REPRO_OPERATION=moments in a fresh process for
each cold moment call. The driver logs the completed read time, source entity,
start of each call, result and elapsed call time. It does not set a time limit or
rewrite input. I compiled it against both the baseline and candidate, and ran all
six single-face entry-point cases on the candidate successfully.

Use the files as supplied: repairing or re-exporting them may change the geometry
or parameterization that triggers the problem.

Observed behavior
-----------------

On unmodified main, each single-face mass test process was still evaluating the
call when its 60-second observation window ended. These windows include reading
the small file. The processes were observed concurrently; these are not completed
call times or evidence of mathematical non-integrability.

On the candidate, first calls measured after input reading returned:

| Face | face_mass | Cold face_moments | Result |
| --- | ---: | ---: | --- |
| A | 0.601 s | 0.988 s | None |
| B | 1.168 s | 0.872 s | None |
| C | 1.583 s | 0.855 s | None |

These are observations under concurrent machine load, not repeated isolated
benchmarks. The candidate avoids redundant work; it does not make these integrals
succeed. No speedup ratio is inferred from incomplete baseline observations.

I also checked diagnostic builds of unmodified v0.4.9
(17a883aa1e489c2336b5866baf543c647c31969c) and the pinned main baseline. Across
three faces, two entry points and two versions, all 12 extracted-sample witnesses
exactly match those from the original full input, including floating-point bits.
Instrumentation stops immediately after a decisive accepted-leaf update without
changing the preceding computation. These witnesses are not completed original
calls, and the public timing driver does not contain that instrumentation.

Each witness is at halving depth 12 and fails both the applicable relative-error
test and the fixed whole-face floor. The existing integrator accumulates unresolved
leaf errors with E = max(E, e); subsequent work cannot lower E below that fixed
floor. This establishes refusal by the current bounded algorithm. It does not
prove that another inversion or integration algorithm cannot calculate these faces.

The performance patch and the remaining convergence problem are separate work.
