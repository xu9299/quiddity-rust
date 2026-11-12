# Reuse bounded exact samples across fixed-parameter inner integrals

Consecutive inner integrals can evaluate identical nodes while holding the other surface parameter fixed. I retain at most 16,384 finite sample values keyed by exact floating-point bits, clear them when the fixed parameter changes, and avoid allocation on the first visit. Panel weights and every addition retain their original order; panel sums are never subtracted or combined.

The directed regressions compare the affected computations with their original
paths, including refusal and boundary cases. They do not replace full workspace
regression or strict snapshot reproduction on the standard ARM macOS runner.
No expected snapshots, reference bands or inventory floors are changed.
