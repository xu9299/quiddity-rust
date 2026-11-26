/* Start only on the worker thread, immediately before the public read. */
#include <valgrind/callgrind.h>

int mesh_pmu_begin(void) {
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_ZERO_STATS;
    return 1;
}

int mesh_pmu_end(int fd) {
    if (fd != 1)
        return 1;
    CALLGRIND_DUMP_STATS_AT("mesh-workload");
    CALLGRIND_STOP_INSTRUMENTATION;
    return 0;
}
