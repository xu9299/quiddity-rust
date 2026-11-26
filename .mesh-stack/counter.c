/* User-space retired instructions on the thread doing the public STEP workload. */
#include <errno.h>
#include <inttypes.h>
#include <linux/perf_event.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>

int mesh_pmu_begin(void) {
    struct perf_event_attr attr = {0};
    attr.size = sizeof(attr);
    attr.type = PERF_TYPE_HARDWARE;
    attr.config = PERF_COUNT_HW_INSTRUCTIONS;
    attr.disabled = 1;
    attr.exclude_kernel = 1;
    attr.exclude_hv = 1;
    attr.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
    int fd = (int)syscall(SYS_perf_event_open, &attr, 0, -1, -1, 0);
    if (fd == -1) {
        fprintf(stderr, "{\"open_errno\":%d}\n", errno);
        return -1;
    }
    if (ioctl(fd, PERF_EVENT_IOC_RESET, 0) || ioctl(fd, PERF_EVENT_IOC_ENABLE, 0)) {
        perror("enable hardware counter");
        close(fd);
        return -1;
    }
    return fd;
}

int mesh_pmu_end(int fd) {
    uint64_t values[3];
    if (ioctl(fd, PERF_EVENT_IOC_DISABLE, 0) ||
        read(fd, values, sizeof(values)) != (ssize_t)sizeof(values)) {
        perror("read hardware counter");
        close(fd);
        return 1;
    }
    close(fd);
    fprintf(stderr, "{\"instructions\":%" PRIu64 ",\"enabled\":%" PRIu64
                    ",\"running\":%" PRIu64 "}\n", values[0], values[1], values[2]);
    return 0;
}

#ifdef STUDY_CALIBRATION
#include <pthread.h>
#include <stdlib.h>

static void *worker(void *arg) {
    uint64_t n = *(uint64_t *)arg;
    volatile uint64_t value = 1;
    int fd = mesh_pmu_begin();
    if (fd == -1)
        return (void *)(intptr_t)1;
    for (uint64_t i = 0; i < n; ++i)
        value = value * UINT64_C(1664525) + UINT64_C(1013904223);
    int status = mesh_pmu_end(fd);
    printf("%" PRIu64 "\n", value);
    return (void *)(intptr_t)status;
}

int main(int argc, char **argv) {
    if (argc != 2)
        return 2;
    uint64_t n = strtoull(argv[1], NULL, 10);
    pthread_t thread;
    void *result = NULL;
    if (pthread_create(&thread, NULL, worker, &n) || pthread_join(thread, &result))
        return 3;
    return (int)(intptr_t)result;
}
#endif
