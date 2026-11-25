/* Check whether this Linux runner exposes a usable user-mode retired-instruction counter. */
#include <errno.h>
#include <inttypes.h>
#include <linux/perf_event.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>

int main(int argc, char **argv) {
    if (argc != 2)
        return 2;
    const uint64_t iterations = strtoull(argv[1], NULL, 10);
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
        int error = errno;
        printf("{\"iterations\":%" PRIu64 ",\"open_errno\":%d}\n", iterations, error);
        return 0;
    }
    if (ioctl(fd, PERF_EVENT_IOC_RESET, 0) || ioctl(fd, PERF_EVENT_IOC_ENABLE, 0))
        return 3;
    volatile uint64_t word = 1;
    for (uint64_t i = 0; i < iterations; ++i)
        word = word * UINT64_C(1664525) + UINT64_C(1013904223);
    if (ioctl(fd, PERF_EVENT_IOC_DISABLE, 0))
        return 4;
    uint64_t data[3];
    if (read(fd, data, sizeof(data)) != sizeof(data))
        return 5;
    close(fd);
    printf("{\"iterations\":%" PRIu64 ",\"value\":%" PRIu64
           ",\"enabled\":%" PRIu64 ",\"running\":%" PRIu64 ",\"checksum\":%" PRIu64 "}\n",
           iterations, data[0], data[1], data[2], word);
    return 0;
}
