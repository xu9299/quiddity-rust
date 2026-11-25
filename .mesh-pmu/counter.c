/* Count user-space hardware instructions for one child and its worker threads. */
#include <errno.h>
#include <inttypes.h>
#include <linux/perf_event.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <unistd.h>

static void *calibration_worker(void *argument) {
    const uint64_t count = *(uint64_t *)argument;
    volatile uint64_t word = 1;
    for (uint64_t i = 0; i < count; ++i)
        word = word * UINT64_C(1664525) + UINT64_C(1013904223);
    printf("%" PRIu64 "\n", word);
    return NULL;
}

static int count_child(char **command) {
    int ready[2];
    if (pipe(ready)) {
        perror("pipe");
        return 2;
    }
    pid_t child = fork();
    if (child == -1) {
        perror("fork");
        return 2;
    }
    if (child == 0) {
        char byte;
        close(ready[1]);
        if (read(ready[0], &byte, 1) != 1)
            _exit(125);
        close(ready[0]);
        execv(command[0], command);
        perror("execv");
        _exit(126);
    }
    close(ready[0]);
    struct perf_event_attr attr = {0};
    attr.size = sizeof(attr);
    attr.type = PERF_TYPE_HARDWARE;
    attr.config = PERF_COUNT_HW_INSTRUCTIONS;
    attr.disabled = 1;
    attr.enable_on_exec = 1;
    attr.inherit = 1;
    attr.exclude_kernel = 1;
    attr.exclude_hv = 1;
    attr.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
    int fd = (int)syscall(SYS_perf_event_open, &attr, child, -1, -1, 0);
    if (fd == -1) {
        perror("perf_event_open on child");
        close(ready[1]);
        waitpid(child, NULL, 0);
        return 3;
    }
    if (write(ready[1], "1", 1) != 1)
        return 4;
    close(ready[1]);
    int status;
    while (waitpid(child, &status, 0) < 0) {
        if (errno != EINTR)
            return 5;
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        return 6;
    uint64_t data[3];
    if (read(fd, data, sizeof(data)) != (ssize_t)sizeof(data))
        return 7;
    close(fd);
    fprintf(stderr, "{\"instructions\":%" PRIu64 ",\"enabled\":%" PRIu64
                    ",\"running\":%" PRIu64 "}\n", data[0], data[1], data[2]);
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 3 && strcmp(argv[1], "--thread-work") == 0) {
        uint64_t iterations = strtoull(argv[2], NULL, 10);
        pthread_t thread;
        if (pthread_create(&thread, NULL, calibration_worker, &iterations))
            return 8;
        return pthread_join(thread, NULL) != 0;
    }
    if (argc >= 3 && strcmp(argv[1], "--") == 0)
        return count_child(argv + 2);
    return 1;
}
