/* Optional Windows benchmark instrumentation; not linked into quiddity.
 * Count executed application instructions in one selected DLL, between explicit
 * markers. Private code caches make each embedded counter pointer thread-local.
 */
#include "dr_api.h"
#include "drmgr.h"
#include "drwrap.h"
#include <intrin.h>

#pragma intrinsic(_InterlockedExchangeAdd64, _InterlockedCompareExchange64)

typedef struct thread_count {
    volatile __int64 instructions;
    uint64 baseline;
    struct thread_count *next;
} thread_count;

static int tls_index;
static thread_count *threads;
static void *lock;
static app_pc module_start, module_end;
static app_pc marker_pc;
static const char *module_name;
static file_t output;
static bool active;
static char region[128];
static uint regions;

static uint64 read_count(thread_count *state)
{
    return (uint64)_InterlockedCompareExchange64(&state->instructions, 0, 0);
}

static void fail(const char *message)
{
    dr_fprintf(STDERR, "instruction_counter: %s\n", message);
    dr_abort_with_code(2);
}

static bool equal(const char *a, const char *b)
{
    while (*a && *a == *b) { ++a; ++b; }
    return *a == *b;
}

static void thread_init(void *drcontext)
{
    thread_count *state = dr_global_alloc(sizeof(*state));
    state->instructions = 0;
    state->baseline = 0;
    dr_mutex_lock(lock);
    state->next = threads;
    threads = state;
    dr_mutex_unlock(lock);
    drmgr_set_tls_field(drcontext, tls_index, state);
}

static void marker(void *wrapcxt, void **unused)
{
    char label[128];
    const char *source = drwrap_get_arg(wrapcxt, 1);
    ptr_uint_t operation = (ptr_uint_t)drwrap_get_arg(wrapcxt, 0);
    bool begin = operation == 1;
    size_t bytes;
    (void)unused;
    if (operation > 1) fail("invalid marker operation");
    for (size_t n = 0; ; ++n) {
        if (n == sizeof(label)) fail("marker label too long");
        if (!dr_safe_read(source + n, 1, &label[n], &bytes) || bytes != 1)
            fail("cannot read marker label");
        if (label[n] == 0) break;
    }
    if (module_start == NULL) fail("target DLL has not loaded");
    dr_mutex_lock(lock);
    if (begin) {
        size_t n = 0;
        if (active) fail("nested measurement");
        while (label[n]) {
            char c = label[n];
            if (n >= sizeof(region) - 1 || !((c >= 'a' && c <= 'z') ||
                (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') ||
                c == '_' || c == '-' || c == '.')) fail("invalid region label");
            region[n] = c;
            ++n;
        }
        region[n] = 0;
        for (thread_count *t = threads; t; t = t->next) t->baseline = read_count(t);
        active = true;
    } else {
        uint64 total = 0;
        uint contributing = 0;
        if (!active || !equal(region, label)) {
            dr_fprintf(STDERR, "instruction_counter: end=%s active=%d region=%s\n",
                       label, active, region);
            fail("unmatched measurement end");
        }
        for (thread_count *t = threads; t; t = t->next) {
            uint64 delta = read_count(t) - t->baseline;
            total += delta;
            if (delta) ++contributing;
        }
        dr_fprintf(output, "{\"region\":\"%s\",\"instructions\":%llu,\"threads\":%u}\n",
                   region, total, contributing);
        ++regions;
        active = false;
    }
    dr_mutex_unlock(lock);
    /* Let the no-op return normally; its address is excluded from instrumentation. */
}

static void module_load(void *drcontext, const module_data_t *info, bool loaded)
{
    const char *name = info->full_path;
    for (const char *p = name; p && *p; ++p) {
        if (*p == '/' || *p == '\\') name = p + 1;
    }
    (void)drcontext;
    (void)loaded;
    if (name && equal(name, module_name)) {
        if (module_start != NULL) fail("multiple target modules");
        module_start = info->start;
        module_end = info->end;
    }
    if (name && equal(name, "icount_calibration.dll")) {
        marker_pc = (app_pc)dr_get_proc_address(info->handle, "icount_marker");
        if (!marker_pc || !drwrap_wrap(marker_pc, marker, NULL)) fail("cannot wrap marker");
    }
}

/* This short clean call is eligible for DynamoRIO's inline optimization. */
static void count_block(volatile __int64 *counter, uint amount)
{
    /* A boundary may sample another thread while its final worker code is running. */
    _InterlockedExchangeAdd64(counter, amount);
}

static dr_emit_flags_t analyse(void *drcontext, void *tag, instrlist_t *bb,
                              bool for_trace, bool translating, void **data)
{
    instr_t *first = instrlist_first_app(bb);
    app_pc pc = first ? instr_get_app_pc(first) : NULL;
    uint count = 0;
    (void)drcontext; (void)tag; (void)for_trace; (void)translating;
    if (module_start && pc >= module_start && pc < module_end && pc != marker_pc) {
        for (instr_t *i = first; i; i = instr_get_next_app(i)) ++count;
    }
    *data = (void *)(ptr_uint_t)count;
    return DR_EMIT_DEFAULT;
}

static dr_emit_flags_t instrument(void *drcontext, void *tag, instrlist_t *bb,
                                 instr_t *where, bool for_trace, bool translating, void *data)
{
    (void)tag; (void)for_trace; (void)translating;
    if (data && drmgr_is_first_instr(drcontext, where)) {
        thread_count *state = drmgr_get_tls_field(drcontext, tls_index);
        if (!state) fail("missing thread counter");
        dr_insert_clean_call(drcontext, bb, where, (void *)count_block, false, 2,
                             OPND_CREATE_INTPTR(&state->instructions),
                             OPND_CREATE_INT32((uint)(ptr_uint_t)data));
    }
    return DR_EMIT_DEFAULT;
}

static void finish(void)
{
    if (active) fail("unterminated measurement");
    dr_fprintf(output, "{\"complete\":true,\"regions\":%u}\n", regions);
    dr_close_file(output);
    drwrap_exit();
    drmgr_exit();
    while (threads) {
        thread_count *next = threads->next;
        dr_global_free(threads, sizeof(*threads));
        threads = next;
    }
    dr_mutex_destroy(lock);
}

DR_EXPORT void dr_client_main(client_id_t id, int argc, const char *argv[])
{
    (void)id;
    if (argc != 3) fail("usage: client.dll OUTPUT.jsonl MODULE.dll");
    if (!dr_using_all_private_caches()) fail("requires drrun -thread_private");
    module_name = argv[2];
    output = dr_open_file(argv[1], DR_FILE_WRITE_REQUIRE_NEW);
    if (output == INVALID_FILE) fail("output must be a new file");
    if (!drmgr_init() || !drwrap_init()) fail("extension initialization failed");
    lock = dr_mutex_create();
    tls_index = drmgr_register_tls_field();
    if (tls_index < 0 || !drmgr_register_thread_init_event(thread_init) ||
        !drmgr_register_module_load_event(module_load) ||
        !drmgr_register_bb_instrumentation_event(analyse, instrument, NULL))
        fail("event registration failed");
    if (!drmgr_register_exit_event(finish)) fail("exit event registration failed");
}
