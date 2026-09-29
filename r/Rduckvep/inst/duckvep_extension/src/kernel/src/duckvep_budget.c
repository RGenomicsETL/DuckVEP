#include "duckvep_budget.h"

#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
#define DUCKVEP_TLS __declspec(thread)
#else
#define DUCKVEP_TLS _Thread_local
#endif

#define HEADER_BYTES ((size_t)16)
#define HEADER_MAGIC UINT32_C(0xB0D6E7A1)
#define PAGE_BYTES ((size_t)4096)
#define MMAP_THRESHOLD ((size_t)131072)

typedef struct header {
    uint64_t charged;
    uint32_t owner;
    uint32_t magic;
} header_t;

static _Atomic uint64_t g_limit = DUCKVEP_BUDGET_DEFAULT_BYTES;
static _Atomic uint64_t g_total;
static _Atomic uint64_t g_total_high;
static _Atomic uint64_t g_refusals;
static _Atomic uint64_t g_charges;
static _Atomic uint64_t g_current[DUCKVEP_OWNER_COUNT];
static _Atomic uint64_t g_high[DUCKVEP_OWNER_COUNT];
static _Atomic uint32_t g_workers = DUCKVEP_BUDGET_DEFAULT_WORKERS;
static _Atomic uint64_t g_scratch = DUCKVEP_BUDGET_DEFAULT_SCRATCH_BYTES;
static _Atomic uint64_t g_emit = DUCKVEP_BUDGET_DEFAULT_EMIT_BYTES;
static _Atomic uint64_t g_idle = DUCKVEP_BUDGET_DEFAULT_IDLE_BYTES;

static DUCKVEP_TLS struct {
    int pending;
    char what[64];
    uint64_t requested, in_use, limit;
} t_failure;

#ifdef DUCKVEP_FAULT_INJECTION
#if defined(__SANITIZE_ADDRESS__)
void __sanitizer_print_stack_trace(void);
#endif
static _Atomic uint64_t g_fault_count;
static _Atomic uint64_t g_fault_target;
static _Atomic uint64_t g_fault_fired;
#endif

const char *
duckvep_budget_owner_name(duckvep_budget_owner_t owner)
{
    static const char *const names[DUCKVEP_OWNER_COUNT] = {
        "model", "index", "reference", "workspace", "scratch", "emit",
        "control"
    };
    return (unsigned)owner < DUCKVEP_OWNER_COUNT ? names[owner] : "unknown";
}

size_t
duckvep_budget_charge_for(size_t size)
{
    size_t total;

    if (size > SIZE_MAX - HEADER_BYTES - 2 * PAGE_BYTES)
        return SIZE_MAX;
    total = size + HEADER_BYTES;
    if (total >= MMAP_THRESHOLD)
        return (total + PAGE_BYTES - 1) & ~(PAGE_BYTES - 1);
    /* glibc chunk: 8 bytes of bookkeeping, 16-byte granularity. */
    return (total + 8 + 15) & ~(size_t)15;
}

static void
raise_high(_Atomic uint64_t *high, uint64_t value)
{
    uint64_t seen = atomic_load_explicit(high, memory_order_relaxed);

    while (value > seen &&
        !atomic_compare_exchange_weak_explicit(high, &seen, value,
        memory_order_relaxed, memory_order_relaxed))
        ;
}

void
duckvep_budget_note_failure(const char *what, uint64_t requested,
    uint64_t in_use, uint64_t limit)
{
    t_failure.pending = 1;
    (void)snprintf(t_failure.what, sizeof t_failure.what, "%s",
        what != NULL ? what : "native memory");
    t_failure.requested = requested;
    t_failure.in_use = in_use;
    t_failure.limit = limit;
}

int
duckvep_budget_take_failure(char *buffer, size_t size)
{
    if (!t_failure.pending)
        return 0;
    if (buffer != NULL && size != 0)
        (void)snprintf(buffer, size,
            "capacity error: %s budget exceeded (requested %llu bytes, "
            "%llu in use, limit %llu)", t_failure.what,
            (unsigned long long)t_failure.requested,
            (unsigned long long)t_failure.in_use,
            (unsigned long long)t_failure.limit);
    t_failure.pending = 0;
    return 1;
}

void
duckvep_budget_clear_failure(void)
{
    t_failure.pending = 0;
}

/* Charge before allocating. */
static int
charge(duckvep_budget_owner_t owner, uint64_t bytes)
{
    uint64_t limit, used;

    limit = atomic_load_explicit(&g_limit, memory_order_relaxed);
    used = atomic_load_explicit(&g_total, memory_order_relaxed);
    for (;;) {
        if (bytes > limit || used > limit - bytes) {
            char what[64];

            atomic_fetch_add_explicit(&g_refusals, 1, memory_order_relaxed);
            (void)snprintf(what, sizeof what, "native memory (%s)",
                duckvep_budget_owner_name(owner));
            duckvep_budget_note_failure(what, bytes, used, limit);
            return 0;
        }
        if (atomic_compare_exchange_weak_explicit(&g_total, &used,
            used + bytes, memory_order_relaxed, memory_order_relaxed))
            break;
    }
    atomic_fetch_add_explicit(&g_charges, 1, memory_order_relaxed);
    raise_high(&g_total_high, used + bytes);
    raise_high(&g_high[owner],
        atomic_fetch_add_explicit(&g_current[owner], bytes,
        memory_order_relaxed) + bytes);
    return 1;
}

static void
uncharge(duckvep_budget_owner_t owner, uint64_t bytes)
{
    atomic_fetch_sub_explicit(&g_current[owner], bytes, memory_order_relaxed);
    atomic_fetch_sub_explicit(&g_total, bytes, memory_order_relaxed);
}

#ifdef DUCKVEP_FAULT_INJECTION
static int
fault_fires(void)
{
    uint64_t index, target;

    index = atomic_fetch_add_explicit(&g_fault_count, 1,
        memory_order_relaxed) + 1;
    target = atomic_load_explicit(&g_fault_target, memory_order_relaxed);
    if (target != 0 && index == target) {
        atomic_fetch_add_explicit(&g_fault_fired, 1, memory_order_relaxed);
        if (getenv("DUCKVEP_FAULT_TRACE") != NULL) {
            (void)fprintf(stderr, "DUCKVEP_FAULT_FIRED index=%llu\n",
                (unsigned long long)index);
#if defined(__SANITIZE_ADDRESS__)
            __sanitizer_print_stack_trace();
#endif
            (void)fflush(stderr);
        }
        return 1;
    }
    return 0;
}
void duckvep_fault_arm(uint64_t nth)
{
    atomic_store(&g_fault_count, 0);
    atomic_store(&g_fault_fired, 0);
    atomic_store(&g_fault_target, nth);
}
uint64_t duckvep_fault_allocations(void) { return atomic_load(&g_fault_count); }
uint64_t duckvep_fault_fired(void) { return atomic_load(&g_fault_fired); }
#endif

static void *
allocate(duckvep_budget_owner_t owner, size_t size, int zero)
{
    size_t charged;
    header_t *block;

#ifdef DUCKVEP_FAULT_INJECTION
    if (fault_fires())
        return NULL;
#endif
    charged = duckvep_budget_charge_for(size);
    if (charged == SIZE_MAX || !charge(owner, charged))
        return NULL;
    block = zero ? calloc(1, size + HEADER_BYTES) : malloc(size + HEADER_BYTES);
    if (block == NULL) {
        uncharge(owner, charged);
        return NULL;
    }
    block->charged = charged;
    block->owner = (uint32_t)owner;
    block->magic = HEADER_MAGIC;
    return (char *)block + HEADER_BYTES;
}

void *
duckvep_budget_malloc(duckvep_budget_owner_t owner, size_t size)
{
    return allocate(owner, size, 0);
}

void *
duckvep_budget_calloc(duckvep_budget_owner_t owner, size_t count, size_t width)
{
    if (width != 0 && count > SIZE_MAX / width)
        return NULL;
    return allocate(owner, count * width, 1);
}

static header_t *
header_of(const void *pointer)
{
    header_t *block = (header_t *)((char *)pointer - HEADER_BYTES);

    if (block->magic != HEADER_MAGIC) {
        (void)fprintf(stderr, "duckvep: budget free/realloc of a foreign pointer\n");
        abort();
    }
    return block;
}

void *
duckvep_budget_realloc(duckvep_budget_owner_t owner, void *pointer, size_t size)
{
    header_t *old, *grown;
    size_t charged;
    duckvep_budget_owner_t old_owner;
    uint64_t old_charged;

    if (pointer == NULL)
        return allocate(owner, size, 0);
    old = header_of(pointer);
    old_owner = (duckvep_budget_owner_t)old->owner;
    old_charged = old->charged;
#ifdef DUCKVEP_FAULT_INJECTION
    if (fault_fires())
        return NULL;
#endif
    charged = duckvep_budget_charge_for(size);
    /* The new block is charged in full while the old one is still charged:
     * that is the realloc overlap.  Shrinking pays for nothing extra. */
    if (charged > old_charged) {
        if (!charge(old_owner, charged - old_charged))
            return NULL;
    }
    grown = realloc(old, size + HEADER_BYTES);
    if (grown == NULL) {
        if (charged > old_charged)
            uncharge(old_owner, charged - old_charged);
        return NULL;
    }
    if (charged < old_charged)
        uncharge(old_owner, old_charged - charged);
    grown->charged = charged;
    return (char *)grown + HEADER_BYTES;
}

void
duckvep_budget_free(void *pointer)
{
    header_t *block;
    uint64_t charged;
    duckvep_budget_owner_t owner;

    if (pointer == NULL)
        return;
    block = header_of(pointer);
    charged = block->charged;
    owner = (duckvep_budget_owner_t)block->owner;
    block->magic = 0;
    free(block);
    uncharge(owner, charged);
}

size_t
duckvep_budget_charged(const void *pointer)
{
    return pointer != NULL ? (size_t)header_of(pointer)->charged : 0;
}

char *
duckvep_budget_strndup(duckvep_budget_owner_t owner, const char *text,
    size_t length)
{
    char *copy = duckvep_budget_malloc(owner, length + 1);

    if (copy != NULL) {
        if (length != 0)
            memcpy(copy, text, length);
        copy[length] = '\0';
    }
    return copy;
}

char *
duckvep_budget_strdup(duckvep_budget_owner_t owner, const char *text)
{
    return duckvep_budget_strndup(owner, text, strlen(text));
}

int
duckvep_budget_reserve(duckvep_budget_owner_t owner, uint64_t bytes)
{
#ifdef DUCKVEP_FAULT_INJECTION
    if (fault_fires())
        return 0;
#endif
    return charge(owner, bytes);
}

void
duckvep_budget_unreserve(duckvep_budget_owner_t owner, uint64_t bytes)
{
    uncharge(owner, bytes);
}

int
duckvep_budget_set_limit(uint64_t bytes)
{
    if (bytes == 0 || bytes < atomic_load(&g_total))
        return 0;
    atomic_store(&g_limit, bytes);
    return 1;
}

uint64_t
duckvep_budget_limit(void)
{
    return atomic_load(&g_limit);
}

void
duckvep_budget_stats(duckvep_budget_stats_t *stats)
{
    int owner;

    memset(stats, 0, sizeof *stats);
    stats->limit = atomic_load(&g_limit);
    for (owner = 0; owner < DUCKVEP_OWNER_COUNT; owner++) {
        stats->current[owner] = atomic_load(&g_current[owner]);
        stats->high_water[owner] = atomic_load(&g_high[owner]);
    }
    stats->total_current = atomic_load(&g_total);
    stats->total_high_water = atomic_load(&g_total_high);
    stats->refusals = atomic_load(&g_refusals);
    stats->charges = atomic_load(&g_charges);
}

void
duckvep_budget_reset_high_water(void)
{
    int owner;

    for (owner = 0; owner < DUCKVEP_OWNER_COUNT; owner++)
        atomic_store(&g_high[owner], atomic_load(&g_current[owner]));
    atomic_store(&g_total_high, atomic_load(&g_total));
}

void
duckvep_budget_worker_limits(duckvep_budget_worker_limits_t *limits)
{
    limits->max_workers = atomic_load(&g_workers);
    limits->scratch_bytes = atomic_load(&g_scratch);
    limits->emit_bytes = atomic_load(&g_emit);
    limits->idle_bytes = atomic_load(&g_idle);
}

int
duckvep_budget_set_worker_limits(const duckvep_budget_worker_limits_t *limits)
{
    if (limits->max_workers == 0)
        return 0;
    atomic_store(&g_workers, limits->max_workers);
    atomic_store(&g_scratch, limits->scratch_bytes);
    atomic_store(&g_emit, limits->emit_bytes);
    atomic_store(&g_idle, limits->idle_bytes);
    return 1;
}
