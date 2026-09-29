/*
 * duckvep_budget.h - process-wide native allocation budget.
 *
 * Every first-party native owner (model arrays and indexes, reference caches,
 * per-worker workspaces, result and text arenas) allocates through these
 * functions.  Each block carries a 16-byte header holding the charged size, so
 * the budget is charged before the system allocator is called, uses
 * page-rounded capacity for large blocks, and charges realloc while the old
 * block is still live (the overlap).  A refusal returns NULL and records a
 * thread-local capacity report that duckvep_budget_take_failure() hands to the
 * error path, so callers can raise an explicit capacity error instead of a
 * generic out-of-memory message.
 *
 * Pointers obtained here must be released with duckvep_budget_free() only.
 */
#ifndef DUCKVEP_BUDGET_H
#define DUCKVEP_BUDGET_H

#include <stddef.h>
#include <stdint.h>

#define DUCKVEP_BUDGET_DEFAULT_BYTES (UINT64_C(4) << 30)
#define DUCKVEP_BUDGET_DEFAULT_WORKERS 6u
#define DUCKVEP_BUDGET_DEFAULT_SCRATCH_BYTES (UINT64_C(128) << 20)
#define DUCKVEP_BUDGET_DEFAULT_EMIT_BYTES (UINT64_C(256) << 20)
#define DUCKVEP_BUDGET_DEFAULT_IDLE_BYTES (UINT64_C(64) << 20)

typedef enum duckvep_budget_owner {
    DUCKVEP_OWNER_MODEL = 0,     /* model arrays, sequence pools, kernel caches */
    DUCKVEP_OWNER_INDEX,         /* cgranges interval indexes */
    DUCKVEP_OWNER_REFERENCE,     /* reference readers: byte caches, faidx reservation */
    DUCKVEP_OWNER_WORKSPACE,     /* kernel workspaces, cursors, event scratch */
    DUCKVEP_OWNER_SCRATCH,       /* per-worker scratch lease */
    DUCKVEP_OWNER_EMIT,          /* per-worker result/text arenas */
    DUCKVEP_OWNER_CONTROL,       /* registry, bind/init state, SQL text, options */
    DUCKVEP_OWNER_COUNT
} duckvep_budget_owner_t;

typedef struct duckvep_budget_stats {
    uint64_t limit;
    uint64_t current[DUCKVEP_OWNER_COUNT];
    uint64_t high_water[DUCKVEP_OWNER_COUNT];
    uint64_t total_current;
    uint64_t total_high_water;
    uint64_t refusals;
    uint64_t charges;
} duckvep_budget_stats_t;

const char *duckvep_budget_owner_name(duckvep_budget_owner_t owner);

void *duckvep_budget_malloc(duckvep_budget_owner_t owner, size_t size);
void *duckvep_budget_calloc(duckvep_budget_owner_t owner, size_t count, size_t width);
/* Growing charges the full new block before the old one is released. */
void *duckvep_budget_realloc(duckvep_budget_owner_t owner, void *pointer, size_t size);
void duckvep_budget_free(void *pointer);
/* Bytes charged for a live block (0 for NULL). */
size_t duckvep_budget_charged(const void *pointer);
/* Bytes that would be charged for a block of this size. */
size_t duckvep_budget_charge_for(size_t size);
char *duckvep_budget_strdup(duckvep_budget_owner_t owner, const char *text);
char *duckvep_budget_strndup(duckvep_budget_owner_t owner, const char *text, size_t length);

/* Fixed reservations for memory that first-party code cannot route through
 * duckvep_budget_malloc (third-party transport buffers). Returns 0 if over budget. */
int duckvep_budget_reserve(duckvep_budget_owner_t owner, uint64_t bytes);
void duckvep_budget_unreserve(duckvep_budget_owner_t owner, uint64_t bytes);

/* Limit control.  A limit below the current use is refused (returns 0). */
int duckvep_budget_set_limit(uint64_t bytes);
uint64_t duckvep_budget_limit(void);
void duckvep_budget_stats(duckvep_budget_stats_t *stats);
void duckvep_budget_reset_high_water(void);

/* Bounded-execution parameters (per-worker leases and admission). */
typedef struct duckvep_budget_worker_limits {
    uint32_t max_workers;
    uint64_t scratch_bytes;
    uint64_t emit_bytes;
    uint64_t idle_bytes;
} duckvep_budget_worker_limits_t;
void duckvep_budget_worker_limits(duckvep_budget_worker_limits_t *limits);
int duckvep_budget_set_worker_limits(const duckvep_budget_worker_limits_t *limits);

/* Capacity reports.  The budget (or a lease) records why the last allocation on
 * this thread was refused; the error path consumes it exactly once. */
void duckvep_budget_note_failure(const char *what, uint64_t requested,
    uint64_t in_use, uint64_t limit);
/* Writes "capacity error: ..." into buffer and clears the record. Returns 0
 * (buffer untouched) when the last failure was not a capacity refusal. */
int duckvep_budget_take_failure(char *buffer, size_t size);
void duckvep_budget_clear_failure(void);

#ifdef DUCKVEP_FAULT_INJECTION
/* Test builds only: fail the Nth budget allocation (1-based) after arming.
 * n == 0 disarms.  Returns nothing; duckvep_fault_allocations() counts every
 * budget allocation attempt since the last arm. */
void duckvep_fault_arm(uint64_t nth);
uint64_t duckvep_fault_allocations(void);
uint64_t duckvep_fault_fired(void);
#endif

#endif
