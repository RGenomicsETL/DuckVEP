/* Alignment score/tie and difference-run semantics from Ensembl Variation
 * release/116 (Apache-2.0), credited in duckvep_sequence_diff.h. This bounded
 * implementation stores two score rows and an exact band of traceback cells. */
#include "duckvep_sequence_diff.h"

#include <string.h>

enum { DIAGONAL, DELETE_REFERENCE, INSERT_ALTERNATE, MATCH };

static size_t minimum(size_t a, size_t b) { return a < b ? a : b; }
static size_t maximum(size_t a, size_t b) { return a > b ? a : b; }

/* Every byte must be nonzero ASCII and not the gap character '-'. Eight bytes at a time: a zero byte, a byte with
 * the high bit, or a '-' (0x2d) in any lane rejects the word. */
static int sequence_bytes_valid(const uint8_t *bytes, size_t length) {
    const uint64_t ones = UINT64_C(0x0101010101010101), highs = UINT64_C(0x8080808080808080),
        dashes = UINT64_C(0x2d2d2d2d2d2d2d2d);
    size_t i = 0u;
    for (; i + 8u <= length; i += 8u) {
        uint64_t word, dash;
        memcpy(&word, bytes + i, sizeof word);
        dash = word ^ dashes;
        if (((word - ones) & ~word & highs) || (word & highs) || ((dash - ones) & ~dash & highs)) return 0;
    }
    for (; i < length; i++)
        if (!bytes[i] || bytes[i] >= 128u || bytes[i] == '-') return 0;
    return 1;
}

/* Positional comparison (no alignment): the shorter sequence is padded with gaps, so the columns are the shared
 * positions followed by one gap run. Differences are the maximal runs of mismatching columns and the tail. When
 * out is NULL only the run count is produced. Equal to the traceback of the implicit diagonal path. */
static void positional_differences(const uint8_t *ref, size_t n, const uint8_t *alt, size_t m,
    duckvep_sequence_difference_t *out, duckvep_sequence_diff_result_t *result) {
    size_t shared = minimum(n, m), columns = maximum(n, m), count = 0u, k = 0u;
    while (k < shared) {
        while (k + 8u <= shared) {
            uint64_t a, b;
            memcpy(&a, ref + k, sizeof a); memcpy(&b, alt + k, sizeof b);
            if (a != b) break;
            k += 8u;
        }
        while (k < shared && ref[k] == alt[k]) k++;
        if (k == shared) break;
        size_t first = k;
        while (k < shared && ref[k] != alt[k]) k++;
        if (out) out[count] = (duckvep_sequence_difference_t){first, first, k - first, k - first, first};
        count++;
    }
    if (n != m) {
        if (out) out[count] = (duckvep_sequence_difference_t){shared, shared, n - shared, m - shared, shared};
        count++;
    }
    result->count = count; result->alignment_length = columns;
}

/* This is a feasible alignment, not a trimmed problem: retain the common
 * prefix/suffix and align the middle positionally. NW still sees every byte,
 * so repeat-associated tie placement is unaffected. */
static uint64_t cost_bound(const uint8_t *ref, size_t n, const uint8_t *alt, size_t m) {
    size_t first = 0u, last = 0u, shared = minimum(n, m);
    while (first < shared && ref[first] == alt[first]) first++;
    while (last < shared - first && ref[n - 1u - last] == alt[m - 1u - last]) last++;
    uint64_t cost = (uint64_t)(maximum(n, m) - shared) * 3u;
    for (size_t i = first; i < shared - last; i++) if (ref[i] != alt[i]) cost += 4u;
    return cost;
}

/* Fills the banded score rows and traceback directions and returns the optimal cost within the band.
 *
 * The common prefix of length p is aligned on the main diagonal at zero cost. For every cell with max(i,j) <= p,
 * D(i,j) = 3|j-i|: gaps alone are a lower bound and diagonal-then-gap reaches it. The tie rules then always give
 * DIAGONAL on the diagonal, INSERT_ALTERNATE above it and DELETE_REFERENCE below it, in every cell whose whole band
 * neighborhood lies inside the prefix. Rows up to prefix_rows = p - band are therefore closed-form: they are neither
 * computed nor stored, and traceback_direction() regenerates their directions. Rows after them are computed from
 * the closed-form row, so the stored directions and the cost equal those of a full-prefix computation. */
static uint64_t build_trace(const uint8_t *ref, size_t n, const uint8_t *alt, size_t m,
                            size_t band, size_t stride, size_t prefix_rows,
                            const duckvep_sequence_diff_scratch_t *b) {
    const uint64_t infinity = UINT64_MAX / 2u;
    uint64_t *previous = b->scores, *current = b->scores + m + 1u;
    for (size_t j = 0u; j <= m; j++) previous[j] = infinity;
    if (prefix_rows == 0u) {
        for (size_t j = 0u; j <= minimum(m, band); j++) {
            previous[j] = (uint64_t)j * 3u; b->trace[j] = INSERT_ALTERNATE;
        }
    } else {
        size_t low = prefix_rows > band ? prefix_rows - band : 0u;
        for (size_t j = low; j <= minimum(m, prefix_rows + band); j++)
            previous[j] = (uint64_t)(j > prefix_rows ? j - prefix_rows : prefix_rows - j) * 3u;
    }
    for (size_t i = prefix_rows + 1u; i <= n; i++) {
        size_t start = i > band ? i - band : 0u;
        size_t end = minimum(m, i + band);
        size_t row = (i - prefix_rows - (prefix_rows ? 1u : 0u)) * stride;
        if (start) current[start - 1u] = infinity;
        if (end > minimum(m, i - 1u + band)) previous[end] = infinity;
        /* Keep the three neighbours in registers: previous[j-1], previous[j] and the cell to the left. */
        uint8_t *cells = b->trace + row;
        uint64_t up_left = start ? previous[start - 1u] : infinity, left = start ? infinity : 0u;
        for (size_t j = start; j <= end; j++) {
            uint8_t direction = DELETE_REFERENCE;
            uint64_t up = previous[j];
            uint64_t score = (uint64_t)i * 3u;
            if (j) {
                uint64_t diagonal = up_left + (ref[i - 1u] == alt[j - 1u] ? 0u : 4u);
                uint64_t insertion = left + 3u, deletion = up + 3u;
                /* Minimizing 4*substitutions + 3*gaps is exactly maximizing
                 * matches - substitutions - gaps at fixed endpoint lengths. */
                if (diagonal < insertion && diagonal < deletion) {
                    score = diagonal; direction = DIAGONAL;
                } else if (insertion < deletion) {
                    score = insertion; direction = INSERT_ALTERNATE;
                } else score = deletion;
            }
            current[j] = score;
            cells[j - start] = direction;
            up_left = up; left = score;
        }
        uint64_t *swap = previous; previous = current; current = swap;
    }
    return previous[m];
}

/* The stored direction of cell (i,j); rows inside the closed-form prefix are regenerated (see build_trace). */
static uint8_t traceback_direction(const uint8_t *trace, size_t band, size_t stride,
                                   size_t prefix_rows, size_t i, size_t j) {
    if (prefix_rows && i <= prefix_rows) return j == i ? DIAGONAL : j > i ? INSERT_ALTERNATE : DELETE_REFERENCE;
    size_t start = i > band ? i - band : 0u;
    size_t row = prefix_rows ? i - prefix_rows - 1u : i;
    return trace[row * stride + j - start];
}

/* Walk once to count complete runs and once to fill. Reverse traversal groups
 * exactly the same adjacent gap classes as upstream's forward difference scan.
 * Reaching the diagonal inside the closed-form prefix leaves only matching columns. */
static void traceback(const uint8_t *ref, size_t n, const uint8_t *alt, size_t m,
    const uint8_t *trace, size_t band, size_t stride, size_t prefix_rows,
    duckvep_sequence_difference_t *out, duckvep_sequence_diff_result_t *result) {
    size_t i = n, j = m, columns = 0u, count = 0u;
    uint8_t previous = MATCH;
    while (i || j) {
        uint8_t direction;
        if (trace) {
            if (prefix_rows && i <= prefix_rows && i == j) { columns += i; break; }
            direction = traceback_direction(trace, band, stride, prefix_rows, i, j);
        } else {
            direction = i == j ? DIAGONAL : i > j ? DELETE_REFERENCE : INSERT_ALTERNATE;
        }
        size_t rn = direction != INSERT_ALTERNATE, an = direction != DELETE_REFERENCE;
        i -= rn; j -= an;
        uint8_t kind = direction == DIAGONAL && ref[i] == alt[j] ? MATCH : direction;
        if (kind != MATCH) {
            if (kind != previous) {
                count++;
                if (out) out[count - 1u] = (duckvep_sequence_difference_t){i, j, 0u, 0u, columns};
            }
            if (out) {
                duckvep_sequence_difference_t *d = &out[count - 1u];
                d->ref_start0 = i; d->alt_start0 = j;
                d->ref_length += rn; d->alt_length += an;
            }
        }
        previous = kind; columns++;
    }
    if (out) {
        for (size_t k = 0u; k < count; k++)
            out[k].alignment_start0 = columns - out[k].alignment_start0 -
                maximum(out[k].ref_length, out[k].alt_length);
        for (size_t k = 0u; k < count / 2u; k++) {
            duckvep_sequence_difference_t swap = out[k];
            out[k] = out[count - 1u - k]; out[count - 1u - k] = swap;
        }
    }
    result->count = count; result->alignment_length = columns;
}

duckvep_sequence_diff_status_t duckvep_sequence_differences(
    const uint8_t *ref, size_t n, const uint8_t *alt, size_t m,
    int align_indels, const duckvep_sequence_diff_scratch_t *b,
    duckvep_sequence_difference_t *out, size_t capacity, duckvep_sequence_diff_result_t *result) {
    if (!result) return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
    memset(result, 0, sizeof(*result));
    if ((!ref && n) || (!alt && m) || (!out && capacity) ||
        (align_indels != 0 && align_indels != 1) ||
        n > SIZE_MAX / 8u || m > SIZE_MAX / 8u)
        return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
    if (!sequence_bytes_valid(ref, n) || !sequence_bytes_valid(alt, m)) return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
    size_t band = 0u, stride = 0u, prefix_rows = 0u;
    const uint8_t *trace = NULL;
    if (align_indels && n && m && (n != m || memcmp(ref, alt, n))) {
        /* Every optimum costs at most this feasible path. Leaving |i-j| <= U/3
         * needs more than U in gap cost alone, so no optimal traceback is lost. */
        band = (size_t)minimum(maximum(n, m), (size_t)(cost_bound(ref, n, alt, m) / 3u));
        /* That band is a bound, not a need: an optimum that stays within |i-j| <= w is exact whenever every path
         * leaving that band costs more. Reaching offset w+1 and returning to the final offset costs at least
         * 3*(2(w+1)-|m-n|) in gaps, so a banded optimum below that value has every optimal path inside the band,
         * every traceback direction on those paths is unchanged, and the traceback is identical to the full band.
         * Start at the length change and widen (never beyond the bound). Capacity is checked against the band of
         * each attempt, (n+1) x stride cells, so a long sequence with a short indel no longer needs the trace of its
         * feasible bound; when the band that is still to be tried does not fit, that band's cell count is reported. */
        size_t gap = n > m ? n - m : m - n, width = minimum(band, gap + 3u);
        size_t prefix = 0u, shared = minimum(n, m);
        while (prefix < shared && ref[prefix] == alt[prefix]) prefix++;
        for (;;) {
            size_t narrow = width >= m / 2u ? m + 1u : width * 2u + 1u;
            if (n + 1u > SIZE_MAX / narrow) return DUCKVEP_SEQUENCE_DIFF_INVALID_ARG;
            result->trace_cells = (n + 1u) * narrow;
            if (!b || !b->trace || result->trace_cells > b->trace_capacity)
                return DUCKVEP_SEQUENCE_DIFF_TRACE_FULL;
            if (!b->scores || m + 1u > b->score_capacity / 2u)
                return DUCKVEP_SEQUENCE_DIFF_SCORE_FULL;
            size_t rows = prefix > width ? prefix - width : 0u;
            uint64_t cost = build_trace(ref, n, alt, m, width, narrow, rows, b);
            if (width >= band || cost < (uint64_t)3u * (2u * (width + 1u) - gap)) {
                band = width; stride = narrow; prefix_rows = rows; break;
            }
            width = minimum(band, width * 2u + 1u);
        }
        trace = b->trace;
    }
    if (!trace) {
        positional_differences(ref, n, alt, m, NULL, result);
        if (result->count > capacity) return DUCKVEP_SEQUENCE_DIFF_OUTPUT_FULL;
        positional_differences(ref, n, alt, m, out, result);
        return DUCKVEP_SEQUENCE_DIFF_OK;
    }
    /* Every difference occupies at least one alignment column and an alignment has at most n+m columns, so with
     * that much capacity the counting pass cannot report OUTPUT_FULL and is skipped. */
    if (capacity < n + m) {
        traceback(ref, n, alt, m, trace, band, stride, prefix_rows, NULL, result);
        if (result->count > capacity) return DUCKVEP_SEQUENCE_DIFF_OUTPUT_FULL;
    }
    traceback(ref, n, alt, m, trace, band, stride, prefix_rows, out, result);
    return DUCKVEP_SEQUENCE_DIFF_OK;
}
