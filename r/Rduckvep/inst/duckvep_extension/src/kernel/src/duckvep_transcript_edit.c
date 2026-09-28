/*
 * duckvep_transcript_edit.c — shared transcript edit projection.
 */
#include "duckvep_transcript_edit.h"

#include "duckvep_projection.h"

#include <limits.h>
#include <string.h>

static uint64_t transcript_projection_abs_diff(uint32_t left, uint32_t right) {
    return left >= right ? (uint64_t)left - (uint64_t)right
                         : (uint64_t)right - (uint64_t)left;
}

static void transcript_projection_order_pair(
    uint32_t first, uint32_t last, uint32_t *lower, uint32_t *upper) {

    if (first != 0u && last != 0u && first > last) {
        *lower = last;
        *upper = first;
    } else {
        *lower = first;
        *upper = last;
    }
}

static int transcript_edit_model_slice_ok(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    size_t                           *offset_out,
    size_t                           *count_out) {

    size_t offset;
    size_t count;

    if (transcripts == NULL || exons == NULL || offset_out == NULL ||
        count_out == NULL || tx_idx >= transcripts->transcript_count ||
        transcripts->exon_offset == NULL || transcripts->exon_count == NULL ||
        transcripts->start1 == NULL || transcripts->end1 == NULL ||
        transcripts->strand == NULL || exons->start1 == NULL ||
        exons->end1 == NULL || exons->cdna_start1 == NULL ||
        exons->cdna_end1 == NULL) {
        return 0;
    }
    offset = (size_t)transcripts->exon_offset[tx_idx];
    count = (size_t)transcripts->exon_count[tx_idx];
    if (count == 0u || offset > exons->exon_count ||
        count > exons->exon_count - offset ||
        transcripts->start1[tx_idx] == 0u ||
        transcripts->end1[tx_idx] < transcripts->start1[tx_idx] ||
        (transcripts->strand[tx_idx] != 1 &&
         transcripts->strand[tx_idx] != -1)) {
        return 0;
    }
    *offset_out = offset;
    *count_out = count;
    return 1;
}

static int transcript_edit_project_exon_coordinate(
    const duckvep_exon_model_t          *exons,
    uint32_t                             exon_hint,
    uint32_t                             genomic_pos1,
    int8_t                               strand,
    duckvep_transcript_coordinate_t     *out) {

    uint64_t cdna;

    if (exons == NULL || out == NULL || exon_hint >= exons->exon_count ||
        genomic_pos1 < exons->start1[exon_hint] ||
        genomic_pos1 > exons->end1[exon_hint] ||
        exons->cdna_start1[exon_hint] == 0u ||
        exons->cdna_end1[exon_hint] < exons->cdna_start1[exon_hint]) {
        return 0;
    }
    cdna = strand > 0
        ? (uint64_t)exons->cdna_start1[exon_hint] +
          (uint64_t)(genomic_pos1 - exons->start1[exon_hint])
        : (uint64_t)exons->cdna_start1[exon_hint] +
          (uint64_t)(exons->end1[exon_hint] - genomic_pos1);
    if (cdna > (uint64_t)exons->cdna_end1[exon_hint] || cdna > UINT32_MAX) {
        return 0;
    }
    memset(out, 0, sizeof *out);
    out->genomic_pos1 = genomic_pos1;
    out->cdna_anchor1 = (uint32_t)cdna;
    out->exon_idx = exon_hint;
    out->exonic = 1u;
    return 1;
}

static duckvep_transcript_edit_status_t transcript_edit_project_pair(
    const duckvep_transcript_model_t    *transcripts,
    const duckvep_exon_model_t          *exons,
    size_t                               tx_idx,
    uint32_t                             genomic_first1,
    uint32_t                             genomic_last1,
    int8_t                               strand,
    uint32_t                             exon_hint,
    duckvep_transcript_coordinate_t     *first_out,
    duckvep_transcript_coordinate_t     *last_out) {

    duckvep_transcript_coordinate_t genomic_first;
    duckvep_transcript_coordinate_t genomic_last;
    duckvep_transcript_edit_status_t status;

    if (first_out == NULL || last_out == NULL || genomic_first1 == 0u ||
        genomic_last1 == 0u) {
        return DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT;
    }
    if (exon_hint != UINT32_MAX &&
        transcript_edit_project_exon_coordinate(
            exons, exon_hint, genomic_first1, strand, &genomic_first)) {
        status = DUCKVEP_TRANSCRIPT_EDIT_OK;
    } else {
        status = duckvep_project_transcript_coordinate(
            transcripts, exons, tx_idx, genomic_first1, &genomic_first);
    }
    if (status != DUCKVEP_TRANSCRIPT_EDIT_OK) return status;
    if (genomic_first1 == genomic_last1) {
        genomic_last = genomic_first;
    } else if (exon_hint != UINT32_MAX &&
               transcript_edit_project_exon_coordinate(
                   exons, exon_hint, genomic_last1, strand, &genomic_last)) {
        status = DUCKVEP_TRANSCRIPT_EDIT_OK;
    } else {
        status = duckvep_project_transcript_coordinate(
            transcripts, exons, tx_idx, genomic_last1, &genomic_last);
        if (status != DUCKVEP_TRANSCRIPT_EDIT_OK) return status;
    }
    if (strand > 0) {
        *first_out = genomic_first;
        *last_out = genomic_last;
    } else {
        *first_out = genomic_last;
        *last_out = genomic_first;
    }
    return DUCKVEP_TRANSCRIPT_EDIT_OK;
}

duckvep_transcript_edit_status_t duckvep_project_transcript_coordinate(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          genomic_pos1,
    duckvep_transcript_coordinate_t  *out) {

    duckvep_transcript_coordinate_t result;
    size_t offset;
    size_t count;
    size_t i;
    size_t lower_idx = SIZE_MAX;
    size_t upper_idx = SIZE_MAX;
    uint32_t lower_end = 0u;
    uint32_t upper_start = UINT32_MAX;
    uint32_t lower_distance;
    uint32_t upper_distance;
    uint32_t cdna;
    uint32_t exon_idx;
    int8_t strand;

    if (out == NULL) return DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (!transcript_edit_model_slice_ok(transcripts, exons, tx_idx,
                                        &offset, &count) ||
        genomic_pos1 == 0u) {
        return DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG;
    }
    if (genomic_pos1 < transcripts->start1[tx_idx] ||
        genomic_pos1 > transcripts->end1[tx_idx]) {
        return DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT;
    }
    memset(&result, 0, sizeof result);
    result.genomic_pos1 = genomic_pos1;
    if (duckvep_project_genomic_to_cdna(transcripts, exons, tx_idx,
                                        genomic_pos1, &cdna, &exon_idx)) {
        result.cdna_anchor1 = cdna;
        result.exon_idx = exon_idx;
        result.exonic = 1u;
        *out = result;
        return DUCKVEP_TRANSCRIPT_EDIT_OK;
    }

    for (i = offset; i < offset + count; i++) {
        if (exons->end1[i] < genomic_pos1 &&
            (lower_idx == SIZE_MAX || exons->end1[i] > lower_end)) {
            lower_idx = i;
            lower_end = exons->end1[i];
        }
        if (exons->start1[i] > genomic_pos1 &&
            (upper_idx == SIZE_MAX || exons->start1[i] < upper_start)) {
            upper_idx = i;
            upper_start = exons->start1[i];
        }
    }
    if (lower_idx == SIZE_MAX || upper_idx == SIZE_MAX) {
        return DUCKVEP_TRANSCRIPT_EDIT_INVALID_PROJECTION;
    }
    lower_distance = genomic_pos1 - lower_end;
    upper_distance = upper_start - genomic_pos1;
    if (lower_distance == 0u || upper_distance == 0u ||
        lower_distance > INT32_MAX || upper_distance > INT32_MAX) {
        return DUCKVEP_TRANSCRIPT_EDIT_OUT_OF_RANGE;
    }
    strand = transcripts->strand[tx_idx];
    if (lower_distance < upper_distance ||
        (lower_distance == upper_distance && strand > 0)) {
        result.exon_idx = (uint32_t)lower_idx;
        if (strand > 0) {
            result.cdna_anchor1 = exons->cdna_end1[lower_idx];
            result.intron_offset = (int32_t)lower_distance;
        } else {
            result.cdna_anchor1 = exons->cdna_start1[lower_idx];
            result.intron_offset = -(int32_t)lower_distance;
        }
    } else {
        result.exon_idx = (uint32_t)upper_idx;
        if (strand > 0) {
            result.cdna_anchor1 = exons->cdna_start1[upper_idx];
            result.intron_offset = -(int32_t)upper_distance;
        } else {
            result.cdna_anchor1 = exons->cdna_end1[upper_idx];
            result.intron_offset = (int32_t)upper_distance;
        }
    }
    if (result.cdna_anchor1 == 0u) {
        return DUCKVEP_TRANSCRIPT_EDIT_INVALID_PROJECTION;
    }
    *out = result;
    return DUCKVEP_TRANSCRIPT_EDIT_OK;
}

int duckvep_transcript_projection_facts_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_variant_batch_t    *variants,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                          region_mask,
    const duckvep_coding_context_t   *coding_context,
    duckvep_transcript_projection_facts_t *out) {

    const uint8_t *feature_ref;
    const uint8_t *feature_alt;
    uint16_t feature_ref_length;
    uint16_t feature_alt_length;
    uint32_t feature_start1;
    uint32_t feature_end1;
    uint32_t first_cdna = 0u;
    uint32_t last_cdna = 0u;
    uint32_t raw_cdna_start = 0u;
    uint32_t raw_cdna_end = 0u;
    uint32_t raw_cds_start = 0u;
    uint32_t raw_cds_end = 0u;
    uint32_t coding_start_cdna;
    uint32_t coding_end_cdna;
    uint8_t first_exon_phase;
    size_t exon_offset;
    size_t exon_count;
    size_t i;
    int insertion;
    int within_transcript;

    if (out == NULL) return 0;
    memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || variants == NULL ||
        event == NULL || tx_idx >= transcripts->transcript_count ||
        (size_t)variant_idx >= variants->count ||
        transcripts->flags == NULL ||
        !transcript_edit_model_slice_ok(transcripts, exons, tx_idx,
                                        &exon_offset, &exon_count) ||
        !duckvep_event_feature_alleles(
            variants, (size_t)variant_idx, event,
            &feature_ref, &feature_ref_length,
            &feature_alt, &feature_alt_length)) {
        return 0;
    }

    feature_start1 = event->feature_start1;
    feature_end1 = event->feature_end1;
    if (feature_start1 == 0u || feature_end1 == 0u) return 0;
    insertion = feature_start1 > feature_end1;
    within_transcript = feature_end1 >= transcripts->start1[tx_idx] &&
                        feature_start1 <= transcripts->end1[tx_idx];

    out->output_allele = feature_alt;
    out->output_allele_length = (size_t)feature_alt_length;
    out->feature_ref_length = feature_ref_length;
    out->feature_alt_length = feature_alt_length;
    out->interbase = (uint8_t)insertion;
    out->exon_total = (uint32_t)exon_count;
    out->intron_total = exon_count > 0u ? (uint32_t)(exon_count - 1u) : 0u;
    out->cds_start_nf = (uint8_t)((transcripts->flags[tx_idx] &
        (uint64_t)DUCKVEP_TX_CDS_START_NF) != 0u);
    out->cds_end_nf = (uint8_t)((transcripts->flags[tx_idx] &
        (uint64_t)DUCKVEP_TX_CDS_END_NF) != 0u);

    if (transcripts->strand[tx_idx] > 0) {
        (void)duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, feature_start1, &first_cdna, NULL);
        (void)duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, feature_end1, &last_cdna, NULL);
    } else {
        (void)duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, feature_end1, &first_cdna, NULL);
        (void)duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, feature_start1, &last_cdna, NULL);
    }
    if (within_transcript) {
        if (insertion) {
            raw_cdna_start = first_cdna;
            if (raw_cdna_start == 0u && last_cdna != 0u &&
                last_cdna != UINT32_MAX) {
                raw_cdna_start = last_cdna + 1u;
            }
            raw_cdna_end = last_cdna;
            if (raw_cdna_end == 0u && first_cdna > 1u) {
                raw_cdna_end = first_cdna - 1u;
            }
        } else {
            raw_cdna_start = first_cdna;
            raw_cdna_end = last_cdna;
        }
    }
    transcript_projection_order_pair(
        raw_cdna_start, raw_cdna_end, &out->cdna_start, &out->cdna_end);

    first_exon_phase = exons->phase != NULL &&
        exons->phase[exon_offset] > 0 ?
        (uint8_t)exons->phase[exon_offset] : 0u;
    if (duckvep_project_coding_cdna_bounds(
            transcripts, exons, tx_idx, &coding_start_cdna,
            &coding_end_cdna, NULL, NULL)) {
        int insertion_range_valid = !insertion ||
            (raw_cdna_start != 0u && raw_cdna_end != 0u &&
             raw_cdna_start <= coding_end_cdna &&
             raw_cdna_end >= coding_start_cdna);
        uint64_t projected;

        if (raw_cdna_start >= coding_start_cdna &&
            raw_cdna_start <= coding_end_cdna && insertion_range_valid) {
            projected = (uint64_t)raw_cdna_start - coding_start_cdna +
                first_exon_phase + 1u;
            if (projected <= UINT32_MAX) raw_cds_start = (uint32_t)projected;
        }
        if (raw_cdna_end >= coding_start_cdna &&
            raw_cdna_end <= coding_end_cdna && insertion_range_valid) {
            projected = (uint64_t)raw_cdna_end - coding_start_cdna +
                first_exon_phase + 1u;
            if (projected <= UINT32_MAX) raw_cds_end = (uint32_t)projected;
        }
    }
    transcript_projection_order_pair(
        raw_cds_start, raw_cds_end, &out->cds_start, &out->cds_end);
    if (raw_cds_start != 0u) {
        out->protein_start = (raw_cds_start - 1u) / 3u + 1u;
    }
    if (raw_cds_end != 0u) {
        out->protein_end = (raw_cds_end - 1u) / 3u + 1u;
    }
    transcript_projection_order_pair(
        out->protein_start, out->protein_end,
        &out->protein_start, &out->protein_end);

    for (i = 0u; i < exon_count; i++) {
        size_t exon_idx = exon_offset + i;
        uint32_t rank = (uint32_t)i + 1u;
        if (feature_start1 <= exons->end1[exon_idx] &&
            feature_end1 >= exons->start1[exon_idx]) {
            if (out->exon_first == 0u) out->exon_first = rank;
            out->exon_last = rank;
        }
        if (i + 1u < exon_count) {
            size_t next_idx = exon_idx + 1u;
            uint32_t high_start = exons->start1[exon_idx] >
                exons->start1[next_idx] ? exons->start1[exon_idx]
                                         : exons->start1[next_idx];
            uint32_t low_end = exons->end1[exon_idx] <
                exons->end1[next_idx] ? exons->end1[exon_idx]
                                       : exons->end1[next_idx];
            uint64_t intron_high = high_start == 0u ? 0u :
                (uint64_t)high_start - 1u;
            uint64_t intron_low = (uint64_t)low_end + 1u;
            if ((uint64_t)feature_start1 <= intron_high &&
                (uint64_t)feature_end1 >= intron_low) {
                if (out->intron_first == 0u) out->intron_first = rank;
                out->intron_last = rank;
            }
        }
    }

    if ((region_mask & ((uint32_t)DUCKVEP_REGION_UPSTREAM |
                        (uint32_t)DUCKVEP_REGION_DOWNSTREAM)) != 0u) {
        uint64_t distance = transcript_projection_abs_diff(
            feature_start1, transcripts->start1[tx_idx]);
        uint64_t candidate = transcript_projection_abs_diff(
            feature_start1, transcripts->end1[tx_idx]);
        if (candidate < distance) distance = candidate;
        candidate = transcript_projection_abs_diff(
            feature_end1, transcripts->start1[tx_idx]);
        if (candidate < distance) distance = candidate;
        candidate = transcript_projection_abs_diff(
            feature_end1, transcripts->end1[tx_idx]);
        if (candidate < distance) distance = candidate;
        out->transcript_distance = distance;
        out->has_transcript_distance = 1u;
    }

    if (coding_context != NULL && raw_cds_start != 0u &&
        raw_cds_end != 0u && duckvep_coding_context_peptide_window_open(
            coding_context, &out->coding_window)) {
        out->coding_context = coding_context;
        out->changed_codon_offset = (uint8_t)((raw_cds_start - 1u) % 3u);
        out->has_coding_window = 1u;
        out->has_amino_acids = (uint8_t)(
            duckvep_feature_allele_peptide_eligible(
                feature_ref, feature_ref_length) &&
            duckvep_feature_allele_peptide_eligible(
                feature_alt, feature_alt_length));
    }
    return 1;
}

uint8_t duckvep_transcript_projection_output_allele_base(
    const duckvep_transcript_projection_facts_t *facts, size_t index) {

    if (facts == NULL || facts->output_allele == NULL ||
        index >= facts->output_allele_length) return 0u;
    return duckvep_event_ascii_upper(facts->output_allele[index]);
}

uint8_t duckvep_transcript_projection_amino_acid_base(
    const duckvep_transcript_projection_facts_t *facts, int alternate,
    size_t index) {

    size_t length;
    if (facts == NULL || !facts->has_coding_window ||
        !facts->has_amino_acids || facts->coding_context == NULL) return 0u;
    length = alternate ? facts->coding_window.alt_length
                       : facts->coding_window.ref_length;
    if (index >= length) return 0u;
    return duckvep_coding_context_peptide_window_base(
        facts->coding_context, &facts->coding_window, alternate, index);
}

uint8_t duckvep_transcript_projection_codon_base(
    const duckvep_transcript_projection_facts_t *facts, int alternate,
    size_t index) {

    size_t length;
    size_t peptide_offset;
    size_t nucleotide_offset;
    size_t changed_length;
    uint8_t base;

    if (facts == NULL || !facts->has_coding_window ||
        facts->coding_context == NULL) return 0u;
    length = alternate ? facts->coding_window.alt_nt_length
                       : facts->coding_window.ref_nt_length;
    peptide_offset = alternate ? facts->coding_window.alt_peptide_offset
                               : facts->coding_window.ref_peptide_offset;
    if (index >= length || peptide_offset > SIZE_MAX / 3u ||
        index > SIZE_MAX - peptide_offset * 3u) return 0u;
    nucleotide_offset = peptide_offset * 3u + index;
    base = (uint8_t)duckvep_coding_context_codon_base(
        facts->coding_context, alternate, nucleotide_offset);
    if (base == 0u) return 0u;
    changed_length = alternate ? (size_t)facts->feature_alt_length
                               : (size_t)facts->feature_ref_length;
    if (index < (size_t)facts->changed_codon_offset ||
        index - (size_t)facts->changed_codon_offset >= changed_length) {
        if (base >= (uint8_t)'A' && base <= (uint8_t)'Z') {
            base = (uint8_t)(base + ((uint8_t)'a' - (uint8_t)'A'));
        }
    }
    return base;
}

duckvep_transcript_edit_status_t duckvep_transcript_edit_project_prepared_hint(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_variant_batch_t    *variants,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                          exon_hint,
    duckvep_transcript_edit_t        *out) {

    duckvep_transcript_edit_t result;
    duckvep_transcript_edit_status_t status;
    uint32_t first_pos1;
    uint32_t last_pos1;
    uint32_t feature_first_pos1;
    uint32_t feature_last_pos1;
    uint32_t transcript_start1;
    uint32_t transcript_end1;
    size_t exon_offset;
    size_t exon_count;
    int8_t strand;

    if (out == NULL) return DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || variants == NULL ||
        event == NULL ||
        tx_idx >= transcripts->transcript_count ||
        transcripts->strand == NULL ||
        transcripts->start1 == NULL ||
        transcripts->end1 == NULL ||
        (uintmax_t)tx_idx > (uintmax_t)UINT32_MAX ||
        (size_t)variant_idx >= variants->count ||
        variants->variant_kind == NULL ||
        variants->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_SV ||
        event->kind != variants->variant_kind[variant_idx] ||
        event->kind == (uint8_t)DUCKVEP_KIND_SV ||
        !duckvep_event_allele_slices_ok(variants, (size_t)variant_idx)) {
        return variants != NULL && variants->variant_kind != NULL &&
               (size_t)variant_idx < variants->count &&
               variants->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_SV
            ? DUCKVEP_TRANSCRIPT_EDIT_UNSUPPORTED_KIND
            : DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG;
    }
    memset(&result, 0, sizeof result);
    result.variant_idx = variant_idx;
    result.tx_idx = (uint32_t)tx_idx;
    result.event = *event;
    if (result.event.ref_diff_offset > variants->ref_length[variant_idx] ||
        result.event.ref_diff_length >
            variants->ref_length[variant_idx] - result.event.ref_diff_offset ||
        result.event.alt_diff_offset > variants->alt_length[variant_idx] ||
        result.event.alt_diff_length >
            variants->alt_length[variant_idx] - result.event.alt_diff_offset ||
        (result.event.interbase != 0u
             ? (result.event.ref_diff_length != 0u ||
                result.event.alt_diff_length == 0u)
             : result.event.ref_diff_length == 0u)) {
        return DUCKVEP_TRANSCRIPT_EDIT_INVALID_EVENT;
    }
    result.raw_ref = variants->allele_bytes + variants->ref_offset[variant_idx];
    result.raw_alt = variants->allele_bytes + variants->alt_offset[variant_idx];
    result.raw_ref_length = variants->ref_length[variant_idx];
    result.raw_alt_length = variants->alt_length[variant_idx];
    if (!duckvep_event_feature_alleles(
            variants, (size_t)variant_idx, &result.event,
            &result.feature_ref, &result.feature_ref_length,
            &result.feature_alt, &result.feature_alt_length)) {
        return DUCKVEP_TRANSCRIPT_EDIT_INVALID_EVENT;
    }
    if (result.event.ref_diff_length == 0u &&
        result.event.alt_diff_length == 0u) {
        return DUCKVEP_TRANSCRIPT_EDIT_INVALID_EVENT;
    }
    result.ref = variants->allele_bytes + variants->ref_offset[variant_idx] +
                 result.event.ref_diff_offset;
    result.alt = variants->allele_bytes + variants->alt_offset[variant_idx] +
                 result.event.alt_diff_offset;
    result.ref_length = result.event.ref_diff_length;
    result.alt_length = result.event.alt_diff_length;
    strand = transcripts->strand[tx_idx];
    result.transcript_strand = strand;
    transcript_start1 = transcripts->start1[tx_idx];
    transcript_end1 = transcripts->end1[tx_idx];
    exon_offset = (size_t)transcripts->exon_offset[tx_idx];
    exon_count = (size_t)transcripts->exon_count[tx_idx];
    if (exon_offset > exons->exon_count ||
        exon_count > exons->exon_count - exon_offset ||
        (size_t)exon_hint < exon_offset ||
        (size_t)exon_hint >= exon_offset + exon_count) {
        exon_hint = UINT32_MAX;
    }

    if (result.event.interbase) {
        first_pos1 = result.event.insertion_boundary0;
        last_pos1 = duckvep_event_right_flank1(&result.event);
        if (first_pos1 == 0u || last_pos1 == 0u ||
            first_pos1 < transcript_start1 ||
            last_pos1 > transcript_end1) {
            return DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT;
        }
    } else {
        if (result.event.end1 < transcript_start1 ||
            result.event.start1 > transcript_end1) {
            /* Allele minimization can move the semantic differing base just
             * outside the transcript while the complete uploaded feature
             * still overlaps its terminal base. VEP clamps that complete
             * feature before deriving transcript-slice REF (terminal CG>CC
             * is the observable c.*10dup case), so retain the terminal
             * projection for the later clamped-feature interpreter. */
            if (result.event.feature_start1 == 0u ||
                result.event.feature_end1 < result.event.feature_start1 ||
                result.event.feature_end1 < transcript_start1 ||
                result.event.feature_start1 > transcript_end1) {
                return DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT;
            }
            first_pos1 = result.event.feature_start1 < transcript_start1
                ? transcript_start1 : result.event.feature_start1;
            last_pos1 = result.event.feature_end1 > transcript_end1
                ? transcript_end1 : result.event.feature_end1;
        } else {
            /* TranscriptVariationAllele::_var2transcript_slice_coords clamps a
             * partially overlapping feature before hgvs_variant_notation reads
             * its transcript-slice REF. Preserve that executable mapper-gap
             * behavior in the shared edit instead of rejecting one endpoint. */
            first_pos1 = result.event.start1 < transcript_start1
                ? transcript_start1 : result.event.start1;
            last_pos1 = result.event.end1 > transcript_end1
                ? transcript_end1 : result.event.end1;
        }
    }
    status = transcript_edit_project_pair(
        transcripts, exons, tx_idx, first_pos1, last_pos1, strand,
        exon_hint, &result.first, &result.last);
    if (status != DUCKVEP_TRANSCRIPT_EDIT_OK) return status;

    if (result.event.interbase) {
        feature_first_pos1 = result.event.insertion_boundary0;
        feature_last_pos1 = duckvep_event_right_flank1(&result.event);
    } else {
        if (result.event.feature_end1 < transcript_start1 ||
            result.event.feature_start1 > transcript_end1) {
            return DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT;
        }
        feature_first_pos1 = result.event.feature_start1 < transcript_start1
            ? transcript_start1 : result.event.feature_start1;
        feature_last_pos1 = result.event.feature_end1 > transcript_end1
            ? transcript_end1 : result.event.feature_end1;
    }
    if (feature_first_pos1 == first_pos1 && feature_last_pos1 == last_pos1) {
        result.feature_first = result.first;
        result.feature_last = result.last;
    } else {
        status = transcript_edit_project_pair(
            transcripts, exons, tx_idx, feature_first_pos1, feature_last_pos1,
            strand, exon_hint, &result.feature_first, &result.feature_last);
        if (status != DUCKVEP_TRANSCRIPT_EDIT_OK) return status;
    }

    *out = result;
    return DUCKVEP_TRANSCRIPT_EDIT_OK;
}

duckvep_transcript_edit_status_t duckvep_transcript_edit_project_prepared(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_variant_batch_t    *variants,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    duckvep_transcript_edit_t        *out) {

    return duckvep_transcript_edit_project_prepared_hint(
        transcripts, exons, variants, variant_idx, tx_idx, event,
        UINT32_MAX, out);
}

duckvep_cds_edit_status_t duckvep_transcript_edit_cds_fill_prepared(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *variants,
    duckvep_haplotype_edit_t         *cds_scratch,
    size_t                            cds_scratch_cap,
    duckvep_transcript_edit_t        *edit) {

    duckvep_cds_edit_status_t status;
    uint32_t cds_exon_hint;

    if (edit == NULL || edit->tx_idx >=
            (transcripts != NULL ? transcripts->transcript_count : 0u) ||
        variants == NULL || (size_t)edit->variant_idx >= variants->count) {
        return DUCKVEP_CDS_EDIT_INVALID_ARG;
    }
    cds_exon_hint = edit->first.exonic && edit->last.exonic &&
            edit->first.exon_idx == edit->last.exon_idx
        ? edit->first.exon_idx : UINT32_MAX;
    status = duckvep_variant_cds_edit_set_build_prepared(
        transcripts, exons, seq, variants, edit->variant_idx, edit->tx_idx,
        edit->transcript_strand, &edit->event, cds_exon_hint, cds_scratch,
        cds_scratch_cap, &edit->cds_edits);
    edit->cds_status = (uint8_t)status;
    edit->cds_built = 1u;
    return status;
}

duckvep_transcript_edit_status_t duckvep_transcript_edit_build_prepared(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *variants,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    duckvep_haplotype_edit_t         *cds_scratch,
    size_t                            cds_scratch_cap,
    duckvep_transcript_edit_t        *out) {

    duckvep_transcript_edit_status_t status;

    status = duckvep_transcript_edit_project_prepared(
        transcripts, exons, variants, variant_idx, tx_idx, event, out);
    if (status != DUCKVEP_TRANSCRIPT_EDIT_OK) return status;
    (void)duckvep_transcript_edit_cds_fill_prepared(
        transcripts, exons, seq, variants, cds_scratch, cds_scratch_cap, out);
    return DUCKVEP_TRANSCRIPT_EDIT_OK;
}

duckvep_transcript_edit_status_t duckvep_transcript_edit_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *variants,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    duckvep_haplotype_edit_t         *cds_scratch,
    size_t                            cds_scratch_cap,
    duckvep_transcript_edit_t        *out) {

    duckvep_event_t event;

    if (out == NULL) return DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (variants == NULL || (size_t)variant_idx >= variants->count ||
        variants->variant_kind == NULL || variants->chrom_id == NULL ||
        variants->pos1 == NULL || variants->end1 == NULL ||
        variants->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_SV ||
        !duckvep_event_allele_slices_ok(variants, (size_t)variant_idx)) {
        return variants != NULL && variants->variant_kind != NULL &&
               (size_t)variant_idx < variants->count &&
               variants->variant_kind[variant_idx] == (uint8_t)DUCKVEP_KIND_SV
            ? DUCKVEP_TRANSCRIPT_EDIT_UNSUPPORTED_KIND
            : DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG;
    }
    duckvep_event_load(variants, (size_t)variant_idx, &event);
    return duckvep_transcript_edit_build_prepared(
        transcripts, exons, seq, variants, variant_idx, tx_idx, &event,
        cds_scratch, cds_scratch_cap, out);
}
