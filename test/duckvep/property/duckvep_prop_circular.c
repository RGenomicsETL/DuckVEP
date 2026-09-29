/* Lifted-interval execution of circular regions.
 *
 * A random circular world (reference, plus- and minus-strand coding and
 * non-coding transcripts, ranked exons, regulatory features, events) is
 * described once in unrolled coordinates. A frame rotates that world by k
 * bases: reference, model and events all move, so objects and events cross the
 * origin at different places. Executing the lifted model must give identical
 * consequences, projected positions, peptides and NMD for every frame once
 * rows are keyed by (event, source object), and one row per pair. Away from the
 * origin the lifted result must equal the ordinary linear kernel. */
#include "duckvep_property.h"
#include "duckvep_lift.h"

#define CP_LMAX      420u
#define CP_MAXTX     5u
#define CP_MAXEXON   4u
#define CP_MAXFEAT   4u
#define CP_MAXEV     28u
#define CP_MAXALLELE 40u
#define CP_HALO      50u
#define CP_MAXROWS   4096u
#define CP_CDNA_MAX  (CP_MAXEXON * 64u)

struct cp_exon { uint32_t u_start; uint32_t length; };

struct cp_tx {
    int8_t   strand;
    uint32_t nex;
    struct cp_exon ex[CP_MAXEXON];
    uint32_t cdna_length;
    uint8_t  cdna[CP_CDNA_MAX];
    int      coding;
    uint32_t cds_a, cds_b; /* 1-based cDNA CDS endpoints */
    uint64_t flags;
    int      mirna;
    int      edit;
};

struct cp_feature { uint32_t u_start; uint32_t length; uint8_t kind; };

struct cp_event {
    uint32_t u_pos;      /* unrolled 1-based start */
    uint32_t ref_length;
    uint32_t alt_length;
    uint8_t  alt[CP_MAXALLELE];
};

struct cp_world {
    uint32_t L;
    uint8_t ref[CP_LMAX];
    uint32_t ntx;
    struct cp_tx tx[CP_MAXTX];
    uint32_t nfeat;
    struct cp_feature feat[CP_MAXFEAT];
    uint32_t nev;
    struct cp_event ev[CP_MAXEV];
    uint64_t seed;
};

struct cp_row {
    uint32_t event, kind, object;
    uint64_t mask;
    uint32_t region, flags;
    uint8_t impact, status, nmd, nmd_escape, aa_ref, aa_alt;
    int32_t cdna, cds, protein;
};

struct cp_frame {
    uint32_t k;
    int guard;
    /* transcripts */
    uint16_t tchrom[CP_MAXTX];
    uint32_t tstart[CP_MAXTX], tend[CP_MAXTX];
    int8_t   tstrand[CP_MAXTX];
    uint64_t tflags[CP_MAXTX];
    uint32_t texoff[CP_MAXTX];
    uint16_t texcnt[CP_MAXTX];
    uint32_t tcds_s[CP_MAXTX], tcds_e[CP_MAXTX];
    uint32_t xstart[CP_MAXTX * CP_MAXEXON], xend[CP_MAXTX * CP_MAXEXON];
    uint32_t xcs[CP_MAXTX * CP_MAXEXON], xce[CP_MAXTX * CP_MAXEXON];
    int8_t   xphase[CP_MAXTX * CP_MAXEXON], xend_phase[CP_MAXTX * CP_MAXEXON];
    uint32_t mirna_off[CP_MAXTX + 1u], mirna_s[CP_MAXTX], mirna_e[CP_MAXTX];
    uint32_t edit_off[CP_MAXTX + 1u], edit_pos[CP_MAXTX];
    uint8_t  edit_alt[CP_MAXTX];
    uint8_t  cds_bytes[CP_MAXTX * CP_CDNA_MAX];
    uint64_t cds_off[CP_MAXTX];
    uint32_t cds_len[CP_MAXTX];
    uint8_t  table[CP_MAXTX];
    uint8_t  flank_bytes[CP_MAXTX * CP_CDNA_MAX];
    uint64_t pre_off[CP_MAXTX], post_off[CP_MAXTX];
    uint32_t pre_len[CP_MAXTX], post_len[CP_MAXTX];
    size_t   cds_total, flank_total;
    uint16_t fchrom[CP_MAXFEAT + 1u];
    uint32_t fstart[CP_MAXFEAT + 1u], fend[CP_MAXFEAT + 1u];
    uint8_t  fkind[CP_MAXFEAT + 1u];
    duckvep_transcript_model_t       tx;
    duckvep_exon_model_t             exons;
    duckvep_sequence_pool_t          seq;
    duckvep_interval_feature_model_t features;
    /* events, sorted by rotated start */
    uint16_t vchrom[CP_MAXEV];
    uint32_t vpos[CP_MAXEV], vend[CP_MAXEV];
    uint32_t vroff[CP_MAXEV], vaoff[CP_MAXEV];
    uint16_t vrlen[CP_MAXEV], valen[CP_MAXEV];
    uint8_t  vkind[CP_MAXEV];
    uint8_t  vbytes[CP_MAXEV * 2u * CP_MAXALLELE];
    uint32_t vid[CP_MAXEV];
    size_t   vbytes_len;
    uint32_t lpos[CP_MAXEV], lend[CP_MAXEV];
};

static uint64_t cp_next(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));

    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static uint32_t cp_range(uint64_t *state, uint32_t lo, uint32_t hi) {
    return lo + (uint32_t)(cp_next(state) % (uint64_t)(hi - lo + 1u));
}

static const char cp_bases[4] = {'A', 'C', 'G', 'T'};

static uint8_t cp_comp(uint8_t base) {
    switch (base) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    default: return 'A';
    }
}

static uint8_t cp_base_at(const struct cp_world *w, uint32_t u) {
    return w->ref[(u - 1u) % w->L];
}

/* Source position of unrolled coordinate u in the frame rotated by k. */
static uint32_t cp_pos(const struct cp_world *w, uint32_t k, uint32_t u) {
    return (u - 1u + k) % w->L + 1u;
}

static uint32_t cp_tx_u(const struct cp_tx *t, uint32_t cdna) {
    uint32_t r, cum = 0u;

    for (r = 0u; r < t->nex; r++) {
        uint32_t j = t->strand > 0 ? r : t->nex - 1u - r;
        uint32_t len = t->ex[j].length;

        if (cdna <= cum + len) {
            uint32_t off = cdna - cum - 1u;
            return t->strand > 0 ? t->ex[j].u_start + off
                                 : t->ex[j].u_start + len - 1u - off;
        }
        cum += len;
    }
    return 0u;
}

/* mode 0: general world; 1: fits one lap so frame 0 needs no lift; 2: a small
 * circle so both images of an object are inside the flank window and ties
 * between images occur. */
static void cp_make_world(struct cp_world *w, uint64_t seed, int mode) {
    int linear_frame0 = mode == 1;
    uint64_t s = seed;
    uint32_t i, j, cursor;

    memset(w, 0, sizeof *w);
    w->seed = seed;
    w->L = linear_frame0 ? CP_LMAX : mode == 2 ? cp_range(&s, 100u, 140u)
                                               : cp_range(&s, 160u, CP_LMAX);
    for (i = 0u; i < w->L; i++) w->ref[i] = (uint8_t)cp_bases[cp_next(&s) & 3u];
    /* Objects tile the circle from a random start, so some cross the cut of
     * the frame they are described in unless linear_frame0 keeps them inside
     * [1, L]. */
    w->ntx = mode == 2 ? cp_range(&s, 1u, 2u) : cp_range(&s, 2u, linear_frame0 ? 3u : CP_MAXTX);
    cursor = linear_frame0 ? 1u + CP_HALO : cp_range(&s, 1u, w->L);
    for (i = 0u; i < w->ntx; i++) {
        struct cp_tx *t = &w->tx[i];
        uint32_t cum = 0u;

        t->strand = (cp_next(&s) & 1u) ? 1 : -1;
        t->nex = cp_range(&s, 1u, mode == 2 ? 2u : 3u);
        for (j = 0u; j < t->nex; j++) {
            t->ex[j].u_start = cursor;
            t->ex[j].length = mode == 2 ? cp_range(&s, 10u, 18u) : cp_range(&s, 14u, 40u);
            cursor += t->ex[j].length + (mode == 2 ? cp_range(&s, 6u, 10u) : cp_range(&s, 10u, 24u));
            cum += t->ex[j].length;
        }
        cursor += cp_range(&s, 8u, 20u);
        t->cdna_length = cum;
        {
            /* transcript-oriented sequence, exons in rank order */
            uint32_t r, out = 0u;

            for (r = 0u; r < t->nex; r++) {
                uint32_t jj = t->strand > 0 ? r : t->nex - 1u - r;
                uint32_t p;

                if (t->strand > 0) {
                    for (p = 0u; p < t->ex[jj].length; p++) {
                        t->cdna[out++] = cp_base_at(w, t->ex[jj].u_start + p);
                    }
                } else {
                    for (p = t->ex[jj].length; p-- > 0u;) {
                        t->cdna[out++] = cp_comp(cp_base_at(w, t->ex[jj].u_start + p));
                    }
                }
            }
        }
        t->coding = (cp_next(&s) % 10u) < 7u && cum >= 24u;
        if (t->coding) {
            uint32_t codons = cp_range(&s, 3u, (cum - 6u) / 3u);

            t->cds_a = cp_range(&s, 1u, cum - codons * 3u + 1u);
            t->cds_b = t->cds_a + codons * 3u - 1u;
            t->flags = (uint64_t)DUCKVEP_TX_HAS_TRANSLATION |
                       (uint64_t)DUCKVEP_TX_BIOTYPE_PROTEIN_CODING;
            if (cp_next(&s) % 5u == 0u) t->flags |= (uint64_t)DUCKVEP_TX_BIOTYPE_NMD;
            t->edit = codons >= 4u && (cp_next(&s) % 3u == 0u);
        } else if (t->nex == 1u && t->ex[0].length >= 16u && cp_next(&s) % 3u == 0u) {
            t->flags = (uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA;
            t->mirna = 1;
        }
    }
    w->nfeat = mode == 2 ? cp_range(&s, 0u, 1u) : cp_range(&s, 1u, 3u);
    for (i = 0u; i < w->nfeat; i++) {
        w->feat[i].u_start = cursor;
        w->feat[i].length = cp_range(&s, 6u, 30u);
        w->feat[i].kind = (uint8_t)((cp_next(&s) & 1u) ? DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE
                                    : DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION);
        cursor += w->feat[i].length + cp_range(&s, 5u, 20u);
    }
    if (mode == 2 && cursor > w->L + 1u) {
        w->L = 0u;
        return;
    }
    if (linear_frame0 && cursor + CP_HALO + CP_MAXALLELE > w->L) {
        /* The generated world does not fit inside one lap of the circle. */
        w->L = 0u;
        return;
    }
    w->nev = cp_range(&s, 12u, CP_MAXEV);
    for (i = 0u; i < w->nev; i++) {
        struct cp_event *e = &w->ev[i];
        uint32_t kind = cp_range(&s, 0u, 5u), n;

        if (cp_next(&s) % 10u < 7u) {
            /* Aim at an exon, including across its boundaries. */
            const struct cp_tx *t = &w->tx[cp_next(&s) % w->ntx];
            const struct cp_exon *x = &t->ex[cp_next(&s) % t->nex];
            uint32_t lo = x->u_start > 6u ? x->u_start - 6u : 1u;

            e->u_pos = cp_range(&s, lo, x->u_start + x->length + 4u);
            if (!linear_frame0) e->u_pos = (e->u_pos - 1u) % w->L + 1u;
        } else {
            e->u_pos = linear_frame0 ? cp_range(&s, 1u + CP_HALO, w->L - CP_HALO - CP_MAXALLELE)
                                     : cp_range(&s, 1u, w->L);
        }
        if (linear_frame0 && (e->u_pos < 1u + CP_HALO ||
                              e->u_pos > w->L - CP_HALO - CP_MAXALLELE)) {
            e->u_pos = 1u + CP_HALO;
        }
        if (kind <= 1u) { /* SNV */
            uint8_t ref = cp_base_at(w, e->u_pos), alt;

            e->ref_length = 1u;
            do { alt = (uint8_t)cp_bases[cp_next(&s) & 3u]; } while (alt == ref);
            e->alt[0] = alt; e->alt_length = 1u;
        } else if (kind == 2u) { /* deletion, anchor base kept */
            e->ref_length = cp_range(&s, 2u, 30u);
            e->alt[0] = cp_base_at(w, e->u_pos); e->alt_length = 1u;
        } else if (kind == 3u) { /* insertion after the anchor */
            e->ref_length = 1u;
            n = cp_range(&s, 1u, 18u);
            e->alt[0] = cp_base_at(w, e->u_pos);
            for (j = 0u; j < n; j++) e->alt[1u + j] = (uint8_t)cp_bases[cp_next(&s) & 3u];
            e->alt_length = n + 1u;
        } else if (kind == 4u) { /* insertion before the anchor */
            e->ref_length = 1u;
            n = cp_range(&s, 1u, 18u);
            for (j = 0u; j < n; j++) e->alt[j] = (uint8_t)cp_bases[cp_next(&s) & 3u];
            e->alt[n] = cp_base_at(w, e->u_pos);
            e->alt_length = n + 1u;
        } else { /* MNV or replacement, up to long alleles */
            e->ref_length = cp_range(&s, 2u, 30u);
            e->alt_length = cp_range(&s, 1u, 30u);
            for (j = 0u; j < e->alt_length; j++) {
                e->alt[j] = (uint8_t)cp_bases[cp_next(&s) & 3u];
            }
            if (e->alt_length == e->ref_length) {
                if (e->alt[0] == cp_base_at(w, e->u_pos)) {
                    e->alt[0] = cp_comp(e->alt[0]);
                }
            }
        }
    }
    /* The base-frame oracle needs no event crossing its own window. */
}

static int cp_event_cmp_pos(const void *a, const void *b) {
    const uint32_t *x = a, *y = b;

    return x[0] < y[0] ? -1 : x[0] > y[0];
}

/* Build the borrowed model and event views of one frame. */
static void cp_make_frame(const struct cp_world *w, uint32_t k, int guard,
                          struct cp_frame *f) {
    uint32_t i, j, x = 0u, mirna = 0u, edits = 0u, nf = 0u;
    uint32_t order[CP_MAXEV][2];
    size_t bytes = 0u;

    memset(f, 0, sizeof *f);
    f->k = k;
    f->guard = guard;
    for (i = 0u; i < w->ntx; i++) {
        const struct cp_tx *t = &w->tx[i];
        uint32_t r, cum = 0u, cds_before = 0u;

        f->tchrom[i] = 0u;
        f->tstrand[i] = t->strand;
        f->tflags[i] = t->flags;
        f->texoff[i] = x;
        f->texcnt[i] = (uint16_t)t->nex;
        f->tstart[i] = cp_pos(w, k, t->ex[0].u_start);
        f->tend[i] = cp_pos(w, k, t->ex[t->nex - 1u].u_start +
                                  t->ex[t->nex - 1u].length - 1u);
        for (r = 0u; r < t->nex; r++, x++) {
            uint32_t jj = t->strand > 0 ? r : t->nex - 1u - r;
            uint32_t us = t->ex[jj].u_start, ue = us + t->ex[jj].length - 1u;
            uint32_t cs = cum + 1u, ce = cum + t->ex[jj].length;

            f->xstart[x] = cp_pos(w, k, us);
            f->xend[x] = cp_pos(w, k, ue);
            f->xcs[x] = cs;
            f->xce[x] = ce;
            f->xphase[x] = -1;
            f->xend_phase[x] = -1;
            if (t->coding && ce >= t->cds_a && cs <= t->cds_b) {
                uint32_t first = cs > t->cds_a ? cs : t->cds_a;
                uint32_t last = ce < t->cds_b ? ce : t->cds_b;

                if (cs >= t->cds_a) f->xphase[x] = (int8_t)(cds_before % 3u);
                if (ce < t->cds_b)
                    f->xend_phase[x] = (int8_t)((cds_before + last - first + 1u) % 3u);
                cds_before += last - first + 1u;
            }
            cum = ce;
        }
        if (t->coding) {
            f->tcds_s[i] = cp_pos(w, k, cp_tx_u(t, t->strand > 0 ? t->cds_a : t->cds_b));
            f->tcds_e[i] = cp_pos(w, k, cp_tx_u(t, t->strand > 0 ? t->cds_b : t->cds_a));
            f->cds_off[i] = f->cds_total;
            f->cds_len[i] = t->cds_b - t->cds_a + 1u;
            memcpy(f->cds_bytes + f->cds_total, t->cdna + t->cds_a - 1u, f->cds_len[i]);
            f->cds_total += f->cds_len[i];
            f->pre_off[i] = f->flank_total;
            f->pre_len[i] = t->cds_a - 1u;
            memcpy(f->flank_bytes + f->flank_total, t->cdna, f->pre_len[i]);
            f->flank_total += f->pre_len[i];
            f->post_off[i] = f->flank_total;
            f->post_len[i] = t->cdna_length - t->cds_b;
            memcpy(f->flank_bytes + f->flank_total, t->cdna + t->cds_b, f->post_len[i]);
            f->flank_total += f->post_len[i];
        } else {
            f->cds_off[i] = f->cds_total;
            f->pre_off[i] = f->post_off[i] = f->flank_total;
        }
        f->table[i] = (uint8_t)STD;
        f->mirna_off[i] = mirna;
        if (t->mirna) {
            uint32_t us = t->ex[0].u_start + 2u;

            f->mirna_s[mirna] = cp_pos(w, k, us);
            f->mirna_e[mirna] = cp_pos(w, k, us + 6u);
            mirna++;
        }
        f->edit_off[i] = edits;
        if (t->edit) {
            f->edit_pos[edits] = 2u;
            f->edit_alt[edits] = (uint8_t)'U';
            edits++;
        }
    }
    f->mirna_off[w->ntx] = mirna;
    f->edit_off[w->ntx] = edits;
    for (i = 0u; i < w->nfeat; i++) {
        f->fchrom[nf] = 0u;
        f->fstart[nf] = cp_pos(w, k, w->feat[i].u_start);
        f->fend[nf] = cp_pos(w, k, w->feat[i].u_start + w->feat[i].length - 1u);
        f->fkind[nf] = w->feat[i].kind;
        nf++;
    }
    if (guard) {
        /* A feature covering all but two bases wraps in every frame whose
         * origin is not adjacent to it, so every frame is lifted. */
        f->fchrom[nf] = 0u;
        f->fstart[nf] = cp_pos(w, k, w->L / 2u);
        f->fend[nf] = cp_pos(w, k, w->L / 2u + w->L - 3u);
        f->fkind[nf] = (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION;
        nf++;
    }
    f->tx.chrom_id = f->tchrom; f->tx.start1 = f->tstart; f->tx.end1 = f->tend;
    f->tx.strand = f->tstrand; f->tx.flags = f->tflags;
    f->tx.exon_offset = f->texoff; f->tx.exon_count = f->texcnt;
    f->tx.cds_start1 = f->tcds_s; f->tx.cds_end1 = f->tcds_e;
    f->tx.transcript_count = w->ntx;
    f->tx.mature_mirna_offset = f->mirna_off;
    f->tx.mature_mirna_start1 = f->mirna_s;
    f->tx.mature_mirna_end1 = f->mirna_e;
    f->tx.mature_mirna_count = mirna;
    f->exons.start1 = f->xstart; f->exons.end1 = f->xend;
    f->exons.cdna_start1 = f->xcs; f->exons.cdna_end1 = f->xce;
    f->exons.phase = f->xphase; f->exons.end_phase = f->xend_phase;
    f->exons.exon_count = x;
    f->seq.cds_bytes = f->cds_bytes; f->seq.cds_bytes_len = f->cds_total;
    f->seq.cds_offset = f->cds_off; f->seq.cds_length = f->cds_len;
    f->seq.codon_table = f->table; f->seq.transcript_count = w->ntx;
    f->seq.peptide_edit_offset = f->edit_off;
    f->seq.peptide_edit_position1 = f->edit_pos;
    f->seq.peptide_edit_alt = f->edit_alt;
    f->seq.peptide_edit_count = edits;
    f->seq.flank_bytes = f->flank_bytes; f->seq.flank_bytes_len = f->flank_total;
    f->seq.pre_cds_offset = f->pre_off; f->seq.pre_cds_length = f->pre_len;
    f->seq.post_cds_offset = f->post_off; f->seq.post_cds_length = f->post_len;
    f->seq.flanks_complete = 1u;
    f->features.chrom_id = f->fchrom; f->features.start1 = f->fstart;
    f->features.end1 = f->fend; f->features.kind = f->fkind;
    f->features.feature_count = nf;

    for (i = 0u; i < w->nev; i++) {
        order[i][0] = cp_pos(w, k, w->ev[i].u_pos);
        order[i][1] = i;
    }
    qsort(order, w->nev, sizeof order[0], cp_event_cmp_pos);
    for (i = 0u; i < w->nev; i++) {
        const struct cp_event *e = &w->ev[order[i][1]];

        f->vid[i] = order[i][1];
        f->vchrom[i] = 0u;
        f->vpos[i] = order[i][0];
        f->vend[i] = f->vpos[i] + e->ref_length - 1u;
        f->vroff[i] = (uint32_t)bytes;
        f->vrlen[i] = (uint16_t)e->ref_length;
        for (j = 0u; j < e->ref_length; j++) f->vbytes[bytes++] = cp_base_at(w, e->u_pos + j);
        f->vaoff[i] = (uint32_t)bytes;
        f->valen[i] = (uint16_t)e->alt_length;
        memcpy(f->vbytes + bytes, e->alt, e->alt_length);
        bytes += e->alt_length;
        {
            duckvep_event_t prepared;

            if (duckvep_event_prepare_small(f->vpos[i], f->vbytes + f->vroff[i],
                                            f->vrlen[i], f->vbytes + f->vaoff[i],
                                            f->valen[i], &prepared)) {
                f->vkind[i] = prepared.kind;
            } else {
                f->vkind[i] = (uint8_t)DUCKVEP_KIND_SNV;
            }
        }
    }
    f->vbytes_len = bytes;
}

static int cp_row_cmp(const void *a, const void *b) {
    const struct cp_row *x = a, *y = b;

    if (x->event != y->event) return x->event < y->event ? -1 : 1;
    if (x->kind != y->kind) return x->kind < y->kind ? -1 : 1;
    return x->object < y->object ? -1 : x->object > y->object;
}

static int cp_row_equal(const struct cp_row *x, const struct cp_row *y) {
    return x->event == y->event && x->kind == y->kind && x->object == y->object &&
           x->mask == y->mask && x->region == y->region && x->flags == y->flags &&
           x->impact == y->impact && x->status == y->status && x->nmd == y->nmd &&
           x->nmd_escape == y->nmd_escape && x->aa_ref == y->aa_ref &&
           x->aa_alt == y->aa_alt && x->cdna == y->cdna && x->cds == y->cds &&
           x->protein == y->protein;
}

/* Annotate one frame. `lift_it` runs the lifted model; otherwise the source
 * arrays feed the linear kernel directly. Returns the canonical row count or
 * -1 on any kernel or lift failure; raw_rows receives the pre-resolve count. */
static int cp_annotate(const struct cp_world *w, struct cp_frame *f, int lift_it,
                       struct cp_row *out, size_t *raw_rows, char *why) {
    duckvep_lift_t *lift = NULL;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_error_t error;
    duckvep_options_init_t oi;
    duckvep_variant_batch_t batch;
    duckvep_result_builder_t builder;
    duckvep_consequence_t *rows;
    size_t *order = NULL, kept = 0u, n, i;
    uint32_t length = 0u, base = 0u, virtual_length = 0u;
    int result = -1;
    uint16_t region_id[1] = {0u};
    uint32_t region_length[1];
    uint8_t region_circular[1] = {1u};
    duckvep_lift_regions_t regions;
    const duckvep_transcript_model_t *tx = &f->tx;
    const duckvep_exon_model_t *ex = &f->exons;
    const duckvep_sequence_pool_t *sq = &f->seq;
    const duckvep_interval_feature_model_t *ft = &f->features;

    rows = malloc(CP_MAXROWS * sizeof *rows);
    order = malloc(CP_MAXROWS * sizeof *order);
    if (rows == NULL || order == NULL) { snprintf(why, 200, "oom"); goto done; }
    memset(&error, 0, sizeof error);
    region_length[0] = w->L;
    regions.chrom_id = region_id; regions.length = region_length;
    regions.circular = region_circular; regions.count = 1u;
    if (lift_it) {
        if (duckvep_lift_open(&regions, &f->tx, &f->exons, &f->seq, &f->features,
                              &lift, &error) != DUCKVEP_OK) {
            snprintf(why, 256, "lift_open: %.200s", error.message);
            goto done;
        }
        tx = &lift->transcripts; ex = &lift->exons; sq = &lift->sequences;
        ft = &lift->interval_features;
        if (!duckvep_lift_region(lift, 0u, &length, &base, &virtual_length)) {
            snprintf(why, 200, "region 0 not lifted");
            goto done;
        }
    }
    if (duckvep_model_open(tx, ex, sq, ft, &model, &error) != DUCKVEP_OK) {
        snprintf(why, 256, "model_open: %.200s", error.message);
        goto done;
    }
    memset(&oi, 0, sizeof oi);
    oi.upstream_dist = CP_HALO; oi.downstream_dist = CP_HALO; oi.halo = CP_HALO;
    oi.distances_are_explicit = 1u;
    if (duckvep_options_open(&oi, &options, &error) != DUCKVEP_OK ||
        duckvep_workspace_open(model, &workspace, &error) != DUCKVEP_OK) {
        snprintf(why, 256, "open: %.200s", error.message);
        goto done;
    }
    for (i = 0u; i < w->nev; i++) {
        f->lpos[i] = f->vpos[i] + base;
        f->lend[i] = f->vend[i] + base;
    }
    memset(&batch, 0, sizeof batch);
    batch.chrom_id = f->vchrom; batch.pos1 = f->lpos; batch.end1 = f->lend;
    batch.ref_offset = f->vroff; batch.ref_length = f->vrlen;
    batch.alt_offset = f->vaoff; batch.alt_length = f->valen;
    batch.allele_bytes = f->vbytes; batch.allele_bytes_len = f->vbytes_len;
    batch.variant_kind = f->vkind; batch.count = w->nev;
    duckvep_result_builder_init(&builder, rows, CP_MAXROWS);
    if (duckvep_annotate_tile(model, &batch, options, workspace, &builder, &error) !=
        DUCKVEP_OK) {
        snprintf(why, 256, "annotate: %.200s", error.message);
        goto done;
    }
    n = duckvep_result_builder_count(&builder);
    *raw_rows = n;
    if (lift_it) {
        if (duckvep_lift_resolve(lift, rows, n, f->lpos, f->lend, order, &kept,
                                 &error) != DUCKVEP_OK) {
            snprintf(why, 256, "resolve: %.200s", error.message);
            goto done;
        }
        n = kept;
    }
    for (i = 0u; i < n; i++) {
        struct cp_row *r = &out[i];
        const duckvep_consequence_t *c = &rows[i];

        r->event = f->vid[c->variant_idx];
        r->kind = c->overlap_object_kind;
        r->object = c->overlap_object_kind == (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT
                        ? c->tx_idx : c->interval_feature_idx;
        r->mask = c->consequence_mask; r->region = c->region_mask; r->flags = c->flags;
        r->impact = c->impact; r->status = c->sequence_status;
        r->nmd = c->nmd_prediction; r->nmd_escape = c->nmd_escape_reasons;
        r->aa_ref = c->aa_ref; r->aa_alt = c->aa_alt;
        r->cdna = c->cdna_pos; r->cds = c->cds_pos; r->protein = c->protein_pos;
    }
    qsort(out, n, sizeof *out, cp_row_cmp);
    result = (int)n;
done:
    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    duckvep_lift_close(lift);
    free(rows);
    free(order);
    return result;
}

static int cp_same(const struct cp_row *a, int na, const struct cp_row *b, int nb) {
    int i;

    if (na != nb) return 0;
    for (i = 0; i < na; i++) if (!cp_row_equal(&a[i], &b[i])) return 0;
    return 1;
}

static struct cp_row cp_rows[3][CP_MAXROWS];
static struct cp_row cp_filtered[CP_MAXROWS];
static uint64_t cp_stat_rows, cp_stat_removed, cp_stat_coding, cp_stat_flank,
    cp_stat_nmd, cp_stat_frames, cp_stat_peptide, cp_stat_resolved;

static enum greatest_test_res cp_rotation(int mode, uint64_t trials, uint64_t seed) {
    uint64_t trial;
    static struct cp_world w;
    static struct cp_frame f;

    for (trial = 0u; trial < trials; trial++) {
        uint64_t s = seed + trial * UINT64_C(0x1000193);
        int n[3], frame;
        size_t raw[3];
        char why[256] = "";

        cp_make_world(&w, cp_next(&s), mode);
        if (w.L == 0u) continue;
        for (frame = 0; frame < 3; frame++) {
            uint32_t k = frame == 0 ? 0u : cp_range(&s, 1u, w.L - 1u);

            /* The guard must wrap in every frame so each one is lifted. */
            while (cp_pos(&w, k, w.L / 2u) <= 3u) k = cp_range(&s, 1u, w.L - 1u);
            if (frame == 0 && cp_pos(&w, k, w.L / 2u) <= 3u) k = 7u;

            cp_make_frame(&w, k, 1, &f);
            n[frame] = cp_annotate(&w, &f, 1, cp_rows[frame], &raw[frame], why);
            if (n[frame] < 0) {
                fprintf(stderr, "seed=%" PRIu64 " trial=%" PRIu64 " frame k=%u: %s\n",
                        seed, trial, k, why);
                FAIL();
            }
            cp_stat_rows += (uint64_t)n[frame];
            cp_stat_removed += raw[frame] - (size_t)n[frame];
            cp_stat_frames++;
            {
                int i;
                for (i = 0; i < n[frame]; i++) {
                    if (cp_rows[frame][i].cds >= 0) cp_stat_coding++;
                    if (cp_rows[frame][i].nmd != 0u) cp_stat_nmd++;
                    if (cp_rows[frame][i].aa_alt != 0u) cp_stat_peptide++;
                    if (cp_rows[frame][i].status == (uint8_t)DUCKVEP_SEQUENCE_RESOLVED) cp_stat_resolved++;
                    if (cp_rows[frame][i].kind == 0u &&
                        (cp_rows[frame][i].region & (DUCKVEP_REGION_UPSTREAM |
                                                     DUCKVEP_REGION_DOWNSTREAM)))
                        cp_stat_flank++;
                    if (i > 0 && cp_row_cmp(&cp_rows[frame][i - 1], &cp_rows[frame][i]) == 0) {
                        fprintf(stderr, "duplicate (event,object) row seed=%" PRIu64
                                " trial=%" PRIu64 "\n", seed, trial);
                        FAIL();
                    }
                }
            }
            if (frame > 0 && !cp_same(cp_rows[0], n[0], cp_rows[frame], n[frame])) {
                int i;
                fprintf(stderr, "rotation mismatch seed=%" PRIu64 " trial=%" PRIu64
                        " L=%u k=%u: %d vs %d rows\n", seed, trial, w.L, k, n[0], n[frame]);
                for (i = 0; i < n[0] && i < n[frame]; i++) {
                    if (!cp_row_equal(&cp_rows[0][i], &cp_rows[frame][i])) {
                        fprintf(stderr, " first diff: event %u kind %u object %u mask %" PRIx64
                                " vs event %u kind %u object %u mask %" PRIx64
                                " cds %d/%d region %x/%x\n",
                                cp_rows[0][i].event, cp_rows[0][i].kind, cp_rows[0][i].object,
                                cp_rows[0][i].mask, cp_rows[frame][i].event,
                                cp_rows[frame][i].kind, cp_rows[frame][i].object,
                                cp_rows[frame][i].mask, cp_rows[0][i].cds,
                                cp_rows[frame][i].cds, cp_rows[0][i].region,
                                cp_rows[frame][i].region);
                        break;
                    }
                }
                FAIL();
            }
        }
    }
    fprintf(stderr, "circular rotation: frames=%" PRIu64 " rows=%" PRIu64
            " duplicate images removed=%" PRIu64 " coding rows=%" PRIu64
            " flank rows=%" PRIu64 " nmd rows=%" PRIu64 " peptide rows=%" PRIu64
            " sequence-resolved rows=%" PRIu64 "\n", cp_stat_frames, cp_stat_rows,
            cp_stat_removed, cp_stat_coding, cp_stat_flank, cp_stat_nmd,
            cp_stat_peptide, cp_stat_resolved);
    ASSERT(cp_stat_removed > 0u);
    ASSERT(cp_stat_coding > 0u);
    ASSERT(cp_stat_flank > 0u);
    ASSERT(cp_stat_peptide > 0u);
    ASSERT(cp_stat_resolved > 0u);
    PASS();
}

/* Rotation equivariance: three random rotations of a random world, all lifted
 * (the guard feature wraps in each), agree on every row. */
TEST circular_lift_rotation_equivariance(void) {
    CHECK_CALL(cp_rotation(0, kprop_env_u64("DUCKVEP_PROP_TRIALS", 2000u),
                           kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED)));
    ASSERT(cp_stat_coding > 0u);
    ASSERT(cp_stat_peptide > 0u);
    PASS();
}

/* A small circle puts both images of every object inside the flank window, so
 * equidistant images and images of one object on both sides of an event occur. */
TEST circular_lift_rotation_equivariance_small_circle(void) {
    CHECK_CALL(cp_rotation(2, kprop_env_u64("DUCKVEP_PROP_TRIALS", 2000u),
                           kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED) ^ UINT64_C(0xa5a5)));
    PASS();
}

/* Away from the origin the lifted model must reproduce the ordinary kernel:
 * frame 0 of a world that fits inside one lap, events kept a halo away from
 * the cut, compared with the linear kernel on the same source arrays. The
 * guard feature only exists in the lifted run and is filtered out. */
TEST circular_lift_matches_linear_kernel_away_from_origin(void) {
    uint64_t trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", 2000u);
    uint64_t seed = kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED) ^ UINT64_C(0x5bd1e995);
    uint64_t trial, compared = 0u;
    static struct cp_world w;
    static struct cp_frame f;

    for (trial = 0u; trial < trials; trial++) {
        uint64_t s = seed + trial * UINT64_C(0x1000193);
        int nl, nc, i, kept = 0;
        size_t raw_l, raw_c;
        char why[256] = "";

        cp_make_world(&w, cp_next(&s), 1);
        if (w.L == 0u) continue;
        cp_make_frame(&w, 0u, 0, &f);
        nl = cp_annotate(&w, &f, 0, cp_rows[0], &raw_l, why);
        if (nl < 0) { fprintf(stderr, "linear: %s\n", why); FAIL(); }
        cp_make_frame(&w, 0u, 1, &f);
        nc = cp_annotate(&w, &f, 1, cp_rows[1], &raw_c, why);
        if (nc < 0) { fprintf(stderr, "lifted: %s\n", why); FAIL(); }
        for (i = 0; i < nc; i++) {
            if (cp_rows[1][i].kind != 0u && cp_rows[1][i].object == w.nfeat) continue;
            cp_filtered[kept++] = cp_rows[1][i];
        }
        if (!cp_same(cp_rows[0], nl, cp_filtered, kept)) {
            fprintf(stderr, "lifted != linear seed=%" PRIu64 " trial=%" PRIu64
                    " L=%u: %d vs %d rows\n", seed, trial, w.L, nl, kept);
            FAIL();
        }
        compared += (uint64_t)nl;
    }
    fprintf(stderr, "circular lift vs linear kernel: rows compared=%" PRIu64 "\n", compared);
    ASSERT(compared > 0u);
    PASS();
}

/* Model contract: a lift needs a circular region with a length, refuses
 * regions that would overflow the lifted interval, and leaves an unwrapped
 * model to the linear kernel. */
TEST circular_lift_contract(void) {
    static struct cp_world w;
    static struct cp_frame f;
    uint16_t id[1] = {0u};
    uint32_t length[1];
    uint8_t circular[1] = {1u};
    duckvep_lift_regions_t regions;
    duckvep_lift_t *lift = NULL;
    duckvep_error_t error;

    memset(&error, 0, sizeof error);
    {
        uint64_t candidate = 12345u;

        do { cp_make_world(&w, candidate++, 1); } while (w.L == 0u);
    }
    cp_make_frame(&w, 0u, 0, &f);
    length[0] = w.L;
    regions.chrom_id = id; regions.length = length; regions.circular = circular;
    regions.count = 1u;
    ASSERT_EQ(DUCKVEP_ERR_UNSUPPORTED,
              duckvep_lift_open(&regions, &f.tx, &f.exons, &f.seq, &f.features,
                                &lift, &error));
    ASSERT(lift == NULL);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_lift_open(NULL, &f.tx, &f.exons, &f.seq, &f.features,
                                &lift, &error));
    cp_make_frame(&w, 5u, 1, &f);
    circular[0] = 0u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_lift_open(&regions, &f.tx, &f.exons, &f.seq, &f.features,
                                &lift, &error));
    circular[0] = 1u;
    length[0] = 1500000000u;
    ASSERT_EQ(DUCKVEP_ERR_OUT_OF_RANGE,
              duckvep_lift_open(&regions, &f.tx, &f.exons, &f.seq, &f.features,
                                &lift, &error));
    length[0] = w.L;
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_lift_open(&regions, &f.tx, &f.exons, &f.seq, &f.features,
                                &lift, &error));
    {
        uint32_t l, b, v;

        ASSERT(duckvep_lift_region(lift, 0u, &l, &b, &v));
        ASSERT_EQ(w.L, l);
        ASSERT_EQ(0u, b % w.L);
        ASSERT(b >= 3u * 1000u);
        ASSERT(v >= 2u * b + 2u * w.L);
        ASSERT(!duckvep_lift_region(lift, 1u, NULL, NULL, NULL));
    }
    ASSERT_EQ(3u * w.ntx, lift->transcripts.transcript_count);
    duckvep_lift_close(lift);
    PASS();
}
