#!/usr/bin/env python3
"""Parity of duckvep_lof_sql against VEP 116 with the LOFTEE plugin (GRCh38, chr21).

Corpus: chr21 records of GIAB HG002 and ClinVar whose DuckVEP annotation touches a
loss-of-function or splice consequence on a protein-coding transcript, plus a
synthetic set on both strands built on real chr21 transcripts. VEP runs from the
pinned image with LOFTEE at a46b502 (scripts/run_vep116_loftee_docker.sh). The
profile leaves the ancestor and PhyloCSF off; GERP is a constant-negative bigwig,
which reduces LOFTEE's GERP-weighted rule to the unweighted 50 bp rule.

The gate is exact agreement of LoF, LoF_filter and LoF_flags for every variant x
transcript pair LOFTEE reports, with missing and extra pairs counted.

Usage: scripts/lof_parity.py --extension build/release/duckvep.duckdb_extension \
  --model-db MODEL.duckdb --fasta GRCh38.fa --cache CACHE_ROOT --loftee LOFTEE_DIR \
  --gerp gerp_const_neg.bw --clinvar CLINVAR.vcf.gz --hg002 HG002.vcf.gz --out DIR
"""
import argparse, csv, gzip, hashlib, json, os, re, subprocess, sys
import duckdb

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHROM = '21'
COMP = str.maketrans('ACGT', 'TGCA')


def rc(s):
    return s.translate(COMP)[::-1]


def read_vcf(path, names):
    out = []
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line[0] == '#':
                continue
            f = line.split('\t', 5)
            if f[0] in names:
                for alt in f[4].split(','):
                    if set(f[3]) <= set('ACGT') and set(alt) <= set('ACGT') and alt != f[3]:
                        out.append((int(f[1]), f[3], alt))
    return out


def chunks(fasta, size=1000000):
    length = offset = line_bases = line_width = None
    for line in open(fasta + '.fai'):
        f = line.split('\t')
        if f[0] == CHROM:
            length, offset, line_bases, line_width = int(f[1]), int(f[2]), int(f[3]), int(f[4])
    with open(fasta, 'rb') as handle:
        handle.seek(offset)
        raw = handle.read(length + (length // line_bases + 1) * (line_width - line_bases))
    seq = raw.replace(b'\n', b'').decode()[:length].upper()
    return seq, [(CHROM, s, min(s + size, length), seq[s:min(s + size, length)]) for s in range(0, length, size)]


def load_model(con, model_db, fasta):
    con.execute("ATTACH '%s' AS m (READ_ONLY)" % model_db)
    con.execute("CREATE VIEW regions AS SELECT * FROM m.model_regions")
    con.execute("CREATE VIEW transcripts AS SELECT * FROM m.model_transcripts")
    con.execute("""CREATE TABLE exons AS SELECT transcript_index, e.exon_start, e.exon_end, e.exon_cdna_start,
        e.exon_cdna_end, e.phase, e.end_phase FROM (SELECT transcript_index, unnest(exons) e FROM m.model_transcripts)""")
    con.execute("""SELECT loaded FROM duckvep_model_load('grch38',
        'SELECT seq_region, sequence_length, seq_region_name FROM regions ORDER BY seq_region',
        'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags,
         cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM transcripts
         ORDER BY seq_region, transcript_start, transcript_index',
        'SELECT * FROM exons ORDER BY transcript_index, exon_start',
        peptide_edit_query := 'SELECT transcript_index, p.protein_position, p.alternate_amino_acid
         FROM (SELECT transcript_index, unnest(peptide_edits) p FROM m.model_transcripts) ORDER BY 1, 2',
        reference_fasta := '%s')""" % fasta)


def synthetic(con, seq, per=3):
    """Variants on real chr21 transcripts. Returns {(pos, ref, alt): tag}."""
    region = con.execute("SELECT seq_region FROM regions WHERE seq_region_name = '%s'" % CHROM).fetchone()[0]
    tx = {}
    for tid, strand, cs, ce, exons in con.execute("""SELECT transcript_stable_id, strand, cds_start, cds_end, exons
            FROM transcripts WHERE seq_region = ? AND transcript_biotype = 'protein_coding' AND cds_start IS NOT NULL
            ORDER BY transcript_stable_id""", [region]).fetchall():
        ex = sorted(((e['exon_start'], e['exon_end']) for e in exons), reverse=strand < 0)
        tx[tid] = (strand, cs, ce, ex)
    out = {}
    B = lambda p: seq[p - 1]

    def put(tag, pos, ref, alt):
        if 1 < pos < len(seq) - 10 and set(ref + alt) <= set('ACGT') and ref != alt:
            out.setdefault((pos, ref, alt), tag)

    def delete(tag, vs, n=1):
        put(tag, vs - 1, seq[vs - 2:vs - 1 + n], B(vs - 1))

    def insert(tag, vs):  # new bases begin at vs
        put(tag, vs - 1, B(vs - 1), B(vs - 1) + ('G' if B(vs - 1) != 'G' else 'T'))

    def snv(tag, pos, alts=None):
        for alt in (alts or [b for b in 'ACGT' if b != B(pos)][:1]):
            put(tag, pos, B(pos), alt)

    def pick(rows, n):
        return rows[::max(1, len(rows) // n)][:n]

    for strand in (1, -1):
        s = '+' if strand > 0 else '-'
        stop = lambda t: tx[t][2] if strand > 0 else tx[t][1]
        # END_TRUNC boundary: variant in the penultimate exon, d bases from its 3' end; stop in the last exon
        rows = [t for t, (st, cs, ce, ex) in tx.items() if st == strand and len(ex) >= 3
                and ex[-1][0] <= stop(t) <= ex[-1][1] and ex[-2][1] - ex[-2][0] > 200
                and cs <= ex[-2][0] and ex[-2][1] <= ce]
        for d in (49, 50, 51):
            for t in pick(rows, per):
                es, ee = tx[t][3][-2]
                vs = ee - d if strand > 0 else es + d
                delete('ENDTRUNC_d%d_%s' % (d, s), vs)
                insert('ENDTRUNC_ins_d%d_%s' % (d, s), vs)
        # last coding exon, and stop exon followed by a 3' UTR exon
        rows = [t for t, (st, cs, ce, ex) in tx.items() if st == strand and len(ex) >= 3
                and ex[-1][0] <= stop(t) <= ex[-1][1] and abs(ex[-1][1] - ex[-1][0]) > 120]
        for t in pick(rows, per):
            es, ee = tx[t][3][-1]
            delete('LAST_EXON_%s' % s, (es + ee) // 2)
            snv('LAST_EXON_snv_%s' % s, (es + ee) // 2)
        rows = [t for t, (st, cs, ce, ex) in tx.items() if st == strand and len(ex) >= 4
                and not ex[-1][0] <= stop(t) <= ex[-1][1]]
        for t in pick(rows, per):
            ex = tx[t][3]
            k = next(i for i, (a, b) in enumerate(ex) if a <= stop(t) <= b)
            if k >= 1:
                a, b = ex[k - 1]
                delete('UTR_EXON_AFTER_STOP_%s' % s, (a + b) // 2)
        # single exon
        rows = [t for t, (st, cs, ce, ex) in tx.items() if st == strand and len(ex) == 1]
        for t in pick(rows, per):
            delete('SINGLE_EXON_%s' % s, (tx[t][1] + tx[t][2]) // 2)
            insert('SINGLE_EXON_ins_%s' % s, (tx[t][1] + tx[t][2]) // 2 + 7)
    # introns, in transcript order with oriented motif
    introns = []
    for t, (strand, cs, ce, ex) in tx.items():
        for k in range(len(ex) - 1):
            lo, hi = (ex[k][1] + 1, ex[k + 1][0] - 1) if strand > 0 else (ex[k + 1][1] + 1, ex[k][0] - 1)
            if hi - lo + 1 < 4:
                continue
            if strand > 0:
                donor, acceptor = lo, hi
                first, last = seq[lo - 1:lo + 1], seq[hi - 2:hi]
            else:
                donor, acceptor = hi, lo
                first, last = rc(seq[hi - 2:hi]), rc(seq[lo - 1:lo + 1])
            introns.append(dict(t=t, strand=strand, lo=lo, hi=hi, size=hi - lo + 1, donor=donor, acceptor=acceptor,
                                first=first, last=last, k=k, n=len(ex)))
    chosen = sorted({size for strand in (1, -1) for size in
                     sorted({i['size'] for i in introns if i['strand'] == strand})[:2]})
    for size in chosen:
        for i in [i for i in introns if i['size'] == size][:4]:
            s = '+' if i['strand'] > 0 else '-'
            snv('INTRON_%d_donor_%s' % (size, s), i['donor'])
            snv('INTRON_%d_acceptor_%s' % (size, s), i['acceptor'])
    for strand in (1, -1):
        s = '+' if strand > 0 else '-'
        mine = [i for i in introns if i['strand'] == strand]
        # NAGNAG: the two acceptor bases; the reference window is oriented to the transcript strand
        def window(p):
            w = seq[p - 5:p + 4]
            return w if strand > 0 else rc(w)
        sites = [(i, p) for i in mine for p in (i['acceptor'], i['acceptor'] - 1 if strand > 0 else i['acceptor'] + 1)]
        hit = [(i, p) for i, p in sites if re.search('AG.AG', window(p))]
        miss = [(i, p) for i, p in sites if not re.search('AG.AG', window(p))]
        for i, p in pick(hit, 6):
            snv('NAGNAG_hit_%s' % s, p)
            delete('NAGNAG_hit_del_%s' % s, p)
            insert('NAGNAG_hit_ins_%s' % s, p)
        for i, p in pick(miss, 6):
            snv('NAGNAG_miss_%s' % s, p)
        # GC donor: the +2 base C>T (transcript strand), other alts and the +1 base
        for i in pick([i for i in mine if i['first'] == 'GC'], 4):
            p = i['donor'] + 1 if strand > 0 else i['donor'] - 1
            snv('GC_DONOR_C2T_%s' % s, p, [b for b in 'ACGT' if b == ('T' if strand > 0 else 'A')])
            snv('GC_DONOR_C2other_%s' % s, p, [b for b in 'ACGT' if b not in (B(p), 'T' if strand > 0 else 'A')])
            snv('GC_DONOR_plus1_%s' % s, i['donor'])
        # non-canonical introns
        for i in pick([i for i in mine if not (i['first'] == 'GT' and i['last'] == 'AG')], 5):
            snv('NONCAN_donor_%s' % s, i['donor'])
            snv('NONCAN_acceptor_%s' % s, i['acceptor'])
        # deletions across an exon-intron boundary, donor side and acceptor side
        for i in pick(mine, 4):
            delete('EXON_EDGE_donor_%s' % s, i['lo'] - 1 if strand > 0 else i['hi'], 2)
            delete('EXON_EDGE_acceptor_%s' % s, i['hi'] if strand > 0 else i['lo'] - 1, 2)
    return out, chosen


def vep_pairs(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            d = json.loads(line)
            f = d['input'].split('\t')
            for t in d.get('transcript_consequences', []):
                if any(k in t for k in ('lof', 'lof_filter', 'lof_flags', 'lof_info')):
                    out[(int(f[1]), f[3], f[4], t['transcript_id'])] = (
                        t.get('lof'), t.get('lof_filter'), t.get('lof_flags'), t.get('lof_info'))
    return out


def run_vep(args, vcf, out_json, params):
    cmd = [os.path.join(ROOT, 'scripts', 'run_vep116_loftee_docker.sh'), args.cache, args.fasta, args.loftee,
           args.gerp, vcf, out_json] + params
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    ap = argparse.ArgumentParser()
    for name in ('extension', 'model-db', 'fasta', 'cache', 'loftee', 'gerp', 'clinvar', 'hg002', 'out'):
        ap.add_argument('--' + name, required=True)
    ap.add_argument('--reuse', action='store_true', help='reuse VEP JSON already in OUT/work')
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    work = os.path.join(args.out, 'work')
    os.makedirs(work, exist_ok=True)

    con = duckdb.connect(config={'allow_unsigned_extensions': 'true'})
    con.execute("LOAD '%s'" % args.extension)
    load_model(con, args.model_db, args.fasta)
    seq, rows = chunks(args.fasta)
    con.execute('CREATE TABLE refchunks(chrom VARCHAR, "start" BIGINT, "end" BIGINT, seq VARCHAR)')
    con.executemany('INSERT INTO refchunks VALUES (?, ?, ?, ?)', rows)
    con.execute("CREATE TABLE gerp_const AS SELECT '%s' AS chrom, 0::BIGINT AS \"start\", %d::BIGINT AS \"end\", -1000.0 AS score"
                % (CHROM, len(seq)))
    region = con.execute("SELECT seq_region FROM regions WHERE seq_region_name = '%s'" % CHROM).fetchone()[0]

    # 1. real candidates: annotate every chr21 record, keep those with a LoF or splice consequence
    real = sorted({v for p, names in ((args.clinvar, ('21',)), (args.hg002, ('chr21',))) for v in read_vcf(p, names)})
    con.execute('CREATE TABLE real_all(position BIGINT, reference VARCHAR, alternate VARCHAR)')
    con.executemany('INSERT INTO real_all VALUES (?, ?, ?)', real)

    def events(table, source):
        con.execute("""CREATE TABLE %s AS SELECT row_number() OVER (ORDER BY position, reference, alternate)::UBIGINT AS event_index,
            %d::UINTEGER AS seq_region, position::UBIGINT AS position, reference, alternate, NULL::UBIGINT AS end_position,
            NULL::VARCHAR AS structural_type, NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region,
            NULL::UBIGINT AS mate_position FROM %s""" % (table, region, source))
    events('real_events', 'real_all')
    con.execute("""CREATE TABLE real_keep AS SELECT DISTINCT e.position, e.reference, e.alternate FROM
        query(duckvep_annotate_sql('real_events', 'grch38', {rich: true})) a JOIN real_events e USING (event_index)
        JOIN transcripts t USING (transcript_index) WHERE t.transcript_biotype = 'protein_coding' AND regexp_matches(a.consequence,
        '(^|&)(stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|start_lost|stop_lost|protein_altering_variant)(&|$)')""")
    corpora = {}
    corpora['real'] = {(p, r, a): 'CORPUS' for p, r, a in con.execute(
        'SELECT position, reference, alternate FROM real_keep').fetchall()}
    corpora['synthetic'], chosen = synthetic(con, seq)
    for name, corpus in corpora.items():
        with open(os.path.join(args.out, 'duckvep_lof_corpus_%s.vcf' % name), 'w') as handle:
            handle.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            for (p, r, a), tag in sorted(corpus.items()):
                handle.write('%s\t%d\t%s\t%s\t%s\t.\t.\t.\n' % (CHROM, p, tag, r, a))
    con.execute('CREATE TABLE all_alleles AS SELECT DISTINCT position, reference, alternate FROM (%s)' % ' UNION ALL '.join(
        "SELECT * FROM (VALUES %s) t(position, reference, alternate)" % ','.join("(%d,'%s','%s')" % k for k in c)
        for c in corpora.values() if c))
    events('events', 'all_alleles')
    con.execute("CREATE TABLE ann AS SELECT * FROM query(duckvep_annotate_projected_sql('events', 'grch38'))")

    # 2. runs: name, corpus, LOFTEE parameters, duckvep_lof_sql options, gate
    runs = [('real', 'real', [], '', True), ('synthetic', 'synthetic', [], '', True)]
    thresholds = sorted({t for size in chosen for t in (size - 1, size, size + 1)})
    for threshold in thresholds:
        touching = ', '.join(str(size) for size in chosen if abs(size - threshold) <= 1)
        runs.append(('synthetic min_intron_size=%d (introns of %s nt)' % (threshold, touching), 'synthetic',
                     ['min_intron_size:%d' % threshold], '{min_intron_size: %d}' % threshold, True))
    runs.append(('real, GERP relation (constant -1000), full info', 'real', [], "{gerp: 'gerp_const'}", True))
    runs.append(('synthetic, GERP relation (constant -1000), full info', 'synthetic', [], "{gerp: 'gerp_const'}", True))
    runs.append(('real, check_complete_cds (informational)', 'real', ['check_complete_cds:true'],
                 '{check_complete_cds: true}', False))

    summary, disagreements, digests, classes, tokens = [], [], [], {}, {}
    for name, corpus, params, options, gate in runs:
        slug = corpus + ''.join('_' + ''.join(ch if ch.isalnum() else '_' for ch in p) for p in params)
        out_json = os.path.join(work, slug + '.json')
        if not (args.reuse and os.path.exists(out_json)):
            run_vep(args, os.path.join(args.out, 'duckvep_lof_corpus_%s.vcf' % corpus), out_json, params)
        vep = vep_pairs(out_json)
        con.execute("""CREATE OR REPLACE TABLE lof AS SELECT * FROM query(duckvep_lof_sql('ann', 'transcripts', 'refchunks'%s))"""
                    % (', ' + options if options else ''))
        keys = corpora[corpus]
        duck = {}
        for p, r, a, tid, lof, flt, flg, info, unc in con.execute("""SELECT e.position, e.reference, e.alternate,
                t.transcript_stable_id, l.lof, l.lof_filter, l.lof_flags, l.lof_info, l.lof_unchecked FROM lof l
                JOIN events e USING (event_index) JOIN transcripts t USING (transcript_index) WHERE l.lof IS NOT NULL""").fetchall():
            if (p, r, a) in keys:
                duck[(p, r, a, tid)] = (lof, flt, flg, info)
        full_info = 'GERP' in name
        def strip(i):
            return None if i is None else (i if full_info else ','.join(x for x in i.split(',') if not x.startswith('GERP_DIST')) or None)
        both = set(vep) & set(duck)
        missing, extra = sorted(set(vep) - set(duck)), sorted(set(duck) - set(vep))
        exact = [k for k in both if vep[k][:3] == duck[k][:3]]
        different = sorted(k for k in both if vep[k][:3] != duck[k][:3])
        info_exact = sum(1 for k in both if strip(vep[k][3]) == strip(duck[k][3]))
        summary.append(dict(run=name, gate='yes' if gate else 'no', variants=len(keys), loftee_pairs=len(vep), duckvep_pairs=len(duck),
                            exact=len(exact), missing=len(missing), extra=len(extra), different=len(different),
                            info_exact=info_exact,
                            lof_hc=sum(1 for v in duck.values() if v[0] == 'HC'), lof_lc=sum(1 for v in duck.values() if v[0] == 'LC')))
        for kind, group in (('missing', missing), ('extra', extra), ('different', different)):
            for k in group:
                v, d = vep.get(k, (None,) * 4), duck.get(k, (None,) * 4)
                if gate:
                    disagreements.append([name, kind, CHROM, k[0], k[1], k[2], k[3],
                                          '|'.join(map(str, v[:3])), '|'.join(map(str, d[:3])), corpora[corpus].get(k[:3], '')])
                split = lambda x: set(x.split(',')) if x else set()
                key = (name, kind, ','.join(sorted(split(v[1]) - split(d[1]))), ','.join(sorted(split(d[1]) - split(v[1]))),
                       ','.join(sorted(split(v[2]) - split(d[2]))), ','.join(sorted(split(d[2]) - split(v[2]))))
                classes[key] = classes.get(key, 0) + 1
        for v in duck.values():
            for field, text in (('filter', v[1]), ('flag', v[2])):
                for token in (text or '').split(','):
                    if token:
                        tokens[(name, field, token)] = tokens.get((name, field, token), 0) + 1
        body = '\n'.join('%s\t%s' % ('\t'.join(map(str, k)), '\t'.join(map(str, duck[k][:3]))) for k in sorted(duck))
        digests.append((name, hashlib.sha256(body.encode()).hexdigest(),
                        hashlib.sha256(open(out_json, 'rb').read()).hexdigest()))
    with open(os.path.join(args.out, 'duckvep_lof_parity_summary.csv'), 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    with open(os.path.join(args.out, 'duckvep_lof_parity_disagreements.tsv'), 'w') as handle:
        handle.write('run\tkind\tchrom\tpos\tref\talt\ttranscript\tloftee_lof_filter_flags\tduckvep_lof_filter_flags\tcorpus_tag\n')
        for row in disagreements:
            handle.write('\t'.join(map(str, row)) + '\n')
    with open(os.path.join(args.out, 'duckvep_lof_parity_classes.tsv'), 'w') as handle:
        handle.write('run\tkind\tloftee_only_filters\tduckvep_only_filters\tloftee_only_flags\tduckvep_only_flags\tpairs\n')
        for key, n in sorted(classes.items()):
            handle.write('\t'.join(key) + '\t%d\n' % n)
    with open(os.path.join(args.out, 'duckvep_lof_parity_tokens.tsv'), 'w') as handle:
        handle.write('run\tfield\ttoken\tpairs\n')
        for key, n in sorted(tokens.items()):
            handle.write('\t'.join(key) + '\t%d\n' % n)
    with open(os.path.join(args.out, 'duckvep_lof_parity_digests.tsv'), 'w') as handle:
        handle.write('run\tduckvep_pairs_sha256\tloftee_json_sha256\n')
        for row in digests:
            handle.write('\t'.join(row) + '\n')
    for row in summary:
        print(row)
    failed = [r for r in summary if r['gate'] == 'yes' and (r['missing'] or r['extra'] or r['different'] or
                                                             (r['info_exact'] != r['exact'] + r['different']))]
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
