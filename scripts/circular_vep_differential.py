#!/usr/bin/env python3
"""Executable VEP differential on public circular genomes with an origin-crossing transcript.

Each target is a public Ensembl Genomes release-63 (Ensembl 116) core database region that has the
`circular_seq` attribute and a transcript stored with seq_region_start > seq_region_end, for which the
matching Ensembl Genomes 63 VEP cache exists, so the pinned VEP image can annotate the same events offline:

  carrot_pt      Daucus carota plastid Pt (155,848 bp); KZM81246, trans-spliced rps12, minus strand
  chlamydia      Chlamydia trachomatis D/UW-3/CX chromosome (1,042,519 bp); AAC68473, plus strand
  nanoarchaeum   Nanoarchaeum equitans Kin4-M chromosome (490,885 bp); AAR38856, minus strand

    scripts/circular_vep_differential.py TARGET WORKDIR EXTENSION

Requires network for the first run (Ensembl FTP), samtools, docker with the pinned VEP image, and the
`duckdb` Python module. The extension is copied to WORKDIR first and loaded from that immutable copy.
Prints one JSON receipt.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import random
import shutil
import subprocess
import sys

import duckdb

EG = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/"
TARGETS = {
    "carrot_pt": dict(
        core=EG + "plants/mysql/daucus_carota_core_63_116_1/",
        fasta=EG + "plants/fasta/daucus_carota/dna/Daucus_carota.ASM162521v1.dna.toplevel.fa.gz",
        cache=EG + "plants/variation/indexed_vep_cache/daucus_carota_vep_63_ASM162521v1.tar.gz",
        species="daucus_carota", assembly="ASM162521v1", region_id=4825, region="Pt", wrapped="KZM81246",
        coord_system_id=2),
    "chlamydia": dict(
        core=EG + "bacteria/mysql/bacteria_0_collection_core_63_116_1/",
        fasta=EG + "bacteria/fasta/bacteria_0_collection/chlamydia_trachomatis_d_uw_3_cx_gca_000008725/dna/"
               "Chlamydia_trachomatis_d_uw_3_cx_gca_000008725.ASM872v1.dna.toplevel.fa.gz",
        cache=EG + "bacteria/variation/indexed_vep_cache/bacteria_0_collection/"
               "chlamydia_trachomatis_d_uw_3_cx_gca_000008725_vep_63_ASM872v1.tar.gz",
        species="chlamydia_trachomatis_d_uw_3_cx_gca_000008725", assembly="ASM872v1", region_id=123,
        region="Chromosome", wrapped="AAC68473", coord_system_id=None),
    "nanoarchaeum": dict(
        core=EG + "bacteria/mysql/bacteria_0_collection_core_63_116_1/",
        fasta=EG + "bacteria/fasta/bacteria_0_collection/nanoarchaeum_equitans_kin4_m_gca_000008085/dna/"
               "Nanoarchaeum_equitans_kin4_m_gca_000008085.ASM808v1.dna.toplevel.fa.gz",
        cache=EG + "bacteria/variation/indexed_vep_cache/bacteria_0_collection/"
               "nanoarchaeum_equitans_kin4_m_gca_000008085_vep_63_ASM808v1.tar.gz",
        species="nanoarchaeum_equitans_kin4_m_gca_000008085", assembly="ASM808v1", region_id=626,
        region="Chromosome", wrapped="AAR38856", coord_system_id=None),
}
TABLES = ["attrib_type", "coord_system", "seq_region", "seq_region_attrib", "gene", "transcript",
          "transcript_attrib", "translation", "translation_attrib", "exon", "exon_transcript"]
HERE = os.path.dirname(os.path.abspath(__file__))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url, path):
    if not os.path.exists(path):
        subprocess.run(["curl", "-sS", "--fail", "--retry", "8", "--retry-delay", "3", "-o", path, url], check=True)
    return path


def stage(target, work, extension):
    """Core tables of one region, its reference and the native builder output."""
    core = target["core"]
    ddl = gzip.open(fetch(core + os.path.basename(core.rstrip("/")) + ".sql.gz", os.path.join(work, "core.sql.gz")),
                    "rt").read()
    columns = {}
    for table in TABLES:
        body = re.search(r"CREATE TABLE `%s` \((.*?)\n\) ENGINE" % table, ddl, re.S).group(1)
        columns[table] = [(m.group(1), m.group(2)) for m in
                          (re.match(r"`(\w+)`\s+(\S+)", line.strip()) for line in body.split("\n")) if m]
    con = duckdb.connect(os.path.join(work, "differential.duckdb"), config={"allow_unsigned_extensions": "true"})
    con.execute("LOAD '%s'" % extension)
    con.execute("CREATE SCHEMA IF NOT EXISTS core")
    region_id = target["region_id"]
    for table, names in columns.items():
        path = fetch(core + table + ".txt.gz", os.path.join(work, table + ".txt.gz"))
        spec = ", ".join("'%s': 'VARCHAR'" % name for name, _ in names)
        con.execute("CREATE OR REPLACE TABLE core.raw AS SELECT * FROM read_csv('%s', delim='\\t', header=false, "
                    "quote='', escape='', nullstr='\\N', columns={%s})" % (path, spec))
        select = ", ".join("TRY_CAST(%s AS BIGINT) AS %s" % (n, n) if t.lower().startswith(("int", "tinyint", "smallint", "mediumint", "bigint"))
                           else n for n, t in names)
        con.execute("CREATE OR REPLACE TABLE core.%s AS SELECT %s FROM core.raw" % (table, select))
    con.execute("DROP TABLE core.raw")
    con.execute("CREATE SCHEMA IF NOT EXISTS region")
    for table in ("attrib_type", "coord_system"):
        con.execute("CREATE OR REPLACE TABLE region.%s AS SELECT * FROM core.%s" % (table, table))
    replace = "* REPLACE (%d AS coord_system_id)" % target["coord_system_id"] if target["coord_system_id"] else "*"
    con.execute("CREATE OR REPLACE TABLE region.seq_region AS SELECT %s FROM core.seq_region WHERE seq_region_id = %d" % (replace, region_id))
    for table in ("seq_region_attrib", "gene", "transcript"):
        con.execute("CREATE OR REPLACE TABLE region.%s AS SELECT * FROM core.%s WHERE seq_region_id = %d" % (table, table, region_id))
    con.execute("CREATE OR REPLACE TABLE region.exon_transcript AS SELECT * FROM core.exon_transcript WHERE transcript_id IN (SELECT transcript_id FROM region.transcript)")
    con.execute("CREATE OR REPLACE TABLE region.exon AS SELECT * FROM core.exon WHERE exon_id IN (SELECT exon_id FROM region.exon_transcript)")
    con.execute("CREATE OR REPLACE TABLE region.translation AS SELECT * FROM core.translation WHERE transcript_id IN (SELECT transcript_id FROM region.transcript)")
    con.execute("CREATE OR REPLACE TABLE region.transcript_attrib AS SELECT * FROM core.transcript_attrib WHERE transcript_id IN (SELECT transcript_id FROM region.transcript)")
    con.execute("CREATE OR REPLACE TABLE region.translation_attrib AS SELECT * FROM core.translation_attrib WHERE translation_id IN (SELECT translation_id FROM region.translation)")
    # Reference: the region's record of the toplevel FASTA, one line per 60 bases, indexed.
    fasta = os.path.join(work, "region.fa")
    if not os.path.exists(fasta):
        gz = fetch(target["fasta"], os.path.join(work, "toplevel.fa.gz"))
        sequence, keep = [], False
        for line in gzip.open(gz, "rt"):
            if line.startswith(">"):
                keep = line.split()[0] == ">" + target["region"]
            elif keep:
                sequence.append(line.strip())
        sequence = "".join(sequence)
        with open(fasta, "w") as handle:
            handle.write(">%s\n" % target["region"])
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i:i + 60] + "\n")
        subprocess.run(["samtools", "faidx", fasta], check=True)
    sequence = "".join(l.strip() for l in open(fasta) if not l.startswith(">"))
    con.execute("CREATE OR REPLACE TABLE ref(chrom VARCHAR, \"start\" BIGINT, \"end\" BIGINT, seq VARCHAR)")
    con.execute("INSERT INTO ref VALUES (?, 0, ?, ?)", [target["region"], len(sequence), sequence])
    # Collection databases hold many species; the region's coordinate system names its species.
    species_id = con.execute("SELECT cs.species_id FROM region.seq_region sr JOIN region.coord_system cs USING (coord_system_id)").fetchone()[0]
    for name, function in (("regions", "duckvep_ensembl_regions_sql"), ("transcripts", "duckvep_ensembl_transcripts_sql")):
        con.execute("CREATE OR REPLACE TABLE %s AS FROM query(%s('region', 'ref', '%s', {species_id: %d}))"
                    % (name, function, target["assembly"], species_id))
    return con, sequence, fasta


def write_events(target, work, sequence, con):
    """SNVs at every base of windows around the origin and the crossing transcript's exon edges, sparse
    tiling elsewhere, and small indels."""
    length = len(sequence)
    rng = random.Random(63116)
    positions = set(list(range(1, 61)) + list(range(length - 59, length + 1)))
    for exon_start, exon_end in con.execute(
            "SELECT exon.exon_start, exon.exon_end FROM transcripts, unnest(exons) u(exon) "
            "WHERE transcript_stable_id = ?", [target["wrapped"]]).fetchall():
        for edge in (exon_start, exon_end):
            positions.update(range(max(1, edge - 120), min(length, edge + 120) + 1))
    positions.update(range(1, length + 1, max(37, length // 4000)))
    step = {"A": "C", "C": "G", "G": "T", "T": "A"}
    rows = []
    for p in sorted(positions):
        ref = sequence[p - 1]
        if ref not in step:
            continue
        rows.append((p, ref, step[ref]))
        if p % 3 == 0 and p + 3 <= length:
            rows.append((p, sequence[p - 1:p + 1 + p % 3 + (p // 3) % 2], sequence[p - 1]))
        if p % 5 == 0:
            rows.append((p, ref, ref + "".join(rng.choice("ACGT") for _ in range(1 + p % 4))))
    rows = sorted(set(rows))
    vcf = os.path.join(work, "events.vcf")
    with open(vcf, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n##contig=<ID=%s,length=%d>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     % (target["region"], length))
        for i, (p, r, a) in enumerate(rows):
            handle.write("%s\t%d\tev%d\t%s\t%s\t.\t.\t.\n" % (target["region"], p, i, r, a))
    return vcf, rows


def baseline_child(args):
    """Annotate the unwrapped transcripts with the pre-lift extension and write the rows to Parquet."""
    target = TARGETS[args.target]
    work = os.path.abspath(args.work)
    con = duckdb.connect(os.path.join(work, "baseline_input.duckdb"), config={"allow_unsigned_extensions": "true"})
    con.execute("LOAD '%s'" % args.extension)
    con.execute("""CREATE OR REPLACE TABLE transcripts_linear AS
      SELECT * REPLACE ((row_number() OVER (ORDER BY seq_region, transcript_start) - 1)::UINTEGER AS transcript_index)
      FROM transcripts WHERE NOT origin_crossing""")
    con.execute("""SELECT loaded FROM duckvep_model_load('baseline',
      'SELECT seq_region, sequence_length, seq_region_name, circular FROM regions',
      'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM transcripts_linear ORDER BY seq_region, transcript_start, transcript_index',
      'SELECT transcript_index, exon.exon_start, exon.exon_end, exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase FROM transcripts_linear, LATERAL unnest(exons) AS u(exon) ORDER BY transcript_index, exon.exon_cdna_start',
      peptide_edit_query := 'SELECT transcript_index, edit.protein_position, edit.alternate_amino_acid FROM transcripts_linear, LATERAL unnest(peptide_edits) AS u(edit) ORDER BY transcript_index, edit.protein_position',
      transcript_coverage_complete := true, reference_fasta := '%s')""" % os.path.join(work, "region.fa")).fetchall()
    con.execute("""COPY (SELECT a.event_index, t.transcript_stable_id AS tid, a.consequence_mask, a.transcript_hgvs, a.protein_hgvs
      FROM query(duckvep_annotate_sql('ev', 'baseline', struct_pack(hgvs := true, upstream_distance := 5000, downstream_distance := 5000))) a
      JOIN transcripts_linear t ON t.transcript_index = a.transcript_index) TO '%s' (FORMAT parquet)""" % os.path.join(work, "baseline.parquet"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", choices=sorted(TARGETS))
    ap.add_argument("work")
    ap.add_argument("extension")
    ap.add_argument("--baseline-extension", help="pre-lift extension: annotate the unwrapped transcripts with it and "
                    "compare, in a separate process because two builds cannot share a registry")
    ap.add_argument("--baseline-child", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.baseline_child:
        return baseline_child(args)
    target = TARGETS[args.target]
    wrapped = target["wrapped"]
    work = os.path.abspath(args.work)
    os.makedirs(work, exist_ok=True)
    extension = os.path.join(work, "duckvep.duckdb_extension")
    if os.path.exists(extension):
        os.chmod(extension, 0o644)
    shutil.copyfile(args.extension, extension)
    os.chmod(extension, 0o444)
    con, sequence, fasta = stage(target, work, extension)
    vcf, rows = write_events(target, work, sequence, con)
    length_of_sequence = len(sequence)
    cache_tar = fetch(target["cache"], os.path.join(work, "cache.tar.gz"))
    cache = os.path.join(work, "cache")
    if not os.path.isdir(cache):
        os.makedirs(cache)
        subprocess.run(["tar", "-xzf", cache_tar, "-C", cache], check=True)
    out = os.path.join(work, "vep")
    os.makedirs(out, exist_ok=True)
    os.chmod(out, 0o777)
    def run_vep(vcf_path, name, *extra):
        result = subprocess.run(["bash", os.path.join(HERE, "run_species_vep116_docker.sh"), target["species"],
                                 target["assembly"], "63", cache, fasta, vcf_path, os.path.join(out, name)] + list(extra),
                                capture_output=True, text=True)
        return result.returncode == 0

    # Consequences for every event. HGVS separately: VEP's transcript mapper can abort on the crossing
    # transcript (Mapper::map_insert), so HGVS runs on the events that miss its exons by more than 10 bases.
    if not run_vep(vcf, "vep.json"):
        raise SystemExit("VEP failed on the consequence run")
    exons = con.execute("SELECT exon.exon_start, exon.exon_end FROM transcripts, unnest(exons) u(exon) "
                        "WHERE transcript_stable_id = ?", [wrapped]).fetchall()
    exon_span = []
    for exon_start, exon_end in exons:
        exon_span.append((exon_start - 10, exon_end + 10) if exon_start <= exon_end
                         else (exon_start - 10, length_of_sequence + 10))
        if exon_start > exon_end:
            exon_span.append((0, exon_end + 10))
    hgvs_rows = [(i, p_, r_, a_) for i, (p_, r_, a_) in enumerate(rows)
                 if not any(lo <= p_ <= hi for lo, hi in exon_span)]
    hgvs_vcf = os.path.join(work, "events_hgvs.vcf")
    with open(hgvs_vcf, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n##contig=<ID=%s,length=%d>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     % (target["region"], length_of_sequence))
        for i, p_, r_, a_ in hgvs_rows:
            handle.write("%s\t%d\tev%d\t%s\t%s\t.\t.\t.\n" % (target["region"], p_, i, r_, a_))
    hgvs_ok = run_vep(hgvs_vcf, "vep_hgvs.json", "--hgvs")
    con.execute("""SELECT loaded FROM duckvep_model_load('target',
      'SELECT seq_region, sequence_length, seq_region_name, circular FROM regions',
      'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM transcripts ORDER BY seq_region, transcript_start, transcript_index',
      'SELECT transcript_index, exon.exon_start, exon.exon_end, exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase FROM transcripts, LATERAL unnest(exons) AS u(exon) ORDER BY transcript_index, exon.exon_cdna_start',
      peptide_edit_query := 'SELECT transcript_index, edit.protein_position, edit.alternate_amino_acid FROM transcripts, LATERAL unnest(peptide_edits) AS u(edit) ORDER BY transcript_index, edit.protein_position',
      transcript_coverage_complete := true, reference_fasta := '%s')""" % fasta).fetchall()
    con.execute("CREATE OR REPLACE TABLE ev_raw(id BIGINT, pos BIGINT, ref VARCHAR, alt VARCHAR)")
    con.executemany("INSERT INTO ev_raw VALUES (?,?,?,?)", [(i, p, r, a) for i, (p, r, a) in enumerate(rows)])
    con.execute("""CREATE OR REPLACE TABLE ev AS SELECT id::UBIGINT AS event_index, (SELECT seq_region FROM regions) AS seq_region,
      pos::UBIGINT AS position, ref AS reference, alt AS alternate, NULL::UBIGINT AS end_position, NULL::VARCHAR AS structural_type,
      NULL::VARCHAR AS copy_change, NULL::UINTEGER AS mate_seq_region, NULL::UBIGINT AS mate_position FROM ev_raw""")
    con.execute("""CREATE OR REPLACE TABLE duck AS
      SELECT e.event_index, e.position, t.transcript_stable_id AS tid, a.consequence_mask, a.transcript_hgvs, a.protein_hgvs,
        (SELECT string_agg(s.consequence, '&' ORDER BY s.consequence) FROM duckvep_so_terms() s WHERE (a.consequence_mask & s.consequence_mask) <> 0) AS terms
      FROM query(duckvep_annotate_sql('ev', 'target', struct_pack(hgvs := true, upstream_distance := 5000, downstream_distance := 5000))) a
      JOIN ev e USING (event_index) LEFT JOIN transcripts t ON t.transcript_index = a.transcript_index
      WHERE t.transcript_stable_id IS NOT NULL""")
    con.execute("""CREATE OR REPLACE TABLE vep AS
      SELECT replace(split_part(input, chr(9), 3), 'ev', '')::UBIGINT AS event_index, tc.transcript_id AS tid,
        (SELECT string_agg(x, '&' ORDER BY x) FROM unnest(tc.consequence_terms) u(x)) AS terms
      FROM read_json('%s', format='newline_delimited', maximum_object_size=67108864), unnest(transcript_consequences) u(tc)""" % os.path.join(out, "vep.json"))
    con.execute("CREATE OR REPLACE TABLE vep_hgvs(event_index UBIGINT, tid VARCHAR, hgvsc VARCHAR, hgvsp VARCHAR)")
    if hgvs_ok:
        con.execute("""INSERT INTO vep_hgvs
          SELECT replace(split_part(input, chr(9), 3), 'ev', '')::UBIGINT, tc.transcript_id,
            split_part(tc.hgvsc, ':', 2), split_part(tc.hgvsp, ':', 2)
          FROM read_json('%s', format='newline_delimited', maximum_object_size=67108864), unnest(transcript_consequences) u(tc)""" % os.path.join(out, "vep_hgvs.json"))
    con.execute("""CREATE OR REPLACE TABLE pairs AS
      SELECT coalesce(v.event_index, d.event_index) AS event_index, coalesce(v.tid, d.tid) AS tid, d.position AS position,
        v.terms AS vep_terms, d.terms AS duck_terms, h.hgvsc, d.transcript_hgvs AS duck_hgvsc, h.hgvsp, d.protein_hgvs AS duck_hgvsp,
        h.event_index IS NOT NULL AS has_hgvs_run
      FROM vep v FULL OUTER JOIN duck d ON v.event_index = d.event_index AND v.tid = d.tid
      LEFT JOIN vep_hgvs h ON h.event_index = coalesce(v.event_index, d.event_index) AND h.tid = coalesce(v.tid, d.tid)""")
    length = len(sequence)
    receipt = {
        "target": args.target, "species": target["species"], "region": target["region"], "region_length": length,
        "extension_sha256": sha256(extension),
        "events": len(rows), "events_in_hgvs_run": len(hgvs_rows), "vep_hgvs_run_succeeded": hgvs_ok,
        "region_transcripts": con.execute("SELECT count(*) FROM transcripts").fetchone()[0],
        "origin_crossing_transcripts": con.execute("SELECT string_agg(transcript_stable_id || ':' || strand::VARCHAR, ',') FROM transcripts WHERE origin_crossing").fetchone()[0],
        "fasta_sha256": sha256(fasta), "cache_sha256": sha256(cache_tar),
        "transcript_dump_sha256": sha256(os.path.join(work, "transcript.txt.gz")),
    }
    other = "tid <> '%s'" % wrapped
    receipt["unwrapped"] = dict(zip(
        ("pairs", "identical_terms", "duckvep_only", "vep_only", "different_terms", "hgvsc_pairs", "hgvsc_identical", "hgvsp_pairs", "hgvsp_identical"),
        con.execute("""SELECT count(*), count(*) FILTER (WHERE vep_terms = duck_terms), count(*) FILTER (WHERE vep_terms IS NULL),
          count(*) FILTER (WHERE duck_terms IS NULL), count(*) FILTER (WHERE vep_terms <> duck_terms),
          count(*) FILTER (WHERE has_hgvs_run AND vep_terms = duck_terms AND hgvsc IS NOT NULL), count(*) FILTER (WHERE has_hgvs_run AND vep_terms = duck_terms AND hgvsc = duck_hgvsc),
          count(*) FILTER (WHERE has_hgvs_run AND vep_terms = duck_terms AND hgvsp IS NOT NULL), count(*) FILTER (WHERE has_hgvs_run AND vep_terms = duck_terms AND hgvsp = duck_hgvsp)
          FROM pairs WHERE %s""" % other).fetchone()))
    # Every DuckVEP-only row of an unwrapped transcript must be a flank row that exists only through the origin.
    receipt["unwrapped"]["duckvep_only_across_origin"] = con.execute(f"""
      WITH d AS (SELECT p.position::BIGINT AS a, (r.pos + length(r.ref) - 1)::BIGINT AS b, t.transcript_start::BIGINT AS s,
          t.transcript_end::BIGINT AS e, p.duck_terms
        FROM pairs p JOIN transcripts t ON t.transcript_stable_id = p.tid JOIN ev_raw r ON r.id = p.event_index
        WHERE p.vep_terms IS NULL AND {other})
      SELECT count(*) FILTER (WHERE greatest(0, s - b, a - e) > 5000
          AND least(greatest(0, s + {length} - b, a - e - {length}), greatest(0, s - {length} - b, a - e + {length})) <= 5000
          AND duck_terms IN ('upstream_gene_variant', 'downstream_gene_variant')) FROM d""").fetchone()[0]
    # HGVS 3' shifting looks up to 1000 bases each way. VEP clips that window at the sequence ends; the
    # circular model does not, so HGVS differences are expected within 1100 bases of the origin only.
    receipt["unwrapped"]["hgvs_differences_near_origin_vs_interior"] = dict(zip(
        ("hgvsc_near_origin", "hgvsc_interior", "hgvsp_near_origin", "hgvsp_interior"),
        con.execute(f"""SELECT count(*) FILTER (WHERE hgvsc IS DISTINCT FROM duck_hgvsc AND (position <= 1100 OR position > {length} - 1100)),
          count(*) FILTER (WHERE hgvsc IS DISTINCT FROM duck_hgvsc AND position > 1100 AND position <= {length} - 1100),
          count(*) FILTER (WHERE hgvsp IS DISTINCT FROM duck_hgvsp AND (position <= 1100 OR position > {length} - 1100)),
          count(*) FILTER (WHERE hgvsp IS DISTINCT FROM duck_hgvsp AND position > 1100 AND position <= {length} - 1100)
          FROM pairs WHERE {other} AND has_hgvs_run AND vep_terms = duck_terms""").fetchone()))
    receipt["wrapped_transcript"] = dict(zip(
        ("pairs", "identical_terms", "duckvep_only", "vep_only", "different_terms"),
        con.execute("""SELECT count(*), count(*) FILTER (WHERE vep_terms = duck_terms), count(*) FILTER (WHERE vep_terms IS NULL),
          count(*) FILTER (WHERE duck_terms IS NULL), count(*) FILTER (WHERE vep_terms <> duck_terms)
          FROM pairs WHERE tid = ?""", [wrapped]).fetchone()))
    receipt["wrapped_transcript"]["vep_terms_by_duckvep_terms"] = [
        list(r) for r in con.execute("""SELECT coalesce(vep_terms, '(none)'), coalesce(duck_terms, '(none)'), count(*) FROM pairs
          WHERE tid = ? AND vep_terms IS DISTINCT FROM duck_terms GROUP BY ALL ORDER BY 3 DESC LIMIT 6""", [wrapped]).fetchall()]
    if args.baseline_extension:
        con.execute("CHECKPOINT")
        shutil.copyfile(os.path.join(work, "differential.duckdb"), os.path.join(work, "baseline_input.duckdb"))
        subprocess.run([sys.executable, os.path.abspath(__file__), "--baseline-child",
                        args.target, work, args.baseline_extension], check=True)
        base = os.path.join(work, "baseline.parquet")
        receipt["baseline_extension_sha256"] = sha256(args.baseline_extension)
        # The pre-lift extension sees the sequence as linear. Rows may differ from the lifted result only where the
        # circular topology reaches: flanks across the origin and HGVS windows within 1100 bases of it.
        receipt["lifted_vs_prelift_linear"] = dict(zip(
            ("pairs", "identical", "different_interior", "different_near_origin", "lifted_only_rows_near_origin",
             "prelift_only_rows"),
            con.execute(f"""WITH j AS (SELECT coalesce(b.event_index, d.event_index) AS event_index, coalesce(b.tid, d.tid) AS tid,
                d.position, b.consequence_mask AS bm, d.consequence_mask AS dm,
                b.transcript_hgvs AS bh, d.transcript_hgvs AS dh, b.protein_hgvs AS bp, d.protein_hgvs AS dp
              FROM read_parquet('{base}') b FULL OUTER JOIN duck d ON b.event_index = d.event_index AND b.tid = d.tid
              WHERE coalesce(b.tid, d.tid) <> '{wrapped}'),
              p AS (SELECT j.*, e.pos AS position2 FROM j JOIN ev_raw e ON e.id = j.event_index)
              SELECT count(*), count(*) FILTER (WHERE bm = dm AND bh IS NOT DISTINCT FROM dh AND bp IS NOT DISTINCT FROM dp),
                count(*) FILTER (WHERE bm IS NOT NULL AND dm IS NOT NULL AND (bm <> dm OR bh IS DISTINCT FROM dh OR bp IS DISTINCT FROM dp)
                  AND position2 > 1100 AND position2 <= {length} - 1100),
                count(*) FILTER (WHERE bm IS NOT NULL AND dm IS NOT NULL AND (bm <> dm OR bh IS DISTINCT FROM dh OR bp IS DISTINCT FROM dp)
                  AND (position2 <= 1100 OR position2 > {length} - 1100)),
                count(*) FILTER (WHERE bm IS NULL AND (position2 <= 5100 OR position2 > {length} - 5100)),
                count(*) FILTER (WHERE dm IS NULL)
              FROM p""").fetchone()))
    print(json.dumps(receipt, indent=1))


if __name__ == "__main__":
    main()
