#!/bin/bash
# Gate 1: throughput floor on the frozen candidate.
cd /root/duckvep/wt-closure
OUT=/root/duckvep/data/closure/g1
EXT=/root/duckvep/data/closure/ext-660a4ed/duckvep.duckdb_extension
MODEL=/root/duckvep/data/scale5/model.duckdb
CORPUS=/root/duckvep/data/scale5/corpus.duckdb
SRC=/root/.cache/duckhts/benchmarks/variantkey-providers/raw/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
FA=/root/duckvep/data/reference/ensembl-116/Homo_sapiens.GRCh38.dna.primary_assembly.fa
waitload() { while :; do l=$(cut -d' ' -f1 /proc/loadavg); if awk -v l=$l 'BEGIN{exit !(l<=3)}'; then break; fi; sleep 60; done; echo "$1 load1=$(cut -d' ' -f1-3 /proc/loadavg)" >> $OUT/load.log; }
run() { # mode threads cpus run
  waitload "$1 t$2 run$4"
  extra=""; case $1 in hgvs|rich_hgvs) extra="--reference-fasta $FA";; esac
  taskset -c $3 Rscript benchmarks/duckvep_throughput.R --extension $EXT --skip-extension-build \
    --database $MODEL --variants-database $CORPUS --variants-table bench_variants --corpus-source $SRC \
    --workload-name ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation --api-surface public \
    --regulatory --variants 4095611 --passes 5 --warmup 100000 --threads $2 --input-partitions 1 \
    --transcript-distance 5000 --output $1 $extra --history $OUT/history.csv \
    --fingerprint $OUT/fp-$1-t$2-r$4.csv > $OUT/$1-t$2-r$4.log 2>&1
  echo "$1 t$2 r$4 exit=$?" >> $OUT/load.log
}
for r in 1 2 3 4 5; do run compact 1 2 $r; done
for r in 1 2 3; do run rich 1 2 $r; run hgvs 1 2 $r; run rich_hgvs 1 2 $r; done
for r in 1 2 3 4 5; do run compact 4 2,4,6,8 $r; done
for r in 1 2 3; do run rich 4 2,4,6,8 $r; run hgvs 4 2,4,6,8 $r; run rich_hgvs 4 2,4,6,8 $r; done
echo DONE >> $OUT/load.log
