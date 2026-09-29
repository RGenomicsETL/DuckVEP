#!/usr/bin/env bash
# Rebuild every scale panel from the staged gnomAD v4.1 lean Parquet, in one command:
#
#   scripts/gnomad_build_panels.sh
#
# 1M / 5M / 25M / 100M distinct genome panels, the 2M coding/splice exome panel and the 12,000-row
# structural controls, all under $DUCKVEP_GNOMAD_ROOT/panels (never committed), then the committed
# receipts in benchmarks/data/scale_contracts/panels/ and the links the runner resolves.
# Run it again after more genome shards land: the source set changes, so new panels and receipts
# are produced; a finished artifact for an unchanged source set is reused. Directories named
# *.partial are still staging and are never read. The 100M panel stays "unavailable"
# (insufficient_distinct_alleles) until enough genome shards have landed.
#
# Environment: DUCKVEP_GNOMAD_ROOT, DUCKVEP_SCALE_MODEL, DUCKVEP_PANEL_THREADS (4),
# DUCKVEP_PANEL_MEMORY (16GB), DUCKVEP_PANEL_TEMP (12GiB, capped so 50 GiB stay free).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export DUCKVEP_PANEL_THREADS="${DUCKVEP_PANEL_THREADS:-4}"
export DUCKVEP_PANEL_MEMORY="${DUCKVEP_PANEL_MEMORY:-16GB}"
export DUCKVEP_PANEL_TEMP="${DUCKVEP_PANEL_TEMP:-12GiB}"
last() { tail -n 1 | tr -d '[:space:]'; }
genomes="$(Rscript "$here/gnomad_genome_panels.R" | last)"
exomes="$(Rscript "$here/gnomad_exome_panel.R" | last)"
structural="$(Rscript "$here/gnomad_structural_controls.R" | last)"
Rscript "$here/gnomad_panel_receipts.R" "$genomes" "$exomes" "$structural"
