# DuckVEP preprint — skeleton (DRAFT, not submitted)

Status: **outline only.** No manuscript has been written, no DOI exists, and the
data/model archive provider is still undecided (issue #1). Nothing here is a
publication claim.

Target: bioRxiv preprint, then a journal submission (venue undecided).

## Metadata (TODO)

- **Title (working):** DuckVEP: relation-first variant consequence and phased
  haplotype annotation inside DuckDB
- **Authors / affiliations / ORCIDs:** TODO — owner to confirm the author list
- **Corresponding author:** Sounkou Mahamane Toure <sounkoutoure@gmail.com>
- **Preprint DOI:** `10.1101/XXXXXXX` — TODO
- **Zenodo software DOI:** `10.5281/zenodo.XXXXXXX` — TODO (record not minted)
- **Licence:** code GPL-2.0-or-later; manuscript TODO

## Abstract (TODO)

Draft after the Results numbers are frozen. Must state the bounded scope
explicitly: independent-event annotation conformant to Ensembl VEP 116 in its
declared domain, plus a bounded whole-haplotype coding profile; not a claim of
whole-haplotype coverage for every VEP consequence class.

## 1. Introduction

- The cost of moving variants to an annotation service vs. annotating where the
  data already lives (DuckDB relations).
- Consequence engines as oracles: Ensembl VEP, bcftools csq, Haplosaurus,
  Mutalyzer.
- Contribution: execute the oracle semantics in-relation, with an executable
  conformance contract rather than a prose claim of equivalence.

## 2. Implementation

- DuckDB C extension (`duckvep`) + R front end (`Rduckvep`).
- Two contracts, kept separate:
  - **independent-event** per-record annotation (VEP 116 authority);
  - **whole-haplotype** prediction over phased coding records (bounded profile).
- Storage/compute model: relation-first, no per-variant service round trip.
- Model/ABI: transcript model, source-revision pinning, receipt checks.
- Sources: `design/duckvep.md`, `design/duckvep_haplotype_contract.md`,
  `docs/v2-host.md`.

## 3. Methods — executable conformance

- Independent-event oracle: Ensembl VEP 116 recomputed locally; compared field
  by field in the declared comparable domain.
  - Evidence: `test/duckvep/conformance/`.
- Haplotype oracle: bcftools csq for the comparable coding domain (coarse term
  vocabulary, identical transcript geometry); Haplosaurus for whole-haplotype
  consequences; Mutalyzer for HGVS.
- Disagreeing cases are enumerated as a refusal/errata ledger, not silently
  dropped: `ERRATA.md`, `test/duckvep/conformance/README.md`.
- Statistical/property testing: `test/duckvep/property/`.
- Reproducibility: pinned revisions + receipts; results regenerate from tracked
  scripts only.

## 4. Results (TODO — freeze numbers, keep evidence links)

Each number must cite a tracked artifact. Do not restate a figure from an older
revision as a claim for the current one.

| Claim | Evidence artifact | Status |
|---|---|---|
| Throughput at gnomAD scale (≥1M alleles/s) | `benchmarks/duckvep_throughput*.md`, `benchmarks/figures/throughput-*.png` | delivered (#3) |
| Memory bound enforced per job | `benchmarks/figures/throughput-dense-memory-1.png` | delivered (#3) |
| bcftools csq agreement on the comparable coding domain | `benchmarks/benchmark_duckvep_vep_rs.md` | delivered |
| Independent-event VEP 116 agreement | `test/duckvep/conformance/` | delivered |
| Phased haplotype classification timing/workspace | `benchmarks/figures/phase-*.png` | delivered |
| Projection cost | `benchmarks/figures/projection-*.png` | delivered |
| HGVS (DNA + protein) vs Mutalyzer/VEP | `benchmarks/duckvep_conformance.*`, conformance README | delivered |

Headline results from the earlier qualification retain their original scope
labels; they are not new performance claims for later revisions.

## 5. Discussion and limitations

- Scope boundary: which consequence classes the whole-haplotype profile covers,
  and which remain independent-event only.
- `compound_hgvs_certified` stays **FALSE**; state exactly what is and is not
  certified.
- Overlap/refusal classes and why they are refused (issue #11 boundary).
- Open qualification: released DuckDB 2.0.0 (issue #48).

## 6. Availability

- **Code:** https://github.com/RGenomicsETL/DuckVEP (tag TODO);
  R package `Rduckvep` on CRAN (TODO — not yet submitted);
  DuckDB community extension (TODO — submission withdrawn pending CRAN).
- **Model/data archive:** TODO — provider undecided (owner Google Drive vs
  Zenodo). A versioned Zenodo record with a manifest is the plan
  (`design/duckvep.md`).
- **Errata:** `ERRATA.md`.

## 7. References (TODO)

HTSlib (10.1093/gigascience/giab007), VariantKey/RegionKey (10.1101/473744),
Ensembl VEP, bcftools, Haplosaurus, Mutalyzer, DuckDB.

## Before posting — checklist

- [ ] All Results numbers re-run on a frozen revision and linked to artifacts
- [ ] Software DOI minted (Zenodo) and added to `CITATION.cff`, `inst/CITATION`
- [ ] CRAN + community-extension availability state matches the text
- [ ] Figures copied from `benchmarks/figures/` at final resolution
- [ ] Limitations section reviewed against `ERRATA.md`
- [ ] Owner sign-off on authorship and scope claims
