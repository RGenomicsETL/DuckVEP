# Reached-effect equality audit

`test/sql_v2/equality_golden.json` is a cross-implementation baseline, not an effect oracle. Of its 636 cases, 628 are unchanged and eight differ in the qualified reached-effect candidate. All 636 current v1/v2 outcomes agree, including 198 expected errors. The JSON baseline is not rewritten by this audit.

Candidate binary SHA256 values:

- v1: `d075b1503ccab5a1a320103d57acd26710fec0185b9fa804e0c27785a0d0f311`
- v2: `904f3b6f7be77efad68b639c92d5c47235e62d8dd5ecb34cca3517662b697339`

The old and new implementations were also run in separate processes against identical fixtures to record all 1,010 complete `hap_sc`, `hap_ss`, `hap_nm` and `hap_nm_hgvs` rows. Exactly 80 rows differ, only in top-level and per-carrier consequence sets. Raw operands, contributor identity, HGVS, impacts, NMD and every other field compare exactly. Whole-row hashes below are derived checksums, not biological proofs.

| Equality case | Baseline → candidate | Independent evidence |
| --- | --- | --- |
| `haplotypes sc nested carriers` | Six `far_syn_and_missense` tuples gain `synonymous_variant` alongside missense. | Final codon2 `GCT`→`GGT` changes A→G; codon3 `GAA`→`GAG` retains E; both precede unchanged stop. Pinned csq agrees. |
| `haplotypes sc full-row hash` | `(n,total,x,transcripts)` `(144,1278663369435005428560,4495889223471160240,138)` → `(144,1266892723490166823604,12412779472436768854,138)`. | All144 complete rows compared; six consequence-only differences. |
| `haplotypes ss nested carriers` | Six stop-lost lanes gain upstream missense; six stop-retained lanes gain synonymous; six protein-altering lanes gain in-frame insertion. | Codon17 `GTC`→`GAC` precedes terminal `TAA`→`CAA`; codon3 `CCC`→`CCA` precedes terminal `TAA`→`TGA`; a source-anchored `GGG` insertion after base60 supplies G before terminal `TAG`. Pinned csq agrees on the first two pairs; the declared mixed-class divergence for the insertion remains. |
| `haplotypes ss full-row hash` | `(366,3616752480536518701126,8127280110378016880,348)` → `(366,3657813613620810045685,2650681600858943641,348)`. | All366 complete rows compared;18 consequence-only differences. |
| `haplotypes nm nested carriers` |28 lanes gain an upstream in-frame insertion(16) or deletion(12) alongside stop-gained. | `check_nmd_reached_effects.R`: complete `GGG` codons inserted after base30 or complete `CAA` codons deleted after base12/30 precede first translated `TAA`. The independent genome-derived NMD gate retains S/J and class; csq is not an NMD oracle. |
| `haplotypes nm full-row hash` | `(250,2206918717975886787721,9268254449696724081,238)` → `(250,2214137170784349814859,3113897188211712037,238)`. | All250 complete rows compared;28 consequence-only differences. |
| `haplotypes nm hgvs full-row hash` | `(n,total,x)` `(250,2203978465762580879390,5136244193682519752)` → `(250,2195180746304829467013,17558476449629570859)`. | All250 complete HGVS-bearing rows compared;28 consequence-only differences, with exact unchanged HGVS fields. |
| `duplicate deletion immediate stop full row` | `stop_gained` → `stop_gained,inframe_deletion`; `M*` unchanged. | Distinct source assertions delete original codon2 `GCT` and change original codon3 `GAA`→`TAA`. The retained-coordinate map `1..3,7..15` places the reached stop at original bases7..9, after the deleted codon. `effects_ancestry_witness.R` verifies aliases, unique-source and reordered forms on both strands. A single opaque `GCTG>T` assertion makes the same RNA without proving that codon deletion and reports only stop-gained. Deletion-only, stop-only and post-stop-deletion controls distinguish these claims. Every other full-row field is unchanged. |

The twelve frame-restored identical-peptide lanes retain aggregate synonymous/LOW: equality of the complete determinate peptide is a whole-product fact, not local displaced-codon ancestry. Pinned csq's additional protein-altering term remains a declared divergence. The lowercase-overlap case is unchanged; case-only differences do not establish a changed codon or retained start.

Local receipts under `work/remaining-implementation/evidence/`: `effects-corrected-host-compare.log`, `effects-corrected-host-deltas.rds`, `effects-full-rows-{baseline,candidate}.json`, `audit-effects-full-rows.R`, `effects-full-rows-audit.log` and `effects-immediate-stop-ancestry.log`. These qualify the named binaries, not later source additions or the full original #11/#50 scope.
