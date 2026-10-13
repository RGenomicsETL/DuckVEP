# Conditional-programme and reached-effect equality audit

The 636-case reference contains 30 individually adjudicated changes: eight reached-effect outcomes, thirteen Ensembl transcript-builder outcomes, and nine receipt-builder outcomes. The other 606 cases, including all 198 expected errors, are unchanged. This is not permission to refresh unrelated outcomes.

## Reached effects

The eight outcomes and original-allele/codon reasoning are documented in [the effect audit](../data/haplotype/effects_equality_audit.md). Full-row comparisons distinguish added SO terms from raw CDS/protein/HGVS, source/carrier identity, impact and NMD. Whole-product equality does not establish local opaque-ALT ancestry. Unchanged ambiguous codons do not create changed-unknown effects; changed or untraceable ambiguity remains covered by the registered stream properties.

## Ensembl transcript builder

The thirteen `builder ensembl transcripts` cases cover defaults, qualification/quoting, NULL options, species bounds and the multirow batch. Generated SQL includes original RNA coding start/end as UBIGINT and code as UTINYINT. Missing authority stays NULL; compatibility CDS/protein/code fields are not repurposed.

Full fixture relation comparisons retain every pre-existing field except a justified source flag on transcript 32: flags 259→131331 preserves its literal `initial_met` source assertion. That assertion remains a conditional-RNA veto. Populated triples on transcripts 31/32 in the small fixture are respectively 1..9/code1 and 1..9/code2. Species-isolated control rows are exact. `duckvep_ensembl.test`, `duckvep_rna_programme.test` and the paired RNA producer/replay scripts exercise the actual relations and predictions, not only builder text.

## Prepared-model receipt

Nonempty named RNA authority participates in content identity. Older 36-column relations and all-NULL triples retain their complete legacy receipt. Receipt identity commits to supplied content, not biological validity.

The eight simple receipt-builder success cases and the 1,500-row aggregate differ only in the row-payload/hash block. Complete generated SQL comparison covers all 1,508 statements, including quoting and optional regulation forms. Source columns named `model`, `authority` or internal payload names cannot shadow the named projection; column reordering and unrelated-column controls are in `receipt_integrity.sql`. All six receipt error outcomes remain exact.

On the pinned 39-transcript GRCh37/38 relations, only indexes 15, 18, 20, 21, 26, 31, 32 contribute RNA authority. Independent reconstruction preserves the 34-member legacy serialization, conditionally appends the three RNA members, hashes rows and 4096-row blocks, and combines the unchanged provenance, region and regulation hashes. All 25 other receipt fields and all input rows are exact. Removing or NULLing authority reproduces each complete legacy receipt.

| Assembly | Content SHA-256 with RNA authority |
|---|---|
| GRCh37 | `cc6c7fe9de77e61c52011b5c3b68e0017a1112f2a26e63954bde2e7ae70a1b1f` |
| GRCh38 | `ac35b613089f6ba2f9d15ffa2a91c4625cee24d2fc9496c9dd08b5ebb4217d5f` |

The two selected assertions in `duckvep_ensembl.test` use these independently reconstructed hashes. `run_receipt_integrity.sh` checks both hosts; `run_receipt_legacy.sh` compares complete older-schema receipts against an explicitly supplied preserved extension.

## Qualification boundary

The final integrated v1/v2 run agrees on all 636 outcomes. Selecting only the 30 reviewed entries produces a source reference byte-identical to that live-qualified recording, SHA-256 `6914ca324648fdeaadc80164983761dc2698b33a4132852f214c4803ded40214`; the preceding 636-case reference is `d773bcd44fe4a852e8d7cf65639c44bb7d5ead5c8124f60d0544387755eb99af`.

This audit does not establish wider runtime/platform support, naturally expressed hypothetical products, general NMD, observed inter-record phase, or full #11/#50 closure.
