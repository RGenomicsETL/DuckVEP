# Release checklist

This checklist applies to the DuckVEP extension and the Rduckvep package.
Stable-v1 is the default host, also on DuckDB 2.0.
Native-v2 is experimental and uses the version-locked DuckDB 2.0.0 SDK.
A stable-v1 release does not require a change to the default host.

## Owner decisions

- Choose the extension version and the R package version.
- Choose the biological scope for the release.
- Decide whether to archive the software on Zenodo.
- Approve the release candidate commit before any tag or submission.

## Candidate and source checks

- Complete the fixes and regression tests for the accepted release scope.
- Record each unresolved discrepancy and its effect on public output.
- Synchronize the R bundle with the extension sources.
- Generate the R documentation and the authored Markdown reports.
- Run the SQL, R, function-documentation, sanitizer, property, and host compatibility checks.
- Record the candidate commit and artifact checksums with each result.
- Confirm that all distribution jobs complete successfully on the candidate.
- Distinguish a runner acquisition failure from a build or test failure.

## Biological evidence

- Run all nine conformance corpora with executable VEP 116 on the candidate.
- Run the HGVS corpora with the same candidate.
- Record the exact VEP revisions, reference assembly, model identity, options, and fixture checksums.
- Compare consequence sets, HGVSc, HGVSp, and support statuses separately.
- Keep VEP 114-to-116 behavior changes outside DuckVEP defect expectations.
- Identify each reconstructed-model result separately from a full-cache result.
- Explain each remaining difference before a parity claim.
- Generate the README from the committed receipts.

Commit ancestry alone does not establish conformance for a candidate.
A comparison that omits allele identity does not establish allele-level parity.
Historical receipts remain evidence for their recorded revisions.

## Runtime and performance evidence

- Run the whole-HG002 comparison on both hosts with the released DuckDB 2.0.0 runtime.
- Record the runtime, SDK, extension, model, reference, and input identities.
- Confirm equal output rows and exact fingerprints where the comparison requires equality.
- Record model budgets, annotation budgets, elapsed time, and throughput separately.
- Refresh the throughput and vep-rs comparisons with their declared inputs and methods.
- Use the [native-v2 preflight](../test/sql_v2/release-readiness.md) before native-v2 admission.

Native-v2 preflight does not replace biological, platform, or distribution checks.
Linux OOM recovery evidence does not establish recovery on macOS, Windows, or WebAssembly.

## R package and CRAN

- Build the source tarball from the candidate.
- Run `R CMD check --as-cran` with the supported R and DuckDB versions.
- Run win-builder checks for R release and R devel.
- Record all errors, warnings, and notes in `r/Rduckvep/cran-comments.md`.
- Distinguish a check without remote incoming checks from a complete incoming check.
- Confirm that examples and tests use local fixtures without network access.
- Submit to CRAN only after owner approval.

## Metadata and distribution

- Make the manifest, NEWS, package metadata, and citation agree with the selected versions.
- Use only an actual DOI in `CITATION.cff` if the owner chooses Zenodo.
- Record immutable model repository revisions and checksums for public model files.
- Make the model download instructions agree with the chosen public host.
- Update the umbrella tracker and the native-v2 tracker with verified evidence.
- Prepare a community-extensions submission at the approved candidate commit.
- Coordinate community distribution with the declared CRAN release policy.
- Publish the tag, release assets, and submissions only with owner approval.

The model repository is [RgenomicsETL/duckvep-models](https://huggingface.co/datasets/RgenomicsETL/duckvep-models).
[Publishing to Hugging Face](publishing-hugging-face.md) gives the staging, review and upload procedure.
A download URL alone does not identify a model revision.
The closed community submission [#2970](https://github.com/duckdb/community-extensions/pull/2970) is not a current release submission.

## Release record

For each gate, retain the command, source commit, artifact checksums, runtime versions, result, and evidence location.
Mark a gate complete only when its evidence applies to the candidate.
Record any explicit owner-approved scope exclusion beside the affected gate.
