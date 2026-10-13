library(tinytest)
library(DBI)

local({
  artifact <- Sys.getenv("DUCKVEP_PRODUCTS_EXTENSION", "")
  con <- rduckvep_connect(extension_path = if (nzchar(artifact)) artifact else NULL)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  run <- getFromNamespace("rduckvep_product_scenario_v1", "Rduckvep")
  reference <- list(start1 = 1L,
                    dna = "ATGAAACCCGGGAAAGGGTTTCCCTAAGGGATGAAATAAGGGAAA", strand = 1L)
  exons <- data.frame(start1 = c(1L, 13L, 19L, 31L, 43L),
                      end1 = c(9L, 15L, 27L, 39L, 45L))
  sources <- data.frame(source_id = c(11L, 22L), start1 = c(13L, 29L),
                        reference = c("A", "G"), alternate = c("C", "C"))
  programme <- list(id = "include", evidence_id = "declared", assumption_id = "splice",
                    expected_reference_rna = "ATGAAACCCAAATTTCCCTAAATGAAATAAAAA",
                    annotated_start1 = 1L, alternative_start1 = 0L,
                    alternative_trigger_id = "", genetic_code = 1L)
  original <- list(reference, exons, sources, programme)
  included <- run(con, reference, exons, exons, sources, programme)
  expect_identical(included$alternate_rna, "ATGAAACCCCAATTTCCCTAAATGAAATAAAAA")
  expect_identical(included$allele_evidence, "hypothetical_assignment")
  expect_identical(included$product_evidence, "supplied_programme")
  expect_identical(included$conditional_nmd, "unknown")
  expect_equal(included$verified_sources, 2L)
  expect_equal(included$sources[[1L]],
               transform(sources, source_id = as.character(source_id),
                         start1 = as.character(start1)))
  expect_equal(included$source_effects[[1L]]$source_id, as.character(sources$source_id))
  expect_equal(included$source_effects[[1L]]$changed_rna_bases, c(1L, 0L))
  expect_equal(included$source_effects[[1L]]$genomic_start1, as.character(sources$start1))
  skipped <- run(con, reference, exons, exons[c(1L, 3L:5L), ], sources,
                 modifyList(programme, list(id = "skip",
                   expected_reference_rna = "ATGAAACCCTTTCCCTAAATGAAATAAAAA")))
  expect_equal(skipped$source_effects[[1L]]$changed_rna_bases, c(0L, 0L))
  retained <- run(con, reference, exons,
                  data.frame(start1 = c(1L, 10L, 13L, 19L, 31L, 43L),
                             end1 = c(9L, 12L, 15L, 27L, 39L, 45L)),
                  sources, modifyList(programme, list(id = "retain",
                    expected_reference_rna = "ATGAAACCCGGGAAATTTCCCTAAATGAAATAAAAA")))
  expect_equal(retained$source_effects[[1L]]$changed_rna_bases, c(1L, 0L))
  swapped <- run(con, reference, exons, exons, sources[2L:1L, ], programme)
  expect_equal(swapped$source_effects[[1L]]$source_id, c("22", "11"))
  zero_sources <- sources
  zero_sources$source_id[1L] <- 0L
  zero <- run(con, reference, exons, exons, zero_sources, programme)
  expect_identical(zero$sources[[1L]]$source_id, c("0", "22"))
  expect_identical(zero$source_effects[[1L]]$source_id, c("0", "22"))
  expect_error(run(con, reference, exons, exons,
                   transform(zero_sources, source_id = c(0L, 0L)), programme),
               "ambiguous composition")
  shifted_reference <- modifyList(reference, list(start1 = 4294967297))
  shifted_exons <- transform(exons, start1 = start1 + 4294967296,
                            end1 = end1 + 4294967296)
  shifted_sources <- transform(sources, start1 = start1 + 4294967296)
  shifted <- run(con, shifted_reference, shifted_exons, shifted_exons,
                 shifted_sources, programme)
  expect_equal(shifted$source_effects[[1L]]$genomic_start1,
               c("4294967309", "4294967325"))
  termination <- function(dna, position, ref, alt, alternative_start = 0L) {
    reference <- list(start1 = 1L, dna = dna, strand = 1L)
    interval <- data.frame(start1 = 1L, end1 = nchar(dna))
    source <- data.frame(source_id = 0L, start1 = position, reference = ref, alternate = alt)
    declared <- programme
    declared$expected_reference_rna <- dna
    declared$alternative_start1 <- alternative_start
    declared$alternative_trigger_id <- if (alternative_start > 0L) "annotated-start-lost" else ""
    run(con, reference, interval, interval, source, declared)
  }
  post_stop <- termination("ATGAAATAATTT", 10L, "T", "C")
  expect_identical(post_stop$reference_peptide, "MK*")
  expect_identical(post_stop$alternate_peptide, "MK*")
  expect_identical(post_stop$alternate_rna, "ATGAAATAACTT")
  stop_lost <- termination("ATGAAATAATTTTGA", 7L, "T", "C")
  expect_identical(stop_lost$reference_peptide, "MK*")
  expect_identical(stop_lost$alternate_peptide, "MKQF*")
  no_stop <- termination("ATGAAATAATTT", 7L, "T", "C")
  expect_identical(no_stop$alternate_peptide, "MKQF")
  declared_ctg <- termination("CTGAAATAA", 4L, "A", "G")
  expect_identical(declared_ctg$reference_peptide, "MK*")
  expect_identical(declared_ctg$alternate_peptide, "ME*")
  alternative_ctg <- termination("ATGAAACTGAAATAA", 2L, "T", "C", 7L)
  expect_identical(alternative_ctg$reference_peptide, "MKLK*")
  expect_identical(alternative_ctg$alternate_peptide, "MK*")
  expect_equal(alternative_ctg$effective_start1, 7L, tolerance = 0)
  expect_identical(alternative_ctg$conditional_nmd, "unknown")
  decimal_sources <- sources
  decimal_sources$source_id <- c(1e15, 2^53)
  decimal <- run(con, reference, exons, exons, decimal_sources, programme)
  expect_identical(decimal$sources[[1L]]$source_id,
    c("1000000000000000", "9007199254740992"))
  expect_identical(decimal$source_effects[[1L]]$source_id,
    c("1000000000000000", "9007199254740992"))
  expect_error(run(con, reference, exons, exons,
    transform(sources, source_id = c(2^53 + 2, 1)), programme), "exact nonnegative")
  top <- run(con, list(start1 = "18446744073709551607", dna = "ATGAAATAA", strand = 1L),
             data.frame(start1 = "18446744073709551607",
                        end1 = "18446744073709551615"),
             data.frame(start1 = "18446744073709551607",
                        end1 = "18446744073709551615"),
             data.frame(source_id = "11", start1 = "18446744073709551610",
                        reference = "A", alternate = "C"),
             list(id = "edge", evidence_id = "declared", assumption_id = "edge",
                  expected_reference_rna = "ATGAAATAA", annotated_start1 = 1L,
                  alternative_start1 = 0L, alternative_trigger_id = "", genetic_code = 1L))
  expect_identical(top$source_effects[[1L]]$genomic_start1,
                   "18446744073709551610")
  expect_identical(top$reference_start1, "18446744073709551607")
  expect_error(run(con, reference, exons, exons,
                   sources[c(1L, 1L), ], programme), "ambiguous composition")
  expect_error(run(con, reference, exons, exons,
                   transform(sources, reference = c("T", "G")), programme), "reference mismatch")
  expect_error(run(con, reference, exons, exons,
                   transform(sources, alternate = c("CC", "C")), programme), "invalid scenario")
  expect_equal(run(con, reference, exons, exons, sources, programme)$alternate_rna,
               included$alternate_rna)
  expect_identical(list(reference, exons, sources, programme), original)
})
