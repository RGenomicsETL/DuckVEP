library(tinytest)

info <- "END=2008;REF=1;RL=3;RU=CAG;REPID=SNV_AND_STR"
format <- "GT:SO:REPCN:REPCI"
sample <- "1/2:SPANNING/SPANNING:2/10:2-2/10-10"
prepare <- function(...) rduckvep_prepare_expansionhunter(info, format, sample,
  "C", "<STR2>,<STR10>", "CAG", ...)
expect_identical(prepare(alt_index = 1L)$status, "ok")
expect_equal(prepare(alt_index = 1L)$alternate_components$count, 2)
expect_equal(prepare(alt_index = 2L)$alternate_components$count, 10)
expect_identical(rduckvep_prepare_expansionhunter(sub("RL=3", "RL=1", info),
  format, sample, "C", "<STR2>,<STR10>", "CAG")$reason, "reference_mismatch")
expect_identical(rduckvep_prepare_expansionhunter(info, format,
  sub("2/10", "3/10", sample), "C", "<STR2>,<STR10>", "CAG")$reason,
  "count_mismatch")
expect_identical(rduckvep_prepare_expansionhunter(info, format,
  sub("2-2", "1-3", sample), "C", "<STR2>,<STR10>", "CAG")$reason,
  "estimated_count")
expect_identical(rduckvep_prepare_expansionhunter(info, format,
  sub("1/2", "0/2", sample), "C", "<STR2>,<STR10>", "CAG")$reason,
  "uncalled_alt")
expect_identical(rduckvep_prepare_expansionhunter(info, format, sample,
  "C", "<STR2>,<STR10>", "CAT")$reason, "reference_mismatch")
expect_identical(rduckvep_prepare_expansionhunter(info, format, sample,
  "C", "<STR2>,<STR10>", "CAG", alt_index = 3L)$reason, "alt_index")
als <- "END=27573544;REF=3;RL=18;RU=GGCCCC"
als_format <- "GT:SO:CN:CI"
als_sample <- "1/2:SPANNING/INREPEAT:2/349:2-2/323-376"
als_call <- function(i) rduckvep_prepare_expansionhunter(als, als_format,
  als_sample, "C", "<STR2>,<STR349>", strrep("GGCCCC", 3L), alt_index = i)
expect_identical(als_call(1L)$status, "ok")
expect_identical(als_call(2L)$reason, "estimated_count")
expect_identical(rduckvep_prepare_expansionhunter(als,
  paste0(als_format, ":REPCN:REPCI"), paste0(als_sample, ":2/349:2-2/323-376"),
  "C", "<STR2>,<STR349>", strrep("GGCCCC", 3L))$reason,
  "ambiguous_count_fields")
expect_identical(als_call(1e20)$reason, "alt_index")
expect_identical(rduckvep_prepare_expansionhunter(
  "END=6000;REF=6000;RL=6000;RU=A", "GT:SO:CN:CI",
  "1:SPANNING:2:2-2", "A", "<STR2>", strrep("A", 6000L))$reason,
  "allele_capacity")
