library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)
info <- "END=2008;REF=1;RL=3;RU=CAG;REPID=SNV_AND_STR"
format <- "GT:SO:REPCN:REPCI"
sample <- "1/2:SPANNING/SPANNING:2/10:2-2/10-10"
prepare <- function(info_value = info, format_value = format, sample_value = sample, ref = "C",
                    alt = "<STR2>,<STR10>", reference = "CAG", alt_index = 1L) {
  dbWriteTable(con, "eh_input", data.frame(event_index = 1L, info = info_value,
    format = format_value, sample = sample_value, ref, alt, alt_index),
    temporary = TRUE, overwrite = TRUE)
  dbWriteTable(con, "eh_reference", data.frame(event_index = 1L,
    reference_sequence = reference), temporary = TRUE, overwrite = TRUE)
  as.list(rduckvep_prepare_expansionhunter(con, "eh_input", "eh_reference")[1L, ])
}
expect_identical(prepare()$status, "ok")
expect_equal(prepare()$alternate_components$count, 2)
expect_equal(prepare(alt_index = 2L)$alternate_components$count, 10)
expect_identical(prepare(info_value = sub("RL=3", "RL=1", info))$reason, "reference_mismatch")
expect_identical(prepare(sample_value = sub("2/10", "3/10", sample))$reason, "count_mismatch")
expect_identical(prepare(sample_value = sub("2-2", "1-3", sample))$reason, "estimated_count")
expect_identical(prepare(sample_value = sub("1/2", "0/2", sample))$reason, "uncalled_alt")
expect_identical(prepare(reference = "CAT")$reason, "reference_mismatch")
expect_identical(prepare(alt_index = 3L)$reason, "alt_index")
als <- "END=27573544;REF=3;RL=18;RU=GGCCCC"
als_format <- "GT:SO:CN:CI"
als_sample <- "1/2:SPANNING/INREPEAT:2/349:2-2/323-376"
als_call <- function(i) prepare(als, als_format, als_sample, "C", "<STR2>,<STR349>",
  strrep("GGCCCC", 3L), i)
expect_identical(als_call(1L)$status, "ok")
expect_identical(als_call(2L)$reason, "estimated_count")
expect_identical(prepare(als, paste0(als_format, ":REPCN:REPCI"),
  paste0(als_sample, ":2/349:2-2/323-376"), "C", "<STR2>,<STR349>",
  strrep("GGCCCC", 3L))$reason, "ambiguous_count_fields")
expect_identical(als_call(1e20)$reason, "alt_index")
expect_identical(prepare("END=6000;REF=6000;RL=6000;RU=A", "GT:SO:CN:CI",
  "1:SPANNING:2:2-2", "A", "<STR2>", strrep("A", 6000L))$reason,
  "allele_capacity")
dbDisconnect(con, shutdown = TRUE)
