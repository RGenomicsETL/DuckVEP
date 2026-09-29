library(tinytest)
library(DBI)
con <- dbConnect(duckdb::duckdb(shared_home = FALSE,
  config = list(allow_unsigned_extensions = "true")))
rduckvep_load(con)
records <- data.frame(
  event_index = 0:7,
  chrom = c("4", "4", "6", "6", "9", "9", "10", "10"),
  pos = c(1806934, 1727977, 117566854, 117321394, 10, 20, 30, 31),
  id = c("A5", "A3", "B5", "B3", "C1", "C2", "D1", "D2"),
  ref = "N",
  alt = c("N[4:1727977[", "]4:1806934]N", "N]6:117321394]", "N]6:117566854]",
          "N[9:20[", "N[9:10[", "N[10:30[", "]10:30]N"),
  info = c("MATEID=A3;EVENT=e", "MATEID=A5;EVENT=e", "MATEID=B3", "MATEID=B5",
           "MATEID=C2", "MATEID=C1", "MATEID=D2", "MATEID=D1"),
  stringsAsFactors = FALSE)
dbWriteTable(con, "bnd_input", records, temporary = TRUE)
pairs <- rduckvep_prepare_breakend_pairs(con, "bnd_input")
expect_equal(nrow(pairs), 8L)
expect_equal(pairs$event_index, 0:7)
# Two physical records per pair stay separate rows sharing one pair key.
expect_identical(pairs$reason[1:4], rep("reciprocal", 4L)[1:4])
expect_identical(pairs$pair_key[1:2], c("0|1", "0|1"))
expect_identical(pairs$reason[5:6], rep("mate_orientation_conflict", 2L))
expect_identical(pairs$reason[7:8], rep("mate_coordinate_conflict", 2L))
expect_true(all(pairs$fusion_status == "not_asserted"))
expect_true(all(pairs$phase_status == "not_evaluated"))
expect_true(all(!pairs$orientation_reciprocal[5:6]))
expect_true(all(pairs$id_reciprocal[5:8]))

dbWriteTable(con, "bnd_pairs", pairs, temporary = TRUE)
dbWriteTable(con, "bnd_genes", data.frame(event_index = c(0L, 1L, 4L, 5L),
  gene_id = c("FGFR3", "TACC3", "X", "Y")), temporary = TRUE)
fusion <- rduckvep_prepare_breakend_fusion(con, "bnd_pairs", "bnd_genes")
expect_equal(nrow(fusion), 8L)
expect_identical(fusion$status[1:2], rep("candidate_partner_genes", 2L))
expect_identical(fusion$status[5:6], rep("candidate_orientation_conflict", 2L))
expect_identical(fusion$status[7:8], rep("identity_unproven", 2L))
expect_identical(fusion$reason[7], "mate_coordinate_conflict")
expect_true(all(!fusion$fusion_asserted))
expect_identical(fusion$endpoint_genes[[1L]], "FGFR3")
expect_identical(fusion$mate_endpoint_genes[[1L]], "TACC3")
expect_error(rduckvep_prepare_breakend_pairs(con, "bnd_input", bogus = 1L), "unknown option")
dbDisconnect(con, shutdown = TRUE)
