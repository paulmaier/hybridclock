check_snp_data <- function(d, n, m, groups) {
  expect_true(is.integer(d$genotypes))
  expect_equal(dim(d$genotypes), c(n, m))
  expect_length(d$groups, n)
  expect_setequal(unique(d$groups), groups)
  expect_true(all(d$genotypes %in% c(0L, 1L, 2L, NA)))
  expect_false(anyDuplicated(rownames(d$genotypes)) > 0)
  expect_true(is.character(d$description) && length(d$description) == 1)
  expect_setequal(d$truth$group, groups)
  expect_equal(d$truth$n[match(names(table(d$groups)), d$truth$group)], as.vector(table(d$groups)))
}

test_that("wolves load with pedigree truth", {
  check_snp_data(wolves, 56, 486, c("ref1", "ref2", "F1", "BC1_MB", "BC2_MB"))
  expect_equal(nrow(wolves$map), ncol(wolves$genotypes))
  expect_identical(wolves$map$snp, colnames(wolves$genotypes))
  expect_setequal(unique(wolves$map$chromosome), 1:38)
  expect_identical(wolves$animals$id, rownames(wolves$genotypes))
  expect_equal(as.vector(table(wolves$groups)[c("ref1", "ref2", "F1", "BC1_MB", "BC2_MB")]), c(31, 6, 5, 11, 3))
  expect_equal(as.vector(tapply(wolves$animals$hybrid_age, wolves$animals$group, unique)[c("F1", "BC1_MB", "BC2_MB")]),
               c(1, 1.5, 1.75))
  mk <- diagnostic_markers(wolves$genotypes, wolves$groups, "ref1", "ref2")
  expect_equal(nrow(mk), ncol(wolves$genotypes))
  expect_false(any(mk$flip))
})

test_that("chickens load as one F8 group with a map", {
  check_snp_data(chickens, 572, 514, "F8")
  expect_identical(chickens$map$snp, colnames(chickens$genotypes))
  expect_true(is.integer(chickens$map$chromosome))
  expect_equal(chickens$truth$hybrid_age, 8)
  expect_gt(mean(chickens$genotypes, na.rm = TRUE) / 2, 0.4)
})

test_that("bettongs load with founders and yearly groups", {
  check_snp_data(bettongs, 215, 500, c("ref1", "ref2", paste0("Matuwa", 2010:2018)))
  expect_equal(as.vector(table(bettongs$groups)[c("ref1", "ref2")]), c(22, 24))
  expect_identical(bettongs$animals$id, rownames(bettongs$genotypes))
  expect_equal(bettongs$truth$max_hybrid_age, c(0, 0, 0:8))
  mk <- diagnostic_markers(bettongs$genotypes, bettongs$groups, "ref1", "ref2")
  expect_equal(nrow(mk), 500)
})

test_that("possums load as raw microsatellite alleles", {
  a <- possums$alleles
  expect_true(is.integer(a))
  expect_equal(dim(a), c(536, 48))
  expect_equal(length(unique(sub("\\.a[12]$", "", colnames(a)))), 24)
  expect_true(all(grepl("\\.a[12]$", colnames(a))))
  grp <- c("ref1", "ref2", paste0("Buller", 2011:2015), "Higginbotham_males_2011", "TimmsSpur_males_2014")
  expect_setequal(unique(possums$groups), grp)
  expect_equal(as.vector(table(possums$groups)[grp]), c(29, 104, 40, 43, 67, 103, 138, 6, 6))
  expect_identical(possums$animals$id, rownames(a))
  expect_equal(possums$truth$max_hybrid_age[match(paste0("Buller", 2011:2015), possums$truth$group)], 0:4)
  ms <- microsat_markers(a, possums$groups, "ref1", "ref2")
  expect_length(ms$loci, 24)
})

test_that("carp load with native-range references and three sites", {
  check_snp_data(carp, 1354, 57, c("ref1", "ref2", "IMAR", "PL26", "MKY"))
  expect_equal(as.vector(table(carp$groups)[c("ref1", "ref2", "IMAR", "PL26", "MKY")]), c(16, 20, 545, 365, 408))
  expect_identical(carp$animals$id, rownames(carp$genotypes))
  expect_equal(carp$truth$max_hybrid_age, c(0, 0, 12, 12, 12))
  mk <- diagnostic_markers(carp$genotypes, carp$groups, "ref1", "ref2")
  expect_equal(nrow(mk), 57)
})
