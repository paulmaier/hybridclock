test_that("summaries do not depend on which species is called species 2", {
  set.seed(1)
  S <- c(0, 0, runif(8, 0.05, 0.3)); H <- c(0, 0, runif(8, 0.1, 0.5))
  a <- sample_summaries(S, H); b <- sample_summaries(1 - S, H)
  expect_equal(a, b)
  expect_named(a, feature_names())
  expect_equal(unname(a["n_hybrids"]), 8)
  expect_equal(unname(a["share_pure_common"]), 0.2)
})

test_that("inbreeding shows up below the floor", {
  S <- rep(0.25, 10)
  outbred <- sample_summaries(S, rep(2 * 0.25 * 0.75, 10))
  inbred <- sample_summaries(S, rep(0.25, 10))
  expect_equal(unname(outbred["het_floor_ratio"]), 1)
  expect_lt(inbred["het_floor_ratio"], 0.9)
  expect_equal(unname(inbred["below_floor"]), 1)
})

test_that("genotype summaries are computed when genotypes are given", {
  set.seed(2)
  G <- matrix(rbinom(20 * 50, 2, 0.1), 20)
  sh <- matrix(rbinom(20 * 40, 2, 0.5), 20)
  v <- sample_summaries(rowMeans(G) / 2, rowMeans(G == 1), diag = G, shared = sh)
  expect_true(all(is.finite(v[c("marker_spread", "shared_fis", "shared_ld")])))
})

test_that("class shares are added when marker counts are given", {
  S <- c(0, 1, 0.5, 0.25, 0.5, 0.1); H <- c(0, 0, 1, 0.5, 0.5, 0.12)
  v <- sample_summaries(S, H, n_markers = rep(800, 6))
  cls <- v[grep("^class_", names(v))]
  expect_equal(sum(cls), 1)
  expect_equal(unname(v["class_F1"]), 1 / 6)
  expect_true(all(is.na(sample_summaries(S, H)[grep("^class_", names(v))])))
})
