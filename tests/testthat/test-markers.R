ex <- hybridclock_example()

test_that("diagnostic markers are found and oriented", {
  mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
  expect_gte(nrow(mk), 300)
  expect_true(any(mk$flip))
  expect_true(all(mk$freq2 >= 0.99 & mk$freq1 <= 0.01))
  sc <- hybrid_scores(ex$genotypes, mk)
  expect_true(all(sc$S[ex$groups == "species1"] == 0))
  expect_true(all(sc$S[ex$groups == "species2"] == 1))
  expect_true(all(sc$H[ex$groups %in% c("species1", "species2")] == 0))
  expect_true(all(sc$H <= 2 * pmin(sc$S, 1 - sc$S) + 1e-9))
})

test_that("a stricter threshold keeps fewer markers and one SNP per locus", {
  strict <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", threshold = 0)
  loose <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", threshold = 0.05)
  expect_lte(nrow(strict), nrow(loose))
  loc <- rep(paste0("L", 1:245), each = 2)
  one <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = loc)
  expect_false(anyDuplicated(loc[one$column]) > 0)
  expect_error(diagnostic_markers(ex$genotypes[, 400:490], ex$groups, "species1", "species2"), "No SNP")
})

test_that("shared loci are common in both species", {
  sh <- shared_loci(ex$genotypes, ex$groups, "species1", "species2")
  expect_gt(nrow(sh), 50)
  expect_true(all(sh$freq > 0.05 & sh$freq < 0.95))
})

test_that("marker information counts fixed markers as one each", {
  expect_equal(marker_information(data.frame(freq1 = c(0, 0.1), freq2 = c(1, 0.6))), 1.25)
})

test_that("likelihood scores handle markers that are not fixed", {
  set.seed(4)
  L <- 300; f1 <- runif(L, 0, 0.3); f2 <- runif(L, 0.7, 1)
  mk <- data.frame(column = 1:L, name = paste0("m", 1:L), flip = FALSE, freq1 = f1, freq2 = f2)
  g <- matrix(rbinom(30 * L, 1, rep(f1, each = 30)) + rbinom(30 * L, 1, rep(f2, each = 30)), 30)
  colnames(g) <- mk$name
  m <- hybrid_scores(g, mk, method = "likelihood")
  expect_equal(mean(m$S), 0.5, tolerance = 0.03)
  expect_gt(mean(m$H), 0.9)
  expect_lt(mean(hybrid_scores(g, mk)$H), 0.8)
})

test_that("microsatellite scores recover pure animals and F1s", {
  set.seed(6)
  loci <- paste0("L", 1:15); al <- 100 + 2 * (0:7)
  f1 <- lapply(1:15, function(i) { p <- c(runif(4, 1, 3), runif(4, 0, 0.2)); p / sum(p) })
  f2 <- lapply(1:15, function(i) { p <- c(runif(4, 0, 0.2), runif(4, 1, 3)); p / sum(p) })
  draw <- function(n, o1, o2) do.call(cbind, lapply(1:15, function(i) cbind(
    sample(al, n, TRUE, if (o1 == 1) f1[[i]] else f2[[i]]), sample(al, n, TRUE, if (o2 == 1) f1[[i]] else f2[[i]]))))
  a <- rbind(draw(40, 1, 1), draw(40, 2, 2), draw(20, 1, 2))
  colnames(a) <- as.vector(rbind(paste0(loci, ".a1"), paste0(loci, ".a2")))
  g <- rep(c("sp1", "sp2", "F1"), c(40, 40, 20))
  ms <- microsat_markers(a, g, "sp1", "sp2")
  sc <- hybrid_scores_microsat(a, ms)
  expect_lt(mean(sc$S[g == "sp1"]), 0.05)
  expect_gt(mean(sc$S[g == "sp2"]), 0.95)
  expect_gt(mean(sc$H[g == "F1"]), 0.85)
  s <- simulate_hybrid_zone(generations = 2, breeders = 50, microsat = ms, sample_size = 20, seed = 2)
  expect_equal(nrow(s), 2)
  cl <- classify_hybrids_microsat(a, ms)
  expect_gte(mean(cl$class[g == "F1"] == "F1"), 0.8)
  expect_gte(mean(cl$class[g == "sp1"] == "pure1"), 0.9)
})

test_that("check_markers flags markers shared by the breeding parents", {
  set.seed(2)
  n <- 30; L <- 20
  g <- matrix(1L, n, L, dimnames = list(NULL, paste0("m", 1:L)))   # F1s heterozygous everywhere
  g[, 3] <- sample(0:2, n, TRUE, prob = c(0.3, 0.4, 0.3))           # marker 3 not diagnostic in the parents
  mk <- data.frame(column = 1:L, name = colnames(g), flip = FALSE, freq1 = 0, freq2 = 1)
  chk <- check_markers(g, mk, hybrid_scores(g, mk))
  expect_equal(attr(chk, "n_f1"), n)
  expect_equal(chk$name[chk$flagged], "m3")
})
