test_that("first generations carry the right ages", {
  s <- simulate_hybrid_zone(generations = 2, breeders = 200, start_mix = 0.5, n_markers = 100,
                            sample_size = 200, seed = 1)
  expect_equal(s$hybrid_age_both[1], 1)
  expect_equal(s$hybrid_age_rare[1], 1)
  expect_equal(s$backcrosses[1], 0)
  expect_equal(s$recent_share[1], 1)
  # generation 2 hybrids are F2s (age 2) and first backcrosses (age 1.5)
  expect_gt(s$hybrid_age_both[2], 1.5)
  expect_lt(s$hybrid_age_both[2], 2)
  expect_equal(s$recent_share[2], 1)
})

test_that("arriving animals start with pure labels", {
  st <- list(H1 = matrix(0L, 10, 5), H2 = matrix(1L, 10, 5), anc = rep(0.5, 10), age = rep(4, 10), t1 = rep(4, 10),
             t2 = rep(4, 10), bc = rep(2, 10), cls = rep(4L, 10), fam = 1:10, pool = rep(1L, 10))
  set.seed(3)
  out <- immigrate(st, c(0.5, 0.5), 5, NULL, 1)
  new <- out$anc %in% c(0, 1)
  expect_true(any(new))
  expect_true(all(out$age[new] == 0 & out$bc[new] == 0 & out$cls[new] == 0L))
  G <- out$H1 + out$H2
  expect_true(all(rowMeans(G[new, , drop = FALSE]) %in% c(0, 2)))
  expect_true(all(G[new, ] / 2 == out$anc[new]))
})

test_that("backcrossing keeps hybrid age low but ancestry age grows", {
  s <- simulate_hybrid_zone(generations = 25, breeders = 150, start_mix = 0.2, arrivals = c(0.15, 0),
                            n_markers = 200, seed = 4)
  late <- s[s$generation > 15, ]
  expect_gt(mean(late$hybrid_age_rare, na.rm = TRUE), mean(late$hybrid_age_both, na.rm = TRUE))
  expect_gt(mean(late$backcrosses, na.rm = TRUE), 1)
})

test_that("seeds reproduce and leave the global random state alone", {
  set.seed(10); before <- runif(1); set.seed(10)
  a <- simulate_hybrid_zone(generations = 3, breeders = 30, n_markers = 50, seed = 5)
  after <- runif(1)
  b <- simulate_hybrid_zone(generations = 3, breeders = 30, n_markers = 50, seed = 5)
  expect_identical(a, b)
  expect_identical(before, after)
})

test_that("samples can be returned", {
  s <- simulate_hybrid_zone(generations = 3, breeders = 30, n_markers = 50, shared_freq = rep(0.5, 10),
                            sample_at = 3, seed = 2)
  smp <- attr(s, "samples")[["3"]]
  expect_equal(dim(smp$genotypes), c(20, 60))
  expect_length(smp$S, 20)
})

test_that("arrival schedules give the current rate", {
  arr <- cbind(0, c(rep(0.02, 5), rep(0, 5)))
  s <- simulate_hybrid_zone(generations = 10, breeders = 50, start_mix = 0.1, arrivals = arr,
                            n_markers = 50, seed = 7)
  expect_equal(recent_rate(arr[, 2], 1), 0)
  expect_equal(recent_rate(arr[, 2], 4), 0.02)
  expect_equal(recent_rate(arr[, 2], 7), 0.02 * 2 / 3)
  expect_equal(recent_rate(arr[, 2], 10), 0)
  expect_equal(s$arriving[s$rare_species == 2 & s$generation %in% 2:6], rep(1, sum(s$rare_species == 2 & s$generation %in% 2:6)))
  expect_error(simulate_hybrid_zone(generations = 10, arrivals = arr[1:5, ], n_markers = 20))
})


test_that("linkage makes F2s more variable, and no map means unlinked", {
  tight <- data.frame(chromosome = rep(1:2, each = 200), position = rep(seq(0, 50, length.out = 200), 2))
  f2 <- function(map) simulate_hybrid_zone(generations = 1, breeders = 200, founders = "F1", n_markers = 400,
                                           sample_size = 200, map = map, seed = 3)
  expect_gt(f2(tight)$hyb_sd_S, 3 * f2(NULL)$hyb_sd_S)
  expect_equal(f2(tight)$mean_H, 0.5, tolerance = 0.05)
  expect_null(recombination(NULL, 10, 5))
  r <- recombination(data.frame(chromosome = c(1, 1, 2), position = c(0, 10, 0)), 3, 2)
  expect_equal(r, c((1 - exp(-0.2)) / 2, 0.5, 0.5, 0.5))
  m <- random_map(100, chromosomes = 5, length_cM = 500)
  expect_equal(nrow(m), 100)
  expect_false(is.unsorted(m$chromosome))
  expect_true(all(m$position <= 100))
})

test_that("hybrid founders start with the right histories", {
  s <- simulate_hybrid_zone(generations = 1, breeders = 100, founders = "BC1_1", n_markers = 200,
                            sample_size = 100, seed = 4)
  expect_equal(s$hybrid_age_both, 2.5)
  expect_equal(s$mean_S, 0.25, tolerance = 0.02)
})

test_that("heterozygote dropout lowers H without changing S", {
  set.seed(9)
  G <- matrix(rbinom(400 * 300, 2, 0.3), 400)
  o <- observe_genotypes(G, 300, missing = 0, het_dropout = 0.3)
  expect_equal(mean(o) / 2, mean(G) / 2, tolerance = 0.01)
  expect_equal(mean(o == 1L), 0.7 * mean(G == 1L), tolerance = 0.01)
  g <- observe_genotypes(matrix(1L, 5, 50), 50, missing = 0.2, het_dropout = 1)
  expect_true(all(g %in% c(0L, 2L, NA)))
  expect_equal(sum(is.na(g)), 50)
  d <- simulate_hybrid_zone(generations = 3, breeders = 100, n_markers = 100, shared_freq = rep(0.5, 100),
                            sample_size = 100, het_dropout = 0.3, seed = 8)
  expect_gt(d$shared_fis[3], 0.15)
  a <- simulate_hybrid_zone(generations = 3, breeders = 100, n_markers = 100, assort = 0.8, seed = 8)
  expect_equal(nrow(a), 3)
})

test_that("founders can carry the other species' alleles", {
  set.seed(1)
  h <- pure_haplotypes(200, 0L, 50, NULL, impurity = 0.1)
  expect_lt(abs(mean(h$H1) - 0.1), 0.02)
  expect_equal(mean(pure_haplotypes(50, 1L, 20, NULL)$H2), 1)
  f <- list(freq1 = rep(0, 50), freq2 = rep(1, 50))
  expect_lt(abs(mean(pure_haplotypes(200, 1L, 50, NULL, f, impurity = 0.2)$H1) - 0.8), 0.02)
  s <- simulate_hybrid_zone(generations = 3, breeders = 40, n_markers = 100, founders = "F1",
                            founder_impurity = 0.05, seed = 2)
  expect_equal(nrow(s), 3)
  ref <- build_reference(n_sims = 4, generations = 3, breeders = c(20, 40), n_markers = 50,
                         founder_impurity = c(0, 0.1), seed = 3)
  expect_true(all(ref$founder_impurity >= 0 & ref$founder_impurity <= 0.1))
  expect_true(any(ref$founder_impurity > 0))
})

test_that("founder samples give F1 and backcross families", {
  s <- simulate_hybrid_zone(generations = 2, breeders = 30, n_markers = 80, founders = "F1",
                            sample_founders = TRUE, seed = 5)
  expect_equal(s$generation, 0:2)
  expect_equal(s$hybrid_age_both[1], 1)
  expect_equal(s$hybrid_age_both[2], 2)
  b <- simulate_hybrid_zone(generations = 1, breeders = 30, n_markers = 80, founders = "BC1_2",
                            sample_founders = TRUE, seed = 6)
  expect_equal(b$hybrid_age_both[1], 1.5)
  expect_equal(b$backcrosses[1], 1)
  p <- simulate_hybrid_zone(generations = 2, breeders = 30, n_markers = 80, sample_founders = TRUE, seed = 7)
  expect_equal(p$generation, 1:2)
})
