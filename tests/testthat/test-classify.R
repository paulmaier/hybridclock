sim_class <- function(a, b, n_animals = 30, L = 500) {
  g <- matrix(rbinom(n_animals * L, 1, a) + rbinom(n_animals * L, 1, b), n_animals)
  data.frame(S = rowMeans(g) / 2, H = rowMeans(g == 1), n_markers = L)
}

test_that("recent classes are recognized", {
  set.seed(5)
  for (cl in c("F1", "F2", "BC1_species1", "BC2_species2", "pure1")) {
    ab <- hybrid_classes[[cl]]
    out <- classify_hybrids(sim_class(ab[1], ab[2]))
    expect_gte(mean(out$class == cl), 0.9)
  }
})

test_that("later generations and inbred animals do not fit early classes", {
  set.seed(6)
  later <- classify_hybrids(sim_class(0.25, 0.25))
  expect_gte(mean(later$class == "later"), 0.9)
  inbred <- data.frame(S = rep(0.3, 10), H = rep(0.3, 10), n_markers = 500)
  expect_true(all(classify_hybrids(inbred)$class == "later"))
})


test_that("linked markers are not rejected once their effective number is given", {
  set.seed(7)
  # a backcross whose F1 parent passed on ancestry in two big blocks
  g <- c(rep(1L, 160), rep(0L, 40))
  sc <- data.frame(S = mean(g) / 2, H = mean(g == 1L), n_markers = 200)
  expect_equal(classify_hybrids(sc)$class, "later")
  expect_equal(classify_hybrids(sc, effective_markers = effective_markers(2, 100))$class, "BC1_species1")
  expect_equal(effective_markers(37, 2500), 62)
})


test_that("allele dropout can be allowed for", {
  set.seed(8)
  g <- matrix(ifelse(runif(20 * 300) < 0.2, sample(c(0L, 2L), 20 * 300, TRUE), 1L), 20)
  sc <- data.frame(S = rowMeans(g) / 2, H = rowMeans(g == 1L), n_markers = 300)
  expect_true(all(classify_hybrids(sc)$class != "F1"))
  expect_gte(mean(classify_hybrids(sc, het_dropout = 0.2)$class == "F1"), 0.9)
})

test_that("class settings reproduce the default classifier and handle markers that are not fixed", {
  ex <- hybridclock_example()
  mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
  sc <- hybrid_scores(ex$genotypes, mk)
  G <- class_counts(sc$S, sc$H, sc$n_markers, NULL, class_settings())
  expect_equal(classify_counts(G, class_settings()), classify_hybrids(sc)$class)
  # F1s at markers with frequencies 0.1 / 0.9: heterozygous at only about 82% of markers
  set.seed(4); L <- 200
  mk2 <- data.frame(freq1 = rep(0.1, L), freq2 = rep(0.9, L))
  g <- matrix(rbinom(30 * L, 1, 0.9) + rbinom(30 * L, 1, 0.1), 30, L)   # one allele from each species
  cs_fixed <- class_settings(); cs_freq <- class_settings(markers = mk2)
  k_fixed <- classify_counts(class_counts(NULL, NULL, NULL, g, class_settings(markers = mk2)), cs_fixed)
  k_freq <- classify_counts(class_counts(NULL, NULL, NULL, g, cs_freq), cs_freq)
  expect_gt(mean(k_freq == "F1"), 0.9)
  expect_lt(mean(k_fixed == "F1"), 0.5)
  expect_output(print(cs_freq), "not fixed")
})
