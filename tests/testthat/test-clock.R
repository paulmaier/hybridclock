test_that("a clock can be trained and applied", {
  skip_on_cran()
  ex <- hybridclock_example()
  mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
  sc <- hybrid_scores(ex$genotypes, mk)
  ref <- build_reference(n_sims = 60, generations = 12, n_markers = nrow(mk), breeders = c(30, 120),
                         sample_size = 20, seed = 3)
  expect_equal(attr(ref, "settings")$sample_size, 20)
  clock <- train_clock(ref, targets = c("hybrid_age_both", "backcrosses", "arriving"), committees = 2, max_rows = 2000)
  expect_s3_class(clock, "hybridclock_model")
  expect_output(print(clock), "hybrid_age_both")
  hyb <- ex$groups %in% c("site_A", "site_B")
  est <- estimate_history(clock, sc[hyb, ], ex$groups[hyb], diag = orient_genotypes(ex$genotypes, mk)[hyb, ],
                          repeats = 3)
  expect_equal(nrow(est), 2)
  expect_true(all(est$hybrid_age_both >= est$hybrid_age_both_lower & est$hybrid_age_both <= est$hybrid_age_both_upper))
  expect_error(estimate_history(clock, sc[hyb, ], ex$groups[hyb]), "diag")
  small <- ex$groups == "site_A"
  expect_warning(est2 <- estimate_history(clock, sc[small, ][1:10, ], rep("A", 10),
                                          diag = orient_genotypes(ex$genotypes, mk)[small, ][1:10, ]), "fewer")
  expect_true(is.na(est2$hybrid_age_both))
  p <- predict(clock, ref[1:5, ])
  expect_named(p, c("hybrid_age_both", "hybrid_age_both_lower", "hybrid_age_both_upper",
                    "backcrosses", "backcrosses_lower", "backcrosses_upper", "arriving"))
  expect_true(all(p$arriving >= 0 & p$arriving <= 1))
  expect_true(all(est$arriving >= 0 & est$arriving <= 1))
  expect_true(all(est$arriving_note %in% c("", "uncertain: young swarm")))
  expect_true(all(est$method %in% c("individual classes", "population clock")))
  expect_true(all(c("clock_hybrid_age_both", "classes_hybrid_age_both", "early_share", "f2_share") %in% names(est)))
  expect_output(print(clock), "share called correctly")
  expect_true(all(p$hybrid_age_both_lower <= p$hybrid_age_both_upper))
})

test_that("references run in parallel and match serial runs", {
  skip_on_cran()
  skip_on_os("windows")
  a <- build_reference(n_sims = 4, generations = 3, n_markers = 50, seed = 9)
  b <- build_reference(n_sims = 4, generations = 3, n_markers = 50, seed = 9, cores = 2)
  expect_equal(a, b)
  expect_true(all(a$timing1 %in% c("none", "constant", "stopped", "started")))
  expect_setequal(unique(a$kind), intersect(c("closed", "one", "both"), unique(a$kind)))
})


test_that("class-based ages follow the early classes", {
  x <- as.data.frame(as.list(setNames(rep(0, length(class_features)), class_features)))
  x$class_F1 <- 0.5; x$class_BC1_common <- 0.5
  ca <- class_ages(x)
  expect_equal(ca$early_share, 1)
  expect_equal(ca$hybrid_age_both, 1.25)
  expect_equal(ca$backcrosses, 0.5)
  expect_equal(ca$f2_share, 0)
  expect_equal(orient_classes(c("pure1", "BC1_species2", "F2"), "species1"), c("pure_rare", "BC1_common", "F2"))
})

test_that("closed references skip targets that cannot vary", {
  skip_on_cran()
  ref <- build_reference(n_sims = 40, generations = 6, n_markers = 60, breeders = c(20, 60), founders = "F1",
                         kinds = c(closed = 1, one = 0, both = 0), seed = 4)
  expect_warning(clock <- train_clock(ref, committees = 2, max_rows = 1000), "do not vary")
  expect_false(any(c("arriving", "arrival_rate") %in% clock$targets))
  expect_true("hybrid_age_both" %in% clock$targets)
})

test_that("sample types, F1 share and stratified intervals behave", {
  cls <- paste0("class_", c(early_classes$class, "later"))
  x <- as.data.frame(matrix(0, 3, length(cls), dimnames = list(NULL, cls)))
  x$class_F2[1] <- 0.9; x$class_later[1] <- 0.1      # F2-like
  x$class_F1[2] <- 1                                  # F1s
  x$class_later[3] <- 1                               # later generations
  expect_equal(sample_type(x), c("f2_like", "other", "other"))
  expect_equal(class_ages(x)$f1_share[2], 1)
  expect_equal(sample_type(data.frame(a = 1:2)), c("other", "other"))
  calib <- conformal_fit(rnorm(600), rnorm(600, sd = 0.1), 0.9, 4, rep(c("a", "b"), c(400, 200)))
  expect_named(calib$strata, c("a", "b"))
  iv <- apply_conformal(c(0, 0), calib, c("a", "zz"))
  expect_true(all(iv$lower < 0 & iv$upper > 0))
  old <- conformal_bins(rnorm(100), rnorm(100, sd = 0.1), 0.9, 4)   # clocks trained before sample types
  expect_length(apply_conformal(0, old)$lower, 1)
})

test_that("estimates report how much they owe to the data", {
  skip_on_cran()
  ex <- hybridclock_example()
  mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
  ref <- build_reference(n_sims = 60, generations = 12, n_markers = nrow(mk), breeders = c(30, 120), seed = 3)
  clock <- train_clock(ref, targets = "hybrid_age_both", committees = 2, max_rows = 2000)
  expect_length(clock$age_range, 2)
  hyb <- ex$groups %in% c("site_A", "site_B")
  est <- estimate_history(clock, hybrid_scores(ex$genotypes, mk)[hyb, ], ex$groups[hyb],
                          diag = orient_genotypes(ex$genotypes, mk)[hyb, ], repeats = 3)
  expect_true(all(est$information >= 0 & est$information <= 1, na.rm = TRUE))
})
