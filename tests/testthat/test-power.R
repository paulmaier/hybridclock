test_that("marker power follows Boecklen and Howard", {
  p <- marker_power(5, 1)
  expect_equal(p$p_missed, 0.5^5)
  n <- markers_needed(c(1, 5, 10))
  expect_equal(n$n_markers, c(5, 95, 3067))
})
