test_that("the triangle plot draws without error", {
  ex <- hybridclock_example()
  mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f)
  expect_null(plot_triangle(hybrid_scores(ex$genotypes, mk), ex$groups))
  expect_null(plot_triangle(hybrid_scores(ex$genotypes, mk), classes = FALSE))
  grDevices::dev.off(); unlink(f)
})
