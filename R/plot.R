#' Plot animals on the hybrid triangle
#'
#' Plots interspecies heterozygosity (H) against hybrid index (S). Possible
#' values fill a triangle with pure animals at the bottom corners and F1s at
#' the top. The expected positions of first- and second-generation classes
#' and second backcrosses are marked. Animals below the dashed curve
#' `H = 2S(1 - S)` are less heterozygous than any offspring of unrelated
#' parents can be on average, a sign of inbreeding.
#'
#' @param scores Output of [hybrid_scores()].
#' @param group Optional grouping for colors.
#' @param col Colors, one per group.
#' @param classes Mark the expected class positions.
#' @param ... Passed to [graphics::points()].
#' @return `NULL`, invisibly.
#' @examples
#' ex <- hybridclock_example()
#' mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
#' plot_triangle(hybrid_scores(ex$genotypes, mk), ex$groups)
#' @export
plot_triangle <- function(scores, group = NULL, col = NULL, classes = TRUE, ...) {
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1.05), xlab = "Hybrid index (S)",
                 ylab = "Interspecies heterozygosity (H)", las = 1, bty = "n")
  graphics::polygon(c(0, 0.5, 1), c(0, 1, 0), border = "grey60")
  s <- seq(0, 1, length.out = 200)
  graphics::lines(s, 2 * s * (1 - s), lty = 2, col = "grey60")
  if (classes) {
    pos <- t(vapply(hybrid_classes, function(cl) c(S = mean(cl), H = cl[1] + cl[2] - 2 * cl[1] * cl[2]), numeric(2)))
    graphics::points(pos[, "S"], pos[, "H"], pch = 3, col = "grey40")
    lab <- sub("_species", " sp", rownames(pos))
    graphics::text(pos[, "S"], pos[, "H"], lab, pos = 3, cex = 0.7, col = "grey40")
  }
  if (is.null(group)) group <- rep("all", nrow(scores))
  g <- factor(group)
  if (is.null(col)) col <- grDevices::hcl.colors(nlevels(g), "viridis")
  graphics::points(scores$S, scores$H, col = col[as.integer(g)], pch = 19, ...)
  if (nlevels(g) > 1) graphics::legend("topright", levels(g), col = col, pch = 19, bty = "n", cex = 0.8)
  invisible(NULL)
}
