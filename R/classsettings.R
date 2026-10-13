#' Settings for the classes that choose between the two tiers
#'
#' The share of hybrids in each early class is part of the summaries the
#' clock learns from, and it decides whether a population is dated by its
#' individual classes or by the clock (see [estimate_history()]). Those
#' classes must be made the same way in the simulated reference and in the
#' data, so the settings are given once, to [build_reference()], and are
#' then used for both.
#'
#' With markers that are not fixed, class probabilities are averaged over
#' markers using the reference allele frequencies, which needs only a few
#' summary moments of the frequencies and keeps the classification fast
#' enough to run on every simulated sample.
#'
#' @param error Genotyping error rate assumed by the classes.
#' @param het_dropout Share of heterozygotes called homozygous (see
#'   [classify_hybrids()]); for RAD or GBS data where F1s show clearly less
#'   than full heterozygosity.
#' @param effective_markers Optional effective number of independent markers
#'   (see [effective_markers()]), for many linked markers.
#' @param markers Optional output of [diagnostic_markers()] (columns `freq1`,
#'   `freq2`) for markers that are not fixed; with it, animals are classified
#'   from their genotype counts at the markers.
#' @param alpha Significance level of the fit test for "later".
#' @return A list of class `hybridclock_class_settings`.
#' @examples
#' class_settings(het_dropout = 0.2)
#' @export
class_settings <- function(error = 0.005, het_dropout = 0, effective_markers = NULL, markers = NULL,
                           alpha = 0.01) {
  m <- if (is.null(markers)) c(x = 0, d = 1, xx = 0, xd = 0, dd = 1) else {
    x <- pmin(pmax(markers$freq1, 0), 1); d <- pmin(pmax(markers$freq2, 0), 1) - x
    c(x = mean(x), d = mean(d), xx = mean(x * x), xd = mean(x * d), dd = mean(d * d))
  }
  structure(list(error = error, het_dropout = het_dropout, effective_markers = effective_markers,
                 moments = m, counts = !is.null(markers), alpha = alpha),
            class = "hybridclock_class_settings")
}

#' @export
print.hybridclock_class_settings <- function(x, ...) {
  cat("class settings: error", x$error, " het_dropout", x$het_dropout,
      " effective markers", if (is.null(x$effective_markers)) "none" else x$effective_markers,
      if (x$counts) " (markers not fixed: genotype counts, averaged frequencies)" else "", "\n")
  invisible(x)
}

# genotype counts (n0, n1, n2) per animal: from the genotypes when the
# settings ask for counts and genotypes are given, otherwise from S and H
class_counts <- function(S, H, n_markers, diag, cs) {
  if (isTRUE(cs$counts) && !is.null(diag)) {
    return(cbind(rowSums(diag == 0L, na.rm = TRUE), rowSums(diag == 1L, na.rm = TRUE), rowSums(diag == 2L, na.rm = TRUE)))
  }
  n1 <- round(n_markers * H); n2 <- pmax(0, round(n_markers * (2 * S - H) / 2))
  cbind(pmax(0, n_markers - n1 - n2), n1, n2)
}

# classes from genotype counts with marker-averaged class probabilities; the
# same function classifies simulated samples and data
classify_counts <- function(G, cs = class_settings()) {
  if (is.null(cs)) cs <- class_settings()
  m <- cs$moments; e <- cs$error; dr <- cs$het_dropout
  loglik <- function(a, b) {
    p2 <- m[["xx"]] + (a + b) * m[["xd"]] + a * b * m[["dd"]]
    p0 <- 1 - 2 * m[["x"]] - (a + b) * m[["d"]] + p2
    p1 <- 1 - p0 - p2
    p0 <- p0 + p1 * dr / 2; p2 <- p2 + p1 * dr / 2; p1 <- p1 * (1 - dr)
    G %*% log(pmax((1 - e) * rbind(p0, p1, p2) + e / 3, 1e-300))
  }
  ca <- vapply(hybrid_classes, `[`, 0, 1); cb <- vapply(hybrid_classes, `[`, 0, 2)
  ll <- loglik(ca, cb); colnames(ll) <- names(hybrid_classes)
  free <- apply(loglik(free_grid$a, free_grid$b), 1, max)
  if (!is.null(cs$effective_markers)) {
    w <- pmin(1, cs$effective_markers / pmax(rowSums(G), 1)); ll <- ll * w; free <- free * w
  }
  best <- max.col(ll, ties.method = "first")
  fit_p <- stats::pchisq(2 * pmax(free - ll[cbind(seq_len(nrow(ll)), best)], 0), df = 2, lower.tail = FALSE)
  ifelse(fit_p < cs$alpha, "later", colnames(ll)[best])
}
