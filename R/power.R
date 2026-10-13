#' Chance that a backcrossed hybrid looks pure
#'
#' After `k` generations of backcrossing to one species, each fixed
#' diagnostic marker is heterozygous with probability `0.5^k`, so the chance
#' that none of `n` markers is heterozygous (and the hybrid passes as pure)
#' is `(1 - 0.5^k)^n` (Boecklen and Howard 1997).
#'
#' @param n_markers Number of diagnostic markers (vector allowed).
#' @param generations Generations of backcrossing (1 is a first backcross).
#' @return A data.frame of `generations`, `n_markers` and `p_missed`.
#' @references Boecklen WJ, Howard DJ (1997) Genetic analysis of hybrid zones:
#'   numbers of markers and power of resolution. Ecology 78:2611-2616.
#' @examples
#' marker_power(c(100, 1000), generations = 1:8)
#' markers_needed(1:10)
#' @export
marker_power <- function(n_markers, generations = 1:10) {
  g <- expand.grid(generations = generations, n_markers = n_markers)
  g$p_missed <- (1 - 0.5^g$generations)^g$n_markers
  g
}

#' @rdname marker_power
#' @param error Acceptable chance of missing a hybrid.
#' @return For `markers_needed()`, the smallest number of markers that keeps
#'   the chance of missing a hybrid of each backcross generation below `error`.
#' @export
markers_needed <- function(generations = 1:10, error = 0.05) {
  data.frame(generations = generations,
             n_markers = ceiling(log(error) / log(1 - 0.5^generations)))
}

#' Expected accuracy of individual classification
#'
#' Simulates animals of each class from the reference allele frequencies of
#' the diagnostic markers, adds genotyping error, and classifies them with
#' [classify_hybrids()]. The result shows how often each class can be called
#' correctly with these markers, and what it is mistaken for.
#'
#' @param markers Output of [diagnostic_markers()] (columns `freq1`,
#'   `freq2`).
#' @param n Animals simulated per class.
#' @param error Genotyping error rate used to simulate and classify.
#' @param alpha Passed to [classify_hybrids()].
#' @param seed Random seed.
#' @return A matrix of proportions: true classes in rows, called classes in
#'   columns. The diagonal is the expected accuracy.
#' @examples
#' ex <- hybridclock_example()
#' mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
#' round(classification_power(mk[1:15, ], n = 200), 2)
#' @export
classification_power <- function(markers, n = 1000, error = 0.005, alpha = 0.01, seed = 1) {
  with_seed(seed, {
    f1 <- markers$freq1; f2 <- markers$freq2; L <- length(f1)
    gamete <- function(a) matrix(stats::rbinom(n * L, 1, rep(a * f2 + (1 - a) * f1, each = n)), n)
    calls <- lapply(names(hybrid_classes), function(k) {
      ab <- hybrid_classes[[k]]
      g <- gamete(ab[1]) + gamete(ab[2])
      wrong <- matrix(stats::runif(n * L) < error, n)
      g[wrong] <- (g[wrong] + sample.int(2, sum(wrong), TRUE)) %% 3L
      sc <- data.frame(S = rowMeans(g) / 2, H = rowMeans(g == 1), n_markers = L)
      factor(classify_hybrids(sc, alpha = alpha, error = error, genotypes = g, markers = markers)$class,
             levels = c(names(hybrid_classes), "later"))
    })
    out <- t(vapply(calls, function(x) as.vector(table(x)) / n, numeric(length(hybrid_classes) + 1)))
    dimnames(out) <- list(true = names(hybrid_classes), called = c(names(hybrid_classes), "later"))
    out
  })
}
