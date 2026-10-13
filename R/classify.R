#' Classify recent hybrids individually
#'
#' Tests each animal against the pure, first- and second-generation hybrid
#' classes and second backcrosses. Each class is a pair of parents, each
#' passing on species-2 alleles with a set probability: 0 (pure species 1),
#' 1 (pure species 2), 0.5 (F1 or F2 parent), 0.25 or 0.75 (first backcross
#' parent). For each animal the likelihood of its genotypes under every class
#' is compared with the best fit under any pair of parental ancestries. If no
#' class fits (likelihood-ratio test at `alpha`), the animal is called a
#' later-generation hybrid.
#'
#' With `genotypes` and `markers`, the likelihood uses each marker's allele
#' frequencies in the two reference groups, so markers that are nearly but
#' not perfectly fixed are handled correctly, as in NewHybrids (Anderson and
#' Thompson 2002) and HIest (Fitzpatrick 2012). A parent with species-2
#' ancestry `a` passes on the species-2 allele at marker `l` with probability
#' `a * freq2[l] + (1 - a) * freq1[l]`. The genotyping error rate can be
#' estimated from the data (`error = "estimate"`). Without them, markers are
#' treated as fixed and only the genotype counts implied by `S` and `H` are
#' used.
#'
#' Several third-generation classes are indistinguishable from earlier ones
#' at unlinked markers whatever the number of markers: F1 x F2 and F2 x F2
#' look like F2, and F2 x pure looks like a first backcross. Inbreeding
#' (related parents) lowers heterozygosity and moves animals off every class,
#' so they are called later-generation hybrids. Use [estimate_history()] for
#' populations of advanced hybrids. Markers must be diagnostic in the
#' populations that actually hybridized, not only in the reference samples:
#' a marker at which one source population carries the "wrong" allele makes
#' its hybrids look like no class.
#'
#' @param scores Output of [hybrid_scores()] (`S`, `H`, `n_markers`).
#' @param alpha Significance level of the fit test.
#' @param error Genotyping error rate, or `"estimate"` (with `genotypes`).
#' @param genotypes Optional genotypes at the diagnostic markers, counting
#'   species-2 alleles, as from [orient_genotypes()], one row per row of
#'   `scores`.
#' @param markers Output of [diagnostic_markers()] matching the columns of
#'   `genotypes` (columns `freq1` and `freq2` are used).
#' @param het_dropout Probability that a heterozygote is called homozygous
#'   (allele dropout, common in reduced-representation sequencing). It is
#'   not estimated, because heterozygosity lost to dropout cannot be told
#'   apart from heterozygosity lost to inbreeding one animal at a time; set it
#'   from replicate genotyping or known error rates.
#' @param effective_markers Optional number of effectively independent
#'   markers. Markers on the same chromosome are inherited in blocks, so they
#'   carry less independent information than their number suggests. With
#'   this, log-likelihoods are scaled by `effective_markers / n_markers`
#'   (a composite-likelihood adjustment). The best class is unchanged, but
#'   animals are no longer rejected from every class just because their
#'   markers come in blocks. See [effective_markers()].
#' @return `scores` with columns `class` (best class or `"later"`),
#'   `probability` (of the best class among the tested classes) and
#'   `fit_p` (p-value of the fit test). With `error = "estimate"`, attribute
#'   `error` holds the estimated rate.
#' @references Anderson EC, Thompson EA (2002) A model-based method for
#'   identifying species hybrids using multilocus genetic data. Genetics
#'   160:1217-1229.
#' @examples
#' ex <- hybridclock_example()
#' mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
#' sc <- hybrid_scores(ex$genotypes, mk)
#' table(ex$groups, classify_hybrids(sc)$class)
#' cl <- classify_hybrids(sc, error = "estimate", genotypes = orient_genotypes(ex$genotypes, mk),
#'                        markers = mk)
#' attr(cl, "error")
#' @export
classify_hybrids <- function(scores, alpha = 0.01, error = 0.005, genotypes = NULL, markers = NULL,
                             effective_markers = NULL, het_dropout = 0) {
  if (is.null(genotypes)) {
    if (identical(error, "estimate")) stop("Estimating the error rate needs `genotypes` and `markers`.", call. = FALSE)
    n <- scores$n_markers
    n1 <- round(n * scores$H)
    n2 <- pmax(0, round(n * (2 * scores$S - scores$H) / 2))
    n0 <- pmax(0, n - n1 - n2)
    G <- cbind(n0, n1, n2)
    loglik <- function(a, b, e) {
      p2 <- a * b; p0 <- (1 - a) * (1 - b); p1 <- 1 - p0 - p2
      p0 <- p0 + p1 * het_dropout / 2; p2 <- p2 + p1 * het_dropout / 2; p1 <- p1 * (1 - het_dropout)
      lp <- log((1 - e) * rbind(p0, p1, p2) + e / 3)
      G %*% lp
    }
  } else {
    stopifnot(nrow(genotypes) == nrow(scores), ncol(genotypes) == nrow(markers))
    I <- lapply(0:2, function(k) { m <- genotypes == k; m[is.na(m)] <- FALSE; m * 1 })
    f1 <- markers$freq1; f2 <- markers$freq2
    loglik <- function(a, b, e) {
      qa <- outer(f2, a) + outer(f1, 1 - a); qb <- outer(f2, b) + outer(f1, 1 - b)
      p2 <- qa * qb; p0 <- (1 - qa) * (1 - qb); p1 <- 1 - p0 - p2
      p0 <- p0 + p1 * het_dropout / 2; p2 <- p2 + p1 * het_dropout / 2; p1 <- p1 * (1 - het_dropout)
      err <- function(p) log((1 - e) * p + e * (1 - p) / 2)
      I[[1]] %*% err(p0) + I[[2]] %*% err(p1) + I[[3]] %*% err(p2)
    }
  }
  ca <- vapply(hybrid_classes, `[`, 0, 1); cb <- vapply(hybrid_classes, `[`, 0, 2)
  free_ll <- function(e) apply(loglik(free_grid$a, free_grid$b, e), 1, max)
  if (identical(error, "estimate")) {
    grid_e <- c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
    error <- grid_e[which.max(vapply(grid_e, function(e) sum(free_ll(e)), 0))]
  }
  ll <- loglik(ca, cb, error)
  colnames(ll) <- names(hybrid_classes)
  free <- free_ll(error)
  if (!is.null(effective_markers)) {
    w <- pmin(1, effective_markers / pmax(scores$n_markers, 1))
    ll <- ll * w; free <- free * w
  }
  best <- max.col(ll, ties.method = "first")
  best_ll <- ll[cbind(seq_len(nrow(ll)), best)]
  fit_p <- stats::pchisq(2 * pmax(free - best_ll, 0), df = 2, lower.tail = FALSE)
  prob <- exp(best_ll - apply(ll, 1, function(r) max(r) + log(sum(exp(r - max(r))))))
  scores$class <- ifelse(fit_p < alpha, "later", colnames(ll)[best])
  scores$probability <- prob
  scores$fit_p <- fit_p
  attr(scores, "error") <- error
  scores
}

free_grid <- local({
  g <- expand.grid(a = seq(0, 1, by = 0.02), b = seq(0, 1, by = 0.02))
  g[g$a <= g$b, ]
})

hybrid_classes <- list(
  pure1 = c(0, 0), pure2 = c(1, 1), F1 = c(0, 1), F2 = c(0.5, 0.5),
  BC1_species1 = c(0, 0.5), BC1_species2 = c(0.5, 1),
  BC2_species1 = c(0, 0.25), BC2_species2 = c(0.75, 1))

#' Effective number of independent markers
#'
#' Markers on one chromosome travel together through a meiosis, broken only
#' by crossovers. An early-generation hybrid's genotype is shaped by the
#' recombinant gamete(s) of its hybrid parent(s): one for a backcross, two
#' for an F2. The number of independently inherited pieces is about
#' `meioses * (chromosomes + map length in Morgans)`. The default of one
#' meiosis suits backcrosses, the classes most affected; in the bear
#' pedigree truth set it let every known backcross be classified correctly
#' while animals with old introgression were still called later-generation.
#'
#' @param chromosomes Number of chromosomes (haploid).
#' @param length_cM Total map length in centimorgans. If unknown, about 50-150
#'   cM per chromosome is typical.
#' @param meioses Number of recombinant gametes shaping the genotype (1 for
#'   a backcross, 2 for an F2).
#' @return A number to pass to [classify_hybrids()] as `effective_markers`.
#' @examples
#' effective_markers(chromosomes = 37, length_cM = 2500)
#' @export
effective_markers <- function(chromosomes, length_cM = 100 * chromosomes, meioses = 1) {
  meioses * (chromosomes + length_cM / 100)
}
