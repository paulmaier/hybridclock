#' Microsatellite markers
#'
#' Allele frequencies of multi-allelic markers (such as microsatellites) in
#' two reference groups, for [hybrid_scores_microsat()] and for simulating
#' reference histories (`microsat` in [simulate_hybrid_zone()] and
#' [build_reference()]). Microsatellites rarely have alleles fixed between
#' species, so each animal's hybrid index and interspecies heterozygosity
#' are estimated by likelihood from these frequencies.
#'
#' @param alleles Integer matrix of allele sizes, two columns per locus
#'   (named `locus.a1`, `locus.a2`), animals in rows; `NA` for missing.
#' @param groups Group label of each animal.
#' @param ref1,ref2 Labels of the two reference groups.
#' @param prior Pseudo-count added to every allele seen in either group, so
#'   that an allele absent from one reference sample is rare, not
#'   impossible, in that species.
#' @return A list with `loci` (locus names), `alleles` (per locus, the allele
#'   sizes), and `freq1`, `freq2` (per locus, their frequencies in each
#'   reference group).
#' @examples
#' set.seed(1)
#' a <- cbind(L1.a1 = c(rep(100, 10), rep(104, 10)), L1.a2 = c(rep(100, 10), rep(104, 10)))
#' ms <- microsat_markers(a, rep(c("sp1", "sp2"), each = 10), "sp1", "sp2")
#' ms$freq1
#' @export
microsat_markers <- function(alleles, groups, ref1, ref2, prior = 0.5) {
  alleles <- as.matrix(alleles)
  loci <- unique(sub("\\.a[12]$", "", colnames(alleles)))
  out <- list(loci = loci, alleles = list(), freq1 = list(), freq2 = list())
  for (l in loci) {
    a <- alleles[, paste0(l, c(".a1", ".a2")), drop = FALSE]
    all_al <- sort(unique(as.vector(a[groups %in% c(ref1, ref2), ])))
    all_al <- all_al[!is.na(all_al)]
    f <- function(g) {
      x <- as.vector(a[groups == g, ]); x <- x[!is.na(x)]
      cnt <- tabulate(match(x, all_al), length(all_al)) + prior
      cnt / sum(cnt)
    }
    out$alleles[[l]] <- all_al; out$freq1[[l]] <- f(ref1); out$freq2[[l]] <- f(ref2)
  }
  out
}

#' Hybrid index and heterozygosity from microsatellites
#'
#' Maximum-likelihood hybrid index (S) and interspecies heterozygosity (H)
#' for each animal from multi-allelic markers, following Fitzpatrick (2012):
#' at each locus an animal carries two species-1 copies, one of each, or two
#' species-2 copies, with probabilities `1 - S - H/2`, `H` and `S - H/2`, and
#' each copy's allele is drawn from that species' frequencies. Alleles not
#' seen in either reference are treated as missing.
#'
#' @param alleles Allele matrix as in [microsat_markers()].
#' @param ms Output of [microsat_markers()].
#' @return A data.frame with `S`, `H` and `n_markers` (loci genotyped).
#' @examples
#' set.seed(1)
#' a <- cbind(L1.a1 = c(rep(100, 10), rep(104, 10)), L1.a2 = c(rep(100, 10), rep(104, 10)))
#' g <- rep(c("sp1", "sp2"), each = 10)
#' ms <- microsat_markers(a, g, "sp1", "sp2")
#' hybrid_scores_microsat(rbind(a, c(100, 104)), ms)
#' @export
hybrid_scores_microsat <- function(alleles, ms) {
  alleles <- as.matrix(alleles)
  n <- nrow(alleles); grid <- ml_grid
  w11 <- 1 - grid$S - grid$H / 2; w12 <- grid$H; w22 <- grid$S - grid$H / 2
  ll <- matrix(0, n, nrow(grid)); n_loci <- numeric(n)
  for (l in ms$loci) {
    x <- match(alleles[, paste0(l, ".a1")], ms$alleles[[l]])
    y <- match(alleles[, paste0(l, ".a2")], ms$alleles[[l]])
    ok <- !is.na(x) & !is.na(y)
    if (!any(ok)) next
    p1 <- ms$freq1[[l]]; p2 <- ms$freq2[[l]]
    xs <- x[ok]; ys <- y[ok]; het <- xs != ys
    k <- ifelse(het, 2, 1)
    p11 <- k * p1[xs] * p1[ys]
    p22 <- k * p2[xs] * p2[ys]
    p12 <- ifelse(het, p1[xs] * p2[ys] + p1[ys] * p2[xs], p1[xs] * p2[xs])
    ll[ok, ] <- ll[ok, ] + log(pmax(outer(p11, w11) + outer(p12, w12) + outer(p22, w22), 1e-300))
    n_loci <- n_loci + ok
  }
  best <- max.col(ll, ties.method = "first")
  data.frame(S = grid$S[best], H = grid$H[best], n_markers = n_loci, row.names = rownames(alleles))
}

# draw observed microsatellite alleles from each copy's species of origin
# (haplotype columns hold 0 = species 1, 1 = species 2), with missing loci
observe_microsat <- function(H1, H2, ms, missing) {
  n <- nrow(H1); L <- length(ms$loci)
  out <- matrix(NA_integer_, n, 2 * L, dimnames = list(NULL, as.vector(rbind(paste0(ms$loci, ".a1"), paste0(ms$loci, ".a2")))))
  draw <- function(origin, l) {
    a <- ms$alleles[[l]]
    ifelse(origin == 1, a[sample.int(length(a), n, TRUE, ms$freq2[[l]])], a[sample.int(length(a), n, TRUE, ms$freq1[[l]])])
  }
  for (j in seq_len(L)) {
    out[, 2 * j - 1] <- draw(H1[, j], ms$loci[j]); out[, 2 * j] <- draw(H2[, j], ms$loci[j])
  }
  rates <- if (length(missing) == 1) rep(missing, n) else sample(missing, n, TRUE)
  for (i in which(rates > 0)) {
    drop <- which(stats::runif(L) < rates[i])
    if (length(drop)) out[i, c(2 * drop - 1, 2 * drop)] <- NA
  }
  out
}

#' Classify recent hybrids from microsatellites
#'
#' The multi-allelic version of [classify_hybrids()]: a parent with
#' species-2 ancestry `a` passes on allele `x` with probability
#' `a * freq2(x) + (1 - a) * freq1(x)`. Each animal is tested against the
#' pure, F1, F2, first- and second-backcross classes and against the best
#' fit under any pair of parental ancestries; if no class fits
#' (likelihood-ratio test at `alpha`), it is called `"later"`.
#'
#' @param alleles Allele matrix as in [microsat_markers()].
#' @param ms Output of [microsat_markers()].
#' @param alpha Significance level of the fit test. Many later-generation
#'   animals (for example backcrosses of backcrosses) have genotypes close to
#'   an early class and cannot be told apart one at a time. A larger `alpha`
#'   calls more of them "later" at the cost of rejecting some true early
#'   hybrids: in a bighorn sheep pedigree with 187 microsatellites, `alpha =
#'   0.01` kept 91% of early calls right and called 20% of later-generation
#'   sheep "later"; `alpha = 0.1` gave 85% and 44%. Population ages from
#'   [estimate_history()] do not depend on this choice.
#' @param error Probability that a called allele is replaced by a random
#'   allele of that locus (genotyping error, null alleles).
#' @return A data.frame with `class`, `probability` and `fit_p`, one row per
#'   animal.
#' @examples
#' set.seed(1)
#' a <- cbind(L1.a1 = c(rep(100, 10), rep(104, 10)), L1.a2 = c(rep(100, 10), rep(104, 10)))
#' g <- rep(c("sp1", "sp2"), each = 10)
#' ms <- microsat_markers(a, g, "sp1", "sp2")
#' classify_hybrids_microsat(rbind(a, c(100, 104)), ms)
#' @export
classify_hybrids_microsat <- function(alleles, ms, alpha = 0.01, error = 0.01) {
  alleles <- as.matrix(alleles); n <- nrow(alleles)
  ca <- vapply(hybrid_classes, `[`, 0, 1); cb <- vapply(hybrid_classes, `[`, 0, 2)
  pa <- c(ca, free_grid$a); pb <- c(cb, free_grid$b)
  ll <- matrix(0, n, length(pa))
  for (l in ms$loci) {
    x <- match(alleles[, paste0(l, ".a1")], ms$alleles[[l]])
    y <- match(alleles[, paste0(l, ".a2")], ms$alleles[[l]])
    ok <- which(!is.na(x) & !is.na(y))
    if (!length(ok)) next
    k <- length(ms$alleles[[l]])
    f1 <- (1 - error) * ms$freq1[[l]] + error / k; f2 <- (1 - error) * ms$freq2[[l]] + error / k
    xs <- x[ok]; ys <- y[ok]
    qa_x <- outer(f2[xs], pa) + outer(f1[xs], 1 - pa); qb_y <- outer(f2[ys], pb) + outer(f1[ys], 1 - pb)
    qa_y <- outer(f2[ys], pa) + outer(f1[ys], 1 - pa); qb_x <- outer(f2[xs], pb) + outer(f1[xs], 1 - pb)
    p <- qa_x * qb_x
    het <- which(xs != ys)
    if (length(het)) p[het, ] <- (qa_x * qb_y + qa_y * qb_x)[het, , drop = FALSE]
    ll[ok, ] <- ll[ok, ] + log(pmax(p, 1e-300))
  }
  nc <- length(ca)
  cls <- ll[, seq_len(nc), drop = FALSE]; colnames(cls) <- names(hybrid_classes)
  free <- apply(ll[, -seq_len(nc), drop = FALSE], 1, max)
  best <- max.col(cls, ties.method = "first")
  best_ll <- cls[cbind(seq_len(n), best)]
  fit_p <- stats::pchisq(2 * pmax(free - best_ll, 0), df = 2, lower.tail = FALSE)
  prob <- exp(best_ll - apply(cls, 1, function(r) max(r) + log(sum(exp(r - max(r))))))
  data.frame(class = ifelse(fit_p < alpha, "later", colnames(cls)[best]), probability = prob, fit_p = fit_p,
             row.names = rownames(alleles))
}
