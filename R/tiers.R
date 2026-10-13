# Two-tier answer: when most hybrids in a sample fit early classes (F1, F2,
# first and second backcrosses), their ages follow directly from the
# classes; otherwise the population clock is used.

# ages of each early class, after orientation (species 2 = rare): hybrid age
# over both species' lines, over the rare species' lines, and backcrosses
early_classes <- data.frame(
  class = c("F1", "F2", "BC1_common", "BC1_rare", "BC2_common", "BC2_rare"),
  hybrid_age_both = c(1, 2, 1.5, 1.5, 1.75, 1.75),
  hybrid_age_rare = c(1, 2, 2, 4 / 3, 3, 11 / 7),
  backcrosses = c(0, 0, 1, 1, 2, 2))

# class-based estimates from the class shares of sample_summaries(); rows of
# `x` are samples
class_ages <- function(x) {
  x <- as.data.frame(x)
  early <- as.matrix(x[, paste0("class_", early_classes$class), drop = FALSE])
  later <- x$class_later
  n_early <- rowSums(early)
  out <- data.frame(early_share = n_early / pmax(n_early + later, 1e-12),
                    f2_share = x$class_F2 / pmax(n_early, 1e-12),
                    f1_share = x$class_F1 / pmax(n_early, 1e-12))
  for (t in c("hybrid_age_both", "hybrid_age_rare", "backcrosses"))
    out[[t]] <- as.vector(early %*% early_classes[[t]]) / pmax(n_early, 1e-12)
  out[n_early == 0, -(1:3)] <- NA
  out
}

# F2-like animals cannot be told from later generations of a symmetric swarm
# one at a time, so the classes are trusted only when the early hybrids are
# mostly F1s and backcrosses
max_f2_share <- 0.5

# large error (absolute log ratio, `level` quantile) at which class-based
# ages are accurate enough to use even where the clock does better on
# simulated runs
class_tolerance <- 0.3

# smallest share of hybrids in early classes for a class-based answer
min_early_share <- 0.8

# share of F1s among early hybrids above which the F1 rule applies
min_f1_share <- 0.8

# smallest half-width (log scale) of class-based intervals
class_floor <- 0.1

# choose the share of early-class hybrids above which the class-based age is
# more accurate than the clock, on held-out runs; `spread_classes` is the
# `level` quantile of the class-based error in each bin, which sets the
# width of class-based intervals
choose_tier <- function(clock, ref, level = 0.9) {
  ca <- class_ages(ref)
  ok <- is.finite(ca$hybrid_age_both) & is.finite(ref$hybrid_age_both) & ca$f2_share <= max_f2_share
  if (sum(ok) < 200) return(list(cutoff = 0.9, table = NULL))
  pc <- stats::predict(clock$models$hybrid_age_both, prep_features(ref[ok, ], clock$features))
  y <- ref$hybrid_age_both[ok]
  err_clock <- abs(pc - log(y)); err_class <- abs(log(ca$hybrid_age_both[ok] / y))
  edges <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
  b <- cut(ca$early_share[ok], edges, include.lowest = TRUE)
  tab <- do.call(rbind, lapply(levels(b), function(l) {
    i <- which(b == l)
    data.frame(early_share = l, n = length(i), error_clock = if (length(i)) stats::median(err_clock[i]) else NA,
               error_classes = if (length(i)) stats::median(err_class[i]) else NA,
               spread_clock = if (length(i)) stats::quantile(err_clock[i], level, names = FALSE) else NA,
               spread_classes = if (length(i)) stats::quantile(err_class[i], level, names = FALSE) else NA)
  }))
  # classes are acceptable where their large errors (the `level` quantile)
  # are no worse than the clock's, or small in their own right (the clock's
  # edge on simulated runs may not carry over to real data that the
  # reference fits less well); judging by large errors keeps the classes out
  # of bins where a minority of older samples is badly misdated
  better <- with(tab, n >= 20 & (spread_classes <= spread_clock | spread_classes <= class_tolerance))
  # lowest bin from which the classes are acceptable in every higher bin with data
  cutoff <- Inf
  for (k in rev(seq_len(nrow(tab)))) {
    if (tab$n[k] < 20) next
    if (better[k]) cutoff <- edges[k] else break
  }
  # F1s stay recognizable with few markers (heterozygous nearly everywhere),
  # so samples made up mostly of F1s may be dated by their classes even where
  # the early share alone does not qualify; learned the same way
  f1 <- which(ca$f1_share[ok] >= min_f1_share & ca$early_share[ok] >= min_early_share)
  f1_spread <- if (length(f1) >= 20) stats::quantile(err_class[f1], level, names = FALSE) else Inf
  f1_ok <- f1_spread <= class_tolerance
  # class-based ages describe only the early hybrids, so they are never used
  # when many hybrids fit no early class
  list(cutoff = max(cutoff, min_early_share), f1_ok = f1_ok, f1_spread = f1_spread, table = tab)
}

# widen class-based intervals to at least the class-based error seen in the
# reference for samples with a similar early share
class_spread <- function(tier, early_share) {
  tab <- tier$table
  if (is.null(tab) || is.null(tab$spread_classes) || !is.finite(early_share)) return(class_floor)
  edges <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
  k <- findInterval(early_share, edges, rightmost.closed = TRUE, left.open = TRUE)
  q <- tab$spread_classes[max(1, k)]
  if (is.finite(q)) max(q, class_floor) else class_floor
}
