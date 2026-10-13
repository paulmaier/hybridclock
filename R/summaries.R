#' Summarize a population sample of hybrid scores
#'
#' Turns the hybrid index (S) and interspecies heterozygosity (H) of a sample
#' of animals, and optionally their genotypes, into the summary statistics
#' that [train_clock()] learns from. The same function is used for simulated
#' and observed samples, so the two are always comparable.
#'
#' The sample is first oriented so that the rarer species among its hybrids
#' counts as species 2 (S below 0.5). The summaries then describe:
#' * the spread of animals over the S-H triangle (means, percentiles, and
#'   the share of animals in each cell of a 10 x 6 grid);
#' * heterozygosity relative to its maximum for each animal's ancestry
#'   (`2 * min(S, 1 - S)`) and to the minimum possible for unrelated parents
#'   (`2 * S * (1 - S)`). Values below the latter point to inbreeding;
#' * with `diag`, how unevenly the rarer species' alleles are spread across
#'   markers (drift);
#' * with `shared`, inbreeding (`shared_fis`) and excess linkage
#'   disequilibrium (`shared_ld`) at loci that do not differ between species;
#' * with `n_markers`, the share of animals in each class of
#'   [classify_hybrids()] (pure, F1, F2, first and second backcrosses toward
#'   the common or the rare species, and later generations).
#'
#' @param S,H Hybrid index and interspecies heterozygosity of each animal, as
#'   from [hybrid_scores()].
#' @param diag Optional genotype matrix (animals x diagnostic markers) counting
#'   species-2 alleles, as from [orient_genotypes()]. Missing values allowed.
#' @param shared Optional genotype matrix (animals x shared loci) of allele
#'   counts. Missing values allowed.
#' @param n_markers Optional number of markers genotyped in each animal, as
#'   from [hybrid_scores()].
#' @param class_settings Settings of the classes behind the class shares
#'   (see [class_settings()]); the same settings must be used for the
#'   reference and the data.
#' @param positions Optional marker positions matching the columns of
#'   `diag`: a data.frame with `chromosome` and `position`. With them, the
#'   rate at which each hybrid's genotype switches between neighbouring
#'   markers along a chromosome is summarized (`switch_rate`,
#'   `switch_rate_sd`). Switches are ancestry junctions plus genotyping
#'   noise; they accumulate with every generation of recombination (Janzen
#'   et al. 2018), so they date early generations even when animals look
#'   alike in S and H.
#' @return A named numeric vector.
#' @examples
#' set.seed(1)
#' S <- c(0, 0, 0.12, 0.25, 0.3, 0.06, 0.5)
#' H <- c(0, 0, 0.22, 0.48, 0.41, 0.12, 1)
#' round(sample_summaries(S, H)[1:12], 3)
#' @export
sample_summaries <- function(S, H, diag = NULL, shared = NULL, n_markers = NULL, positions = NULL, class_settings = NULL) {
  stopifnot(length(S) == length(H), length(S) >= 2)
  n <- length(S)
  pure1 <- S < 0.005 & H < 0.005
  pure2 <- S > 0.995 & H < 0.005
  hy <- !(pure1 | pure2)
  if (any(hy) && mean(S[hy]) > 0.5) {
    S <- 1 - S; tmp <- pure1; pure1 <- pure2; pure2 <- tmp
    if (!is.null(diag)) diag <- 2L - diag
  }
  Sh <- S[hy]; Hh <- H[hy]
  if (any(hy)) {
    max_ratio <- pmin(Hh / pmax(2 * pmin(Sh, 1 - Sh), 1e-3), 1.5)
    floor_ratio <- Hh / pmax(2 * Sh * (1 - Sh), 1e-3)
  } else max_ratio <- floor_ratio <- NA
  cells <- table(cut(S, s_breaks, include.lowest = TRUE),
                 cut(H, h_breaks, include.lowest = TRUE, right = FALSE)) / n
  spread <- NA
  if (!is.null(diag)) {
    p <- colMeans(diag, na.rm = TRUE) / 2; pb <- mean(p, na.rm = TRUE)
    spread <- if (is.finite(pb) && pb > 0 && pb < 1) stats::var(p, na.rm = TRUE) / (pb * (1 - pb)) else 0
  }
  fis <- ld <- NA
  if (!is.null(shared)) {
    p <- colMeans(shared, na.rm = TRUE) / 2; he <- 2 * p * (1 - p); ok <- which(he > 0.05)
    if (length(ok) > 5) {
      fis <- 1 - sum(colMeans(shared[, ok, drop = FALSE] == 1L, na.rm = TRUE)) / sum(he[ok])
      r <- suppressWarnings(stats::cor(shared[, ok, drop = FALSE], use = "pairwise.complete.obs"))
      ld <- mean(r[upper.tri(r)]^2, na.rm = TRUE) - 1 / n
    }
  }
  sw <- c(switch_rate = NA_real_, switch_rate_sd = NA_real_)
  if (!is.null(positions) && !is.null(diag) && any(hy)) {
    o <- order(positions$chromosome, positions$position)
    D <- diag[hy, o, drop = FALSE]; chr <- positions$chromosome[o]; M <- ncol(D)
    if (M > 1) {
      A <- D[, -M, drop = FALSE]; B <- D[, -1, drop = FALSE]
      same <- matrix(chr[-M] == chr[-1], nrow(A), M - 1, byrow = TRUE)
      valid <- !is.na(A) & !is.na(B) & same
      rate <- rowSums((A != B) & valid, na.rm = TRUE) / pmax(rowSums(valid), 1)
      sw[] <- c(mean(rate), if (length(rate) > 1) stats::sd(rate) else 0)
    }
  }
  cls <- stats::setNames(rep(NA_real_, length(class_features)), class_features)
  if (!is.null(n_markers)) {
    k <- classify_counts(class_counts(S, H, n_markers, diag, class_settings), class_settings)
    cls[] <- vapply(names(class_map), function(c) mean(k == c), 0)
  }
  q <- function(x, p) stats::quantile(x, p, names = FALSE)
  v <- c(n_hybrids = sum(hy), share_pure_common = mean(pure1), share_pure_rare = mean(pure2),
         mean_S = mean(S), sd_S = stats::sd(S), q10_S = q(S, .1), q50_S = q(S, .5), q90_S = q(S, .9),
         mean_H = mean(H), sd_H = stats::sd(H), q10_H = q(H, .1), q50_H = q(H, .5), q90_H = q(H, .9),
         max_H = max(H),
         hyb_mean_S = if (any(hy)) mean(Sh) else NA, hyb_sd_S = if (sum(hy) > 1) stats::sd(Sh) else 0,
         hyb_mean_H = if (any(hy)) mean(Hh) else NA, hyb_sd_H = if (sum(hy) > 1) stats::sd(Hh) else 0,
         het_max_ratio = mean(max_ratio), het_max_ratio_sd = if (sum(hy) > 1) stats::sd(max_ratio) else 0,
         het_max_ratio_min = if (any(hy)) min(max_ratio) else NA,
         het_max_ratio_max = if (any(hy)) max(max_ratio) else NA,
         het_floor_ratio = mean(floor_ratio),
         het_floor_ratio_sd = if (sum(hy) > 1) stats::sd(floor_ratio) else 0,
         below_floor = if (any(hy)) mean(floor_ratio < 0.9) else NA,
         cor_SH = if (sum(hy) > 2 && stats::sd(Sh) > 0 && stats::sd(Hh) > 0) stats::cor(Sh, Hh) else 0,
         marker_spread = spread, shared_fis = fis, shared_ld = ld)
  cell_names <- paste0("cell_S", rep(seq_len(length(s_breaks) - 1), length(h_breaks) - 1),
                       "_H", rep(seq_len(length(h_breaks) - 1), each = length(s_breaks) - 1))
  c(v, stats::setNames(as.vector(cells), cell_names), cls, sw)
}

# after orientation species 1 is the common species and species 2 the rare one
class_map <- c(pure1 = "pure_common", pure2 = "pure_rare", F1 = "F1", F2 = "F2",
               BC1_species1 = "BC1_common", BC1_species2 = "BC1_rare",
               BC2_species1 = "BC2_common", BC2_species2 = "BC2_rare", later = "later")
class_features <- paste0("class_", class_map)

s_breaks <- c(0, .005, .05, .15, .3, .5, .7, .85, .95, .995, 1)
h_breaks <- c(0, .005, .1, .25, .4, .6, 1.01)
