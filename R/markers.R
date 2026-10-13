#' Find diagnostic markers between two reference groups
#'
#' Keeps SNPs that are nearly fixed for different alleles in two reference
#' groups of (presumed) pure animals: at most `threshold` frequency in one
#' group and at least `1 - threshold` in the other. With `locus`, only the
#' first such SNP per locus is kept, so markers are effectively unlinked.
#'
#' @param geno Genotype matrix, animals in rows and SNPs in columns, counting
#'   copies of one allele (0, 1, 2; `NA` for missing).
#' @param groups Group label of each animal.
#' @param ref1,ref2 Labels of the two reference groups. `ref2` defines
#'   species 2, whose alleles the hybrid index counts.
#' @param threshold Largest frequency allowed for the minor allele in each
#'   reference group (0.01 keeps markers within 1% of fixation).
#' @param locus Optional locus name of each SNP.
#' @return A data.frame with one row per marker: `column` (index in `geno`),
#'   `name`, `flip` (`TRUE` when the counted allele belongs to species 1),
#'   and the counted allele's frequency in each reference group (`freq1`,
#'   `freq2`).
#' @examples
#' ex <- hybridclock_example()
#' mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
#' nrow(mk)
#' @export
diagnostic_markers <- function(geno, groups, ref1, ref2, threshold = 0.01, locus = NULL) {
  geno <- as.matrix(geno)
  stopifnot(nrow(geno) == length(groups), ref1 %in% groups, ref2 %in% groups)
  f1 <- colMeans(geno[groups == ref1, , drop = FALSE], na.rm = TRUE) / 2
  f2 <- colMeans(geno[groups == ref2, , drop = FALSE], na.rm = TRUE) / 2
  up <- f1 <= threshold & f2 >= 1 - threshold
  down <- f2 <= threshold & f1 >= 1 - threshold
  keep <- which((up | down) %in% TRUE)
  if (!is.null(locus)) keep <- keep[!duplicated(locus[keep])]
  if (!length(keep)) stop("No SNP passes the threshold; try a larger `threshold`.", call. = FALSE)
  nm <- colnames(geno); if (is.null(nm)) nm <- paste0("snp", seq_len(ncol(geno)))
  flip <- down[keep]
  data.frame(column = keep, name = nm[keep], flip = flip,
             freq1 = ifelse(flip, 1 - f1[keep], f1[keep]), freq2 = ifelse(flip, 1 - f2[keep], f2[keep]),
             row.names = NULL)
}

#' Find loci shared by both species
#'
#' Keeps SNPs that are common in both reference groups and differ little
#' between them. Inbreeding and linkage disequilibrium measured at these
#' loci help separate inbreeding from hybrid age (see [sample_summaries()]).
#'
#' @inheritParams diagnostic_markers
#' @param min_maf Smallest minor allele frequency required in each group.
#' @param max_diff Largest allele frequency difference allowed between groups.
#' @return A data.frame with `column`, `name` and the pooled allele
#'   frequency `freq`.
#' @examples
#' ex <- hybridclock_example()
#' nrow(shared_loci(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus))
#' @export
shared_loci <- function(geno, groups, ref1, ref2, min_maf = 0.05, max_diff = 0.1, locus = NULL) {
  geno <- as.matrix(geno)
  f1 <- colMeans(geno[groups == ref1, , drop = FALSE], na.rm = TRUE) / 2
  f2 <- colMeans(geno[groups == ref2, , drop = FALSE], na.rm = TRUE) / 2
  keep <- which((pmin(f1, 1 - f1) >= min_maf & pmin(f2, 1 - f2) >= min_maf & abs(f1 - f2) < max_diff) %in% TRUE)
  if (!is.null(locus)) keep <- keep[!duplicated(locus[keep])]
  nm <- colnames(geno); if (is.null(nm)) nm <- paste0("snp", seq_len(ncol(geno)))
  data.frame(column = keep, name = nm[keep], freq = (f1[keep] + f2[keep]) / 2, row.names = NULL)
}

#' Count species-2 alleles at diagnostic markers
#'
#' @param geno Genotype matrix as in [diagnostic_markers()].
#' @param markers Output of [diagnostic_markers()].
#' @return An integer matrix (animals x markers) counting species-2 alleles.
#' @examples
#' ex <- hybridclock_example()
#' mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
#' dim(orient_genotypes(ex$genotypes, mk))
#' @export
orient_genotypes <- function(geno, markers) {
  g <- as.matrix(geno)[, markers$column, drop = FALSE]
  g[, markers$flip] <- 2L - g[, markers$flip]
  storage.mode(g) <- "integer"
  g
}

#' Hybrid index and interspecies heterozygosity
#'
#' For each animal, the hybrid index S is the share of species-2 alleles at
#' the diagnostic markers, and the interspecies heterozygosity H is the share
#' of markers carrying one allele from each species. With markers that are
#' fixed or nearly fixed, these simple counts match the maximum-likelihood
#' estimates of Fitzpatrick (2012) closely.
#'
#' @inheritParams orient_genotypes
#' @param method `"count"` (simple counts; for fixed or nearly fixed
#'   markers) or `"likelihood"` (maximum-likelihood estimates that use each
#'   marker's allele frequencies in the two reference groups, for markers that
#'   are not fixed; Fitzpatrick 2012).
#' @return A data.frame with `S`, `H` and `n_markers` (markers genotyped).
#' @references Fitzpatrick BM (2012) Estimating ancestry and heterozygosity of
#'   hybrids using molecular markers. BMC Evolutionary Biology 12:131.
#' @examples
#' ex <- hybridclock_example()
#' mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
#' head(hybrid_scores(ex$genotypes, mk))
#' @export
hybrid_scores <- function(geno, markers, method = c("count", "likelihood")) {
  method <- match.arg(method)
  g <- orient_genotypes(geno, markers)
  if (method == "likelihood") return(ml_scores(g, markers$freq1, markers$freq2, rownames(g)))
  n <- rowSums(!is.na(g))
  data.frame(S = rowSums(g, na.rm = TRUE) / (2 * n), H = rowSums(g == 1L, na.rm = TRUE) / n,
             n_markers = n, row.names = rownames(g))
}

# maximum-likelihood hybrid index and interspecies heterozygosity for markers
# that are not fixed (Fitzpatrick 2012): at each marker an animal carries one
# allele from each species with probability H, two species-2 alleles with
# probability S - H/2, and two species-1 alleles with probability 1 - S - H/2;
# alleles are then drawn from each species' frequencies
ml_scores <- function(g, f1, f2, ids = NULL) {
  f1 <- pmin(pmax(f1, 0.001), 0.999); f2 <- pmin(pmax(f2, 0.001), 0.999)
  # genotype probabilities for each origin: rows markers, columns genotype 0/1/2
  pg <- function(pa, pb) cbind((1 - pa) * (1 - pb), pa * (1 - pb) + pb * (1 - pa), pa * pb)
  o11 <- pg(f1, f1); o12 <- pg(f1, f2); o22 <- pg(f2, f2)
  grid <- ml_grid
  w11 <- 1 - grid$S - grid$H / 2; w12 <- grid$H; w22 <- grid$S - grid$H / 2
  I <- lapply(0:2, function(k) { m <- g == k; m[is.na(m)] <- FALSE; m * 1 })
  ll <- 0
  for (k in 1:3) {
    P <- outer(o11[, k], w11) + outer(o12[, k], w12) + outer(o22[, k], w22)   # markers x grid
    ll <- ll + I[[k]] %*% log(pmax(P, 1e-12))
  }
  best <- max.col(ll, ties.method = "first")
  data.frame(S = grid$S[best], H = grid$H[best], n_markers = rowSums(!is.na(g)), row.names = ids)
}

ml_grid <- local({
  g <- expand.grid(S = seq(0, 1, by = 0.02), H = seq(0, 1, by = 0.02))
  g[g$H <= 2 * pmin(g$S, 1 - g$S) + 1e-9, ]
})

#' Information carried by a set of markers
#'
#' The sum over markers of the squared allele-frequency difference between
#' the two reference groups. A fixed marker contributes 1, so the total is
#' roughly the number of fixed markers the set is worth. In simulations,
#' accuracy of hybrid ages fell gradually with it: typical errors were about
#' 20% at 300, 26% at 126 and 36-39% at 7-12, with honest intervals
#' throughout.
#'
#' @param markers Output of [diagnostic_markers()] (columns `freq1`, `freq2`).
#' @return A number.
#' @examples
#' ex <- hybridclock_example()
#' marker_information(diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2"))
#' @export
marker_information <- function(markers) sum((markers$freq2 - markers$freq1)^2)

#' Check that markers are diagnostic in the parents that actually bred
#'
#' Markers chosen in reference samples can still be shared by the animals
#' that founded a hybrid population, for example when a release strain or
#' cross parent differs from the reference sample. First-generation hybrids
#' expose such markers: at a truly diagnostic marker every F1 is
#' heterozygous. This function finds animals that look like F1s (hybrid
#' index near 0.5 and high interspecies heterozygosity) and reports, for each
#' marker, how often they are heterozygous there.
#'
#' @inheritParams orient_genotypes
#' @param scores Output of [hybrid_scores()] for the same animals.
#' @param S_range Range of hybrid index for likely F1s.
#' @param min_H Smallest interspecies heterozygosity for likely F1s.
#' @param flag Markers heterozygous in fewer than this share of likely F1s
#'   are flagged.
#' @return A data.frame with one row per marker: `name`, `n_f1` (likely F1s
#'   genotyped), `het_f1` (their share of heterozygotes) and `flagged`.
#'   Attribute `n_f1` gives the number of likely F1s; with fewer than about
#'   5, the check says little.
#' @examples
#' ex <- hybridclock_example()
#' mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
#' chk <- check_markers(ex$genotypes, mk, hybrid_scores(ex$genotypes, mk))
#' attr(chk, "n_f1")
#' @export
check_markers <- function(geno, markers, scores, S_range = c(0.4, 0.6), min_H = 0.7, flag = 0.7) {
  g <- orient_genotypes(geno, markers)
  f1 <- which(scores$S >= S_range[1] & scores$S <= S_range[2] & scores$H >= min_H)
  het <- if (length(f1)) colMeans(g[f1, , drop = FALSE] == 1L, na.rm = TRUE) else rep(NA_real_, ncol(g))
  n <- if (length(f1)) colSums(!is.na(g[f1, , drop = FALSE])) else rep(0L, ncol(g))
  out <- data.frame(name = markers$name, n_f1 = n, het_f1 = het, flagged = n >= 5 & het < flag, row.names = NULL)
  attr(out, "n_f1") <- length(f1)
  out
}
