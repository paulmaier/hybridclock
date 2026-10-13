#' Example data set
#'
#' Simulates a small data set: 20 animals from each of two pure reference
#' groups and two hybrid sites. Site A is a young hybrid swarm (4
#' generations, no arrivals); site B has been backcrossing into species 1 for
#' 15 generations. Animals carry 300 diagnostic SNPs, 40 SNPs that are not
#' quite fixed, and 150 shared SNPs, on 490 separate loci. The data are
#' simulated afresh with a fixed seed, so every call returns the same data.
#'
#' @return A list with `genotypes` (animals x SNPs, counts of the reference
#'   allele), `groups` (`"species1"`, `"species2"`, `"site_A"`, `"site_B"`),
#'   `locus` (locus of each SNP) and `truth` (the simulated histories of the
#'   two sites at the time of sampling).
#' @examples
#' ex <- hybridclock_example()
#' table(ex$groups)
#' ex$truth
#' @export
hybridclock_example <- function() {
  with_seed(2027, {
    n_diag <- 300; n_near <- 40; n_shared <- 150
    shared_freq <- stats::runif(n_shared, 0.1, 0.9)
    near1 <- stats::runif(n_near, 0, 0.03); near2 <- 1 - stats::runif(n_near, 0, 0.03)
    site <- function(gens, ...) {
      s <- simulate_hybrid_zone(generations = gens, n_markers = n_diag, shared_freq = shared_freq,
                                sample_size = 30, sample_at = gens, ...)
      list(G = attr(s, "samples")[[1]]$genotypes, truth = s[gens, c("generation", "hybrid_age_both",
            "hybrid_age_rare", "backcrosses", "arrival_rate")])
    }
    A <- site(4, breeders = 80, start_mix = 0.4)
    B <- site(15, breeders = 60, start_mix = 0.15, arrivals = c(0.06, 0), pools = 2)
    pure <- function(code, n) cbind(matrix(code, n, n_diag),
                                    matrix(stats::rbinom(n * n_shared, 2, rep(shared_freq, each = n)), n, n_shared))
    G <- rbind(pure(0L, 20), pure(2L, 20), A$G, B$G)
    groups <- rep(c("species1", "species2", "site_A", "site_B"), c(20, 20, nrow(A$G), nrow(B$G)))
    # nearly fixed SNPs follow each animal's ancestry share
    anc <- rowMeans(G[, seq_len(n_diag)]) / 2
    near <- vapply(seq_len(n_near), function(j) stats::rbinom(nrow(G), 2, (1 - anc) * near1[j] + anc * near2[j]), numeric(nrow(G)))
    G <- cbind(G[, seq_len(n_diag)], near, G[, n_diag + seq_len(n_shared)])
    # flip the allele counted at a third of SNPs, as real data do
    flip <- stats::runif(ncol(G)) < 1 / 3
    G[, flip] <- 2L - G[, flip]
    storage.mode(G) <- "integer"
    colnames(G) <- paste0("snp", seq_len(ncol(G)))
    rownames(G) <- paste0(groups, "_", stats::ave(seq_along(groups), groups, FUN = seq_along))
    list(genotypes = G, groups = groups, locus = paste0("locus", seq_len(ncol(G))),
         truth = rbind(cbind(site = "site_A", A$truth), cbind(site = "site_B", B$truth)))
  })
}
