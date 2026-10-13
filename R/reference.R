#' Simulate a reference set of hybrid histories
#'
#' Runs [simulate_hybrid_zone()] many times with settings drawn at random
#' from the ranges given, and stacks every generation of every run. The
#' ranges act as the prior of the estimator trained on the result, so they
#' should cover what is plausible for the study system. The breeding
#' population size matters most: past about five generations of a closed
#' hybrid swarm, drift is the only clock left, and its speed depends on
#' population size.
#'
#' Each run is one of three kinds of history, with probabilities `kinds`:
#' a closed swarm (no arrivals after founding), arrivals of one species
#' only (repeated backcrossing), or arrivals of both species (ongoing
#' hybridization). Arrivals of each species are, with equal probability,
#' constant, stop at a random generation, or start at a random generation,
#' so that the reference includes hybridization that has ended and
#' hybridization that began late.
#'
#' @param n_sims Number of simulated histories.
#' @param generations Generations per history.
#' @param breeders Range of breeding population sizes (drawn on a log scale).
#'   An effective size (Ne) range is usually the most relevant figure.
#' @param start_mix Range of the founding share of species 2.
#' @param arrival_max Largest share of breeders replaced by arrivals per
#'   generation.
#' @param kinds Probabilities of closed, one-species and two-species
#'   arrival histories.
#' @param sib_mating Range of sibling-mating probability; half the runs use 0.
#' @param pools Largest number of breeding pools; half the runs use 1.
#' @param n_markers,shared_freq,sample_size,missing Passed to
#'   [simulate_hybrid_zone()]. Match these to the observed data.
#' @param microsat Optional microsatellite markers (see [microsat_markers()]);
#'   observed scores must then come from [hybrid_scores_microsat()], and
#'   `n_markers` is set to the number of loci.
#' @param freq Optional allele frequencies for markers that are not fixed
#'   (see [simulate_hybrid_zone()]); observed scores must then come from
#'   `hybrid_scores(method = "likelihood")`.
#' @param assort Range of the probability of mating by ancestry (see
#'   [simulate_hybrid_zone()]); half the runs use 0.
#' @param het_dropout Range of the probability that a heterozygote is
#'   called homozygous; half the runs use 0 unless the range excludes 0. Give `c(0, 0)` for genotyping
#'   methods without allele dropout.
#' @param founder_impurity Range of the share of diagnostic alleles that
#'   founders carry from the other species (see [simulate_hybrid_zone()]),
#'   drawn for each history.
#' @param class_settings Settings of the classes behind the class shares in
#'   the summaries (see [class_settings()]): genotyping error, heterozygote
#'   dropout, linked markers, and reference allele frequencies for markers
#'   that are not fixed. They are stored with the reference and used again by
#'   [estimate_history()], so that the data are classified the same way.
#' @param founder_samples With F1 or backcross founders, also include samples
#'   of the founders themselves (generation 0), for experimental families of
#'   F1s or first backcrosses (see [simulate_hybrid_zone()]).
#' @param founders Founding designs to draw from (see
#'   [simulate_hybrid_zone()]): `"pure"` for natural contact, or `"F1"`,
#'   `"BC1_1"`, `"BC1_2"` for experimental or hatchery populations.
#' @param map Linkage: `NULL` for unlinked markers; a marker map (see
#'   [random_map()]) used by every history; or a list with ranges
#'   `chromosomes` and `length_cM`, from which each history draws its own
#'   random map, so that an unknown genome enters the reference as
#'   uncertainty. A fixed map should hold the observed markers' real
#'   positions, in the order of the diagnostic markers; then genotype
#'   switches along chromosomes are summarized too (see
#'   [sample_summaries()]).
#' @param map_scale With a fixed map, a range by which positions are
#'   multiplied in each history, for an uncertain recombination rate (for
#'   example, positions in Mb with `map_scale = c(0.5, 3)` for 0.5-3 cM per
#'   Mb).
#' @param cores Number of processes (forked; 1 on Windows).
#' @param seed Base random seed; run `i` uses `seed + i`.
#' @return A data.frame of all generations of all runs, with the run's
#'   settings (`run`, `breeders`, `start_mix`, `arrival1`, `arrival2` and
#'   their `timing1`, `timing2` and `switch1`, `switch2` generations, `kind`,
#'   `sib_mating`, `pools`), the true history and the summaries.
#'   Attribute `settings` records the arguments.
#' @examples
#' \donttest{
#' ref <- build_reference(n_sims = 20, generations = 10, n_markers = 200)
#' dim(ref)
#' }
#' @export
build_reference <- function(n_sims = 4000, generations = 50, breeders = c(10, 500),
                            start_mix = c(0.01, 0.99), arrival_max = 0.1,
                            kinds = c(closed = 0.2, one = 0.4, both = 0.4),
                            sib_mating = c(0, 0.9), pools = 4, n_markers = 1000,
                            shared_freq = NULL, sample_size = 20, missing = 0,
                            founders = "pure", map = NULL, map_scale = c(1, 1), assort = c(0, 0.8),
                            het_dropout = c(0, 0.3), freq = NULL, microsat = NULL, founder_impurity = c(0, 0),
                            founder_samples = FALSE, class_settings = NULL, cores = 1, seed = 1) {
  if (is.null(class_settings)) class_settings <- hybridclock::class_settings()
  if (!is.null(microsat)) n_markers <- length(microsat$loci)
  founders <- match.arg(founders, c("pure", "F1", "BC1_1", "BC1_2"), several.ok = TRUE)
  one <- function(i) with_seed(seed + i, {
    N <- round(exp(stats::runif(1, log(breeders[1]), log(breeders[2]))))
    mix <- stats::runif(1, start_mix[1], start_mix[2])
    kind <- sample(c("closed", "one", "both"), 1, prob = kinds)
    arr <- switch(kind, closed = c(0, 0),
                  one = { a <- c(0, 0); a[sample.int(2, 1)] <- stats::runif(1, 0, arrival_max); a },
                  both = exp(stats::runif(2, log(0.002), log(arrival_max))))
    timing <- ifelse(arr > 0, sample(c("constant", "stopped", "started"), 2, TRUE), "none")
    switch_at <- ifelse(arr > 0 & timing != "constant", 1L + sample.int(max(1L, generations - 2L), 2, TRUE), NA)
    sched <- vapply(1:2, function(k) arr[k] * switch(timing[k], none = rep(0, generations),
      constant = rep(1, generations), stopped = as.numeric(seq_len(generations) < switch_at[k]),
      started = as.numeric(seq_len(generations) >= switch_at[k])), numeric(generations))
    sib <- if (stats::runif(1) < 0.5) 0 else stats::runif(1, sib_mating[1], sib_mating[2])
    np <- if (pools < 2 || stats::runif(1) < 0.5) 1L else sample(2:pools, 1)
    cross <- stats::runif(1, 0, 0.2)
    fd <- if (length(founders) == 1) founders else sample(founders, 1)
    chrom <- NA; len <- NA; mp <- NULL
    if (is.data.frame(map)) {
      len <- stats::runif(1, map_scale[1], map_scale[length(map_scale)])
      mp <- data.frame(chromosome = map$chromosome, position = map$position * len)
    }
    if (is.list(map) && !is.data.frame(map)) {
      chrom <- round(stats::runif(1, map$chromosomes[1], map$chromosomes[length(map$chromosomes)]))
      len <- stats::runif(1, map$length_cM[1], map$length_cM[length(map$length_cM)])
      mp <- random_map(n_markers, chrom, len)
    }
    as_ <- if (assort[2] <= 0 || stats::runif(1) < 0.5) 0 else stats::runif(1, assort[1], assort[2])
    # dropout is absent in half the runs unless the range excludes zero
    dr <- if (het_dropout[2] <= 0 || (het_dropout[1] <= 0 && stats::runif(1) < 0.5)) 0 else
      stats::runif(1, het_dropout[1], het_dropout[2])
    imp <- if (founder_impurity[2] <= 0) 0 else stats::runif(1, founder_impurity[1], founder_impurity[2])
    s <- simulate_hybrid_zone(generations, N, mix, matrix(sched, ncol = 2), sib, np, cross, n_markers,
                              shared_freq, sample_size, missing, founders = fd, map = mp,
                              assort = as_, het_dropout = dr, junctions = is.data.frame(map), freq = freq,
                              microsat = microsat, founder_impurity = imp,
                              sample_founders = founder_samples, class_settings = class_settings)
    cbind(run = i, breeders = N, start_mix = mix, arrival1 = arr[1], arrival2 = arr[2],
          timing1 = match(timing[1], timings), timing2 = match(timing[2], timings),
          switch1 = switch_at[1], switch2 = switch_at[2],
          kind = match(kind, c("closed", "one", "both")), sib_mating = sib, pools = np,
          founders = match(fd, c("pure", "F1", "BC1_1", "BC1_2")), chromosomes = chrom, map_cM = len,
          assort = as_, het_dropout = dr, founder_impurity = imp, s)
  })
  res <- if (cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n_sims), one, mc.cores = cores)
  } else lapply(seq_len(n_sims), one)
  bad <- vapply(res, inherits, NA, "try-error")
  if (any(bad)) stop("Some simulations failed: ", as.character(res[[which(bad)[1]]]), call. = FALSE)
  out <- do.call(rbind, res)
  out$kind <- c("closed", "one", "both")[out$kind]
  out$timing1 <- timings[out$timing1]; out$timing2 <- timings[out$timing2]
  out$founders <- c("pure", "F1", "BC1_1", "BC1_2")[out$founders]
  attr(out, "settings") <- list(generations = generations, breeders = breeders, start_mix = start_mix,
                                arrival_max = arrival_max, kinds = kinds, sib_mating = sib_mating,
                                pools = pools, n_markers = n_markers, n_shared = length(shared_freq),
                                sample_size = sample_size, founders = founders, assort = assort,
                                het_dropout = het_dropout, founder_impurity = founder_impurity,
                                founder_samples = founder_samples, class_settings = class_settings,
                                map = if (is.data.frame(map)) "fixed" else map)
  out
}

timings <- c("none", "constant", "stopped", "started")
