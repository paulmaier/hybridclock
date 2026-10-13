#' Simulate a hybrid population forward in time
#'
#' Simulates a population founded by two species that interbreed for a
#' number of generations, with optional ongoing arrival of pure animals of
#' either species, mating between siblings, and separate breeding pools.
#' Each generation, a sample of offspring (for example tadpoles) is drawn
#' from that generation's matings, scored at unlinked diagnostic markers,
#' and summarized with [sample_summaries()]. The true history of every
#' sampled animal is tracked alongside, so the output can be used to train
#' and test estimators.
#'
#' Genotypes count species-2 alleles: 0 is a pure species-1 animal and 2 a
#' pure species-2 animal at every diagnostic marker. Shared loci (optional)
#' are polymorphic in both species at the same frequencies and are used to
#' measure inbreeding and linkage disequilibrium.
#'
#' Two hybrid ages and a backcross count are tracked for each hybrid:
#' * `hybrid_age_both` (hybrid age, both species): every line of the pedigree
#'   is followed back to its first pure ancestor, of either species, and the
#'   number of generations is averaged over lines. An F1 is 1, a first
#'   backcross 1.5, an F2 2. Recent backcrossing keeps it low.
#' * `hybrid_age_rare` (hybrid age, rare species): only the lines leading
#'   back to pure animals of the rarer species in the sample, weighted by how
#'   much ancestry each contributes: how long ago, on average, the rarer
#'   species' genes entered.
#' * `backcrosses`: the number of matings with a pure animal since the first
#'   cross (averaged over the two parents when both are hybrids).
#'
#' @param generations Number of generations to simulate.
#' @param breeders Number of breeding adults per generation.
#' @param start_mix Share of species 2 among the founders (pure founders
#'   only).
#' @param map Optional marker map for linkage: a data.frame with columns
#'   `chromosome` and `position` (cM), one row per diagnostic marker, as from
#'   [random_map()]. Markers are simulated in map order, and recombination
#'   between neighbours on a chromosome follows Haldane's map function.
#'   `NULL` (the default) treats all markers as unlinked. Shared loci are
#'   always unlinked.
#' @param freq Optional allele frequencies for markers that are not fixed:
#'   a list or data.frame with `freq1` and `freq2`, the frequency of the
#'   counted allele at each diagnostic marker in species 1 and species 2 (as
#'   from [diagnostic_markers()]). Pure animals then draw their alleles from
#'   these frequencies, and S and H are maximum-likelihood estimates (see
#'   [hybrid_scores()]).
#' @param founder_impurity Share of diagnostic alleles that founders carry
#'   from the other species, for founders that were not quite pure (for
#'   example, cross parents with some old introgression). Arrivals are
#'   always pure.
#' @param class_settings Settings of the classes in the summaries (see
#'   [class_settings()]).
#' @param sample_founders With F1 or backcross founders, also sample the
#'   founders themselves (generation 0), so that the output covers families
#'   of F1s or first backcrosses as well as their descendants.
#' @param microsat Optional output of [microsat_markers()]: the diagnostic
#'   markers are then microsatellites. Each copy's species of origin is
#'   simulated, observed alleles are drawn from that species' frequencies, and
#'   S and H are estimated with [hybrid_scores_microsat()]. `n_markers` is set
#'   to the number of loci.
#' @param junctions Also summarize genotype switches along chromosomes
#'   (requires `map`; see [sample_summaries()]). Use this only when the
#'   observed data have marker positions too.
#' @param founders Who founds the population: `"pure"` animals of both
#'   species (in proportions `start_mix`), or all `"F1"`, or all first
#'   backcrosses to species 1 (`"BC1_1"`) or species 2 (`"BC1_2"`), as in
#'   experimental and hatchery populations.
#' @param arrivals Share of breeders replaced each generation by newly
#'   arrived pure animals of species 1 and species 2: a length-2 vector for
#'   constant rates, or a matrix with one row per generation and two columns
#'   for rates that change over time (arrivals that start or stop).
#' @param sib_mating Probability that an animal mates with a sibling (shared
#'   mother) when one is available.
#' @param pools Number of breeding pools. Animals mate within their pool
#'   except with probability `pool_crossing`.
#' @param pool_crossing Probability that a mating crosses pools.
#' @param assort Probability that a mating is assortative by ancestry: the
#'   mate is the closest in ancestry among five random candidates (for
#'   example through mate choice or different breeding times).
#' @param het_dropout Probability that a heterozygous genotype is called as
#'   homozygous (either allele), as with allele dropout in reduced-
#'   representation sequencing. Applied to the observed genotypes.
#' @param n_markers Number of diagnostic markers.
#' @param shared_freq Allele frequencies of shared loci, or `NULL` for none.
#' @param sample_size Number of offspring sampled each generation.
#' @param missing Per-animal missing-data rates to draw from (a single value
#'   or a vector of observed rates).
#' @param sample_at Generations at which to return the sampled animals'
#'   genotypes (attribute `samples`).
#' @param seed Optional random seed. The global random state is restored
#'   afterwards.
#' @return A data.frame with one row per generation: `generation`, the true
#'   history of the sampled hybrids (`hybrid_age_both`, `hybrid_age_rare`,
#'   `backcrosses`, `recent_share` for the share that are F1, F2 or first
#'   backcrosses, `rare_species`, the current `arrival_rate` of that
#'   species (its mean over the three generations before sampling) and
#'   `arriving` (whether that rate is above zero)), and the summaries from
#'   [sample_summaries()].
#' @examples
#' sim <- simulate_hybrid_zone(generations = 8, breeders = 60, start_mix = 0.2,
#'                             arrivals = c(0.05, 0), n_markers = 300, seed = 1)
#' sim[, c("generation", "hybrid_age_both", "hybrid_age_rare", "backcrosses", "mean_S", "mean_H")]
#' @export
simulate_hybrid_zone <- function(generations = 20, breeders = 100, start_mix = 0.5,
                                 arrivals = c(0, 0), sib_mating = 0, pools = 1,
                                 pool_crossing = 0.1, n_markers = 1000, shared_freq = NULL,
                                 sample_size = 20, missing = 0, sample_at = NULL, seed = NULL,
                                 founders = c("pure", "F1", "BC1_1", "BC1_2"), map = NULL,
                                 assort = 0, het_dropout = 0, junctions = FALSE, freq = NULL,
                                 microsat = NULL, founder_impurity = 0, sample_founders = FALSE,
                                 class_settings = NULL) {
  founders <- match.arg(founders)
  if (!is.null(seed)) return(with_seed(seed, simulate_hybrid_zone(
    generations, breeders, start_mix, arrivals, sib_mating, pools, pool_crossing,
    n_markers, shared_freq, sample_size, missing, sample_at, seed = NULL, founders = founders, map = map,
    assort = assort, het_dropout = het_dropout, junctions = junctions, freq = freq, microsat = microsat,
    founder_impurity = founder_impurity, sample_founders = sample_founders,
    class_settings = class_settings)))
  if (!is.null(microsat)) n_markers <- length(microsat$loci)
  arr <- if (is.matrix(arrivals)) arrivals else matrix(arrivals, generations, 2, byrow = TRUE)
  stopifnot(breeders >= 4, nrow(arr) >= generations, ncol(arr) == 2, all(arr >= 0), all(rowSums(arr) <= 1),
            start_mix > 0, start_mix < 1, sample_size >= 2)
  N <- as.integer(round(breeders)); L <- as.integer(n_markers); K <- length(shared_freq)
  n2 <- max(1L, min(N - 1L, as.integer(round(start_mix * N))))
  rec <- recombination(map, L, K)
  if (!is.null(freq)) stopifnot(length(freq$freq1) == L, length(freq$freq2) == L)
  h1 <- pure_haplotypes(n2, 1L, L, shared_freq, freq, founder_impurity)
  h2 <- pure_haplotypes(N - n2, 0L, L, shared_freq, freq, founder_impurity)
  st <- list(H1 = rbind(h1$H1, h2$H1), H2 = rbind(h1$H2, h2$H2),
             anc = c(rep(1, n2), rep(0, N - n2)), age = numeric(N), t1 = numeric(N), t2 = numeric(N),
             bc = numeric(N), cls = integer(N), fam = -seq_len(N),
             pool = sample.int(max(1L, pools), N, TRUE))
  if (founders != "pure") st <- hybrid_founders(st, founders, L, rec, freq, founder_impurity)
  # summaries and true history of one sample (tad) at generation g
  record <- function(tad, g) {
    tad$G <- tad$H1 + tad$H2
    Gobs <- observe_genotypes(tad$G, L, missing, het_dropout)
    obs <- if (!is.null(microsat)) {
      a <- observe_microsat(tad$H1[, seq_len(L), drop = FALSE], tad$H2[, seq_len(L), drop = FALSE], microsat, missing)
      m <- hybrid_scores_microsat(a, microsat)
      list(S = m$S, H = m$H, n = m$n_markers)
    } else if (is.null(freq)) {
      list(S = rowMeans(Gobs[, seq_len(L), drop = FALSE], na.rm = TRUE) / 2,
           H = rowMeans(Gobs[, seq_len(L), drop = FALSE] == 1L, na.rm = TRUE),
           n = rowSums(!is.na(Gobs[, seq_len(L), drop = FALSE])))
    } else {
      m <- ml_scores(Gobs[, seq_len(L), drop = FALSE], freq$freq1, freq$freq2)
      list(S = m$S, H = m$H, n = m$n_markers)
    }
    hy <- !(tad$anc %in% c(0, 1))
    rare <- if (any(hy) && mean(tad$anc[hy]) > 0.5) 1L else 2L
    t_rare <- if (rare == 2L) tad$t2 else tad$t1
    feats <- sample_summaries(obs$S, obs$H, diag = if (is.null(microsat)) Gobs[, seq_len(L), drop = FALSE],
                              shared = if (K) Gobs[, L + seq_len(K), drop = FALSE] else NULL,
                              n_markers = obs$n, positions = if (junctions && !is.null(map)) map,
                              class_settings = class_settings)
    if (g %in% sample_at) samples[[as.character(g)]] <<- list(genotypes = tad$G, S = obs$S, H = obs$H,
                                                              class = class_names[tad$cls + 1L])
    c(generation = g, n_hybrids_true = sum(hy),
      hybrid_age_both = if (any(hy)) mean(tad$age[hy]) else NA,
      hybrid_age_rare = if (any(hy)) mean(t_rare[hy]) else NA,
      backcrosses = if (any(hy)) mean(tad$bc[hy]) else NA,
      recent_share = if (any(hy)) mean(tad$cls[hy] %in% 1:3) else NA,
      rare_species = rare, arrival_rate = recent_rate(arr[, rare], g),
      arriving = as.numeric(recent_rate(arr[, rare], g) > 0), feats)
  }
  rows <- vector("list", generations); samples <- list()
  # a sample of the founders themselves (generation 0): F1 or backcross families
  first <- if (sample_founders && founders != "pure") {
    pick <- sample.int(N, sample_size, replace = sample_size > N)
    list(record(lapply(st[c("H1", "H2", "anc", "age", "t1", "t2", "bc", "cls")], function(x)
      if (is.matrix(x)) x[pick, , drop = FALSE] else x[pick]), 0))
  }
  for (g in seq_len(generations)) {
    tad <- offspring(st, mates(sample_size, st, sib_mating, pool_crossing, assort), rec)
    rows[[g]] <- record(tad, g)
    st <- offspring(st, mates(N, st, sib_mating, pool_crossing, assort), rec)
    st <- immigrate(st, arr[g, ], L, shared_freq, pools, freq)
  }
  out <- as.data.frame(do.call(rbind, c(first, rows)))
  if (length(samples)) attr(out, "samples") <- samples
  out
}

class_names <- c("pure", "F1", "F2", "BC1", "later")

# replace the founders by F1s or first backcrosses, with matching histories
hybrid_founders <- function(st, founders, L, rec, freq = NULL, impurity = 0) {
  N <- length(st$anc)
  diag <- seq_len(L)
  hap <- function(code) pure_haplotypes(N, code, L, NULL, freq, impurity)$H1
  if (founders == "F1") {
    st$H1[, diag] <- hap(0L); st$H2[, diag] <- hap(1L)
  } else {
    # one haplotype from the pure parent, one recombinant gamete from an F1
    f1 <- list(H1 = matrix(0L, N, ncol(st$H1)), H2 = matrix(0L, N, ncol(st$H1)))
    f1$H1[, diag] <- hap(0L); f1$H2[, diag] <- hap(1L)
    g <- gametes(f1, seq_len(N), rec)
    st$H1[, diag] <- hap(if (founders == "BC1_1") 0L else 1L)
    st$H2[, diag] <- g[, diag]
  }
  # a BC1 has one F1 parent (2 generations back to each species' pure line)
  # and one pure parent (1 generation back)
  h <- switch(founders, F1 = list(anc = 0.5, t1 = 1, t2 = 1, cls = 1L, age = 1, bc = 0),
              BC1_1 = list(anc = 0.25, t1 = 4 / 3, t2 = 2, cls = 3L, age = 1.5, bc = 1),
              BC1_2 = list(anc = 0.75, t1 = 2, t2 = 4 / 3, cls = 3L, age = 1.5, bc = 1))
  st$anc[] <- h$anc; st$t1[] <- h$t1; st$t2[] <- h$t2; st$cls[] <- h$cls; st$age[] <- h$age; st$bc[] <- h$bc
  st
}

# arrivals at the end of generation j join the breeders of generation j + 1
recent_rate <- function(rate, g, window = 3) {
  if (g < 2) return(0)
  mean(rate[max(1, g - window):(g - 1)])
}

# haplotypes of pure animals: allele `code` (0 or 1) at every diagnostic
# marker, shared loci drawn from their frequencies; with `impurity`, that
# share of diagnostic alleles comes from the other species
pure_haplotypes <- function(n, code, L, shared_freq, freq = NULL, impurity = 0) {
  K <- length(shared_freq)
  diag <- function() {
    if (!is.null(freq)) {
      own <- if (code == 1L) freq$freq2 else freq$freq1; other <- if (code == 1L) freq$freq1 else freq$freq2
      return(matrix(stats::rbinom(n * L, 1, rep((1 - impurity) * own + impurity * other, each = n)), n, L))
    }
    m <- matrix(code, n, L)
    if (impurity > 0) m[stats::runif(n * L) < impurity] <- 1L - code
    m
  }
  one <- function() cbind(diag(),
                          if (K) matrix(stats::rbinom(n * K, 1, rep(shared_freq, each = n)), n, K))
  list(H1 = one(), H2 = one())
}

# recombination fraction between each column and the next: Haldane within a
# chromosome, 0.5 between chromosomes, for unmapped markers and shared loci
recombination <- function(map, L, K) {
  r <- rep(0.5, L + K - 1)
  if (!is.null(map) && L > 1) {
    stopifnot(nrow(map) == L, all(c("chromosome", "position") %in% names(map)))
    same <- map$chromosome[-1] == map$chromosome[-L]
    d <- pmax(diff(map$position), 0) / 100
    r[seq_len(L - 1)] <- ifelse(same, (1 - exp(-2 * d)) / 2, 0.5)
  }
  if (all(r == 0.5)) NULL else r
}

# first parent at random; second from the same pool unless the mating crosses
# pools, or a sibling (same mother) with probability sib_mating
mates <- function(n, st, sib_mating, pool_crossing, assort = 0) {
  N <- length(st$anc)
  p1 <- sample.int(N, n, TRUE)
  p2 <- (p1 + sample.int(N - 1L, n, TRUE) - 1L) %% N + 1L
  if (length(unique(st$pool)) > 1) {
    by_pool <- split(seq_len(N), st$pool)
    for (k in which(stats::runif(n) > pool_crossing)) {
      m <- by_pool[[as.character(st$pool[p1[k]])]]; m <- m[m != p1[k]]
      if (length(m)) p2[k] <- m[sample.int(length(m), 1L)]
    }
  }
  if (assort > 0) {
    multi <- length(unique(st$pool)) > 1
    by_pool <- if (multi) split(seq_len(N), st$pool)
    for (k in which(stats::runif(n) < assort)) {
      m <- if (multi) by_pool[[as.character(st$pool[p1[k]])]] else seq_len(N)
      m <- m[m != p1[k]]
      if (!length(m)) next
      cand <- m[sample.int(length(m), min(5L, length(m)))]
      p2[k] <- cand[which.min(abs(st$anc[cand] - st$anc[p1[k]]))]
    }
  }
  if (sib_mating > 0) {
    fams <- split(seq_len(N), st$fam)
    for (k in which(stats::runif(n) < sib_mating)) {
      m <- fams[[as.character(st$fam[p1[k]])]]; m <- m[m != p1[k]]
      if (length(m)) p2[k] <- m[sample.int(length(m), 1L)]
    }
  }
  cbind(p1, p2)
}

# one gamete from each listed parent: start on a random haplotype and switch
# between neighbouring columns with their recombination fraction
gametes <- function(st, parents, rec) {
  A <- st$H1[parents, , drop = FALSE]; B <- st$H2[parents, , drop = FALSE]
  n <- nrow(A); M <- ncol(A)
  if (is.null(rec)) {
    pick <- matrix(stats::runif(n * M) < 0.5, n, M)
  } else {
    pick <- matrix(FALSE, n, M)
    cur <- stats::runif(n) < 0.5
    pick[, 1] <- cur
    for (j in seq_len(M - 1)) {
      cur <- xor(cur, stats::runif(n) < rec[j])
      pick[, j + 1] <- cur
    }
  }
  out <- A; out[pick] <- B[pick]
  out
}

offspring <- function(st, par, rec) {
  m <- par[, 1]; f <- par[, 2]
  a <- st$anc[m]; b <- st$anc[f]
  H1 <- gametes(st, m, rec); H2 <- gametes(st, f, rec)
  pa <- a %in% c(0, 1); pb <- b %in% c(0, 1)
  same_pure <- pa & pb & a == b
  age <- ifelse(same_pure, 0, (st$age[m] + st$age[f]) / 2 + 1)
  w2 <- a + b; w1 <- 2 - a - b
  t2 <- ifelse(w2 > 0, (a * (st$t2[m] + 1) + b * (st$t2[f] + 1)) / pmax(w2, 1e-12), 0)
  t1 <- ifelse(w1 > 0, ((1 - a) * (st$t1[m] + 1) + (1 - b) * (st$t1[f] + 1)) / pmax(w1, 1e-12), 0)
  bc <- ifelse(!pa & pb, st$bc[m] + 1, ifelse(pa & !pb, st$bc[f] + 1,
         ifelse(!pa & !pb, (st$bc[m] + st$bc[f]) / 2, 0)))
  cm <- st$cls[m]; cf <- st$cls[f]
  cls <- ifelse(pa & pb, ifelse(a == b, 0L, 1L),
         ifelse(cm == 1L & cf == 1L, 2L,
         ifelse((cm == 1L & pb) | (cf == 1L & pa), 3L, 4L)))
  anc <- (a + b) / 2
  pure <- anc %in% c(0, 1)
  t1[pure] <- 0; t2[pure] <- 0; bc[pure] <- 0
  list(H1 = H1, H2 = H2, anc = anc, age = age, t1 = t1, t2 = t2, bc = bc, cls = as.integer(cls),
       fam = m, pool = st$pool[m])
}

# arriving pure animals replace residents and start with pure labels
immigrate <- function(st, arrivals, L, shared_freq, pools, freq = NULL) {
  N <- length(st$anc)
  n1 <- stats::rbinom(1, N, arrivals[1]); n2 <- stats::rbinom(1, N, arrivals[2])
  if (n1 + n2 == 0) return(st)
  idx <- sample.int(N, min(N, n1 + n2))
  sp <- c(rep(0L, min(n1, length(idx))), rep(2L, length(idx) - min(n1, length(idx))))
  for (code in c(0L, 2L)) {
    ii <- idx[sp == code]
    if (!length(ii)) next
    h <- pure_haplotypes(length(ii), code / 2L, L, shared_freq, freq)
    st$H1[ii, ] <- h$H1; st$H2[ii, ] <- h$H2
    st$anc[ii] <- code / 2; st$age[ii] <- 0; st$t1[ii] <- 0; st$t2[ii] <- 0
    st$bc[ii] <- 0; st$cls[ii] <- 0L
    st$fam[ii] <- -(max(abs(st$fam)) + seq_along(ii))
    st$pool[ii] <- sample.int(max(1L, pools), length(ii), TRUE)
  }
  st
}

# what a genotyping run would report: heterozygotes called homozygous with
# probability het_dropout, and per-animal missing data at diagnostic markers
observe_genotypes <- function(G, L, missing, het_dropout = 0) {
  if (het_dropout > 0) {
    het <- which(G == 1L)
    drop <- het[stats::runif(length(het)) < het_dropout]
    G[drop] <- ifelse(stats::runif(length(drop)) < 0.5, 0L, 2L)
  }
  rates <- if (length(missing) == 1) rep(missing, nrow(G)) else sample(missing, nrow(G), TRUE)
  for (i in which(rates > 0)) {
    n_keep <- max(10L, round(L * (1 - rates[i])))
    if (n_keep < L) G[i, sample.int(L, L - n_keep)] <- NA
  }
  G
}

with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  expr
}

#' Random marker map
#'
#' Places markers at random positions on chromosomes, for simulating linkage
#' when no genetic map is available. Use the species' chromosome number and,
#' if known, its map length; otherwise give [build_reference()] a range so
#' that the uncertainty is part of the reference.
#'
#' @param n_markers Number of markers.
#' @param chromosomes Number of chromosomes.
#' @param length_cM Total map length in centimorgans, shared equally among
#'   chromosomes.
#' @return A data.frame with `chromosome` and `position` (cM), sorted.
#' @examples
#' head(random_map(500, chromosomes = 11, length_cM = 1500))
#' @export
random_map <- function(n_markers, chromosomes, length_cM) {
  chrom <- sort(sample.int(chromosomes, n_markers, TRUE))
  pos <- stats::runif(n_markers, 0, length_cM / chromosomes)
  m <- data.frame(chromosome = chrom, position = pos)
  m[order(m$chromosome, m$position), , drop = FALSE]
}
