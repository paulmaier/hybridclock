#' Estimate the hybrid history of sampled populations
#'
#' Gives one answer per population, from one of two tiers:
#' * **Individual classes.** When most hybrids confidently fit an early class
#'   (F1, F2, first or second backcross; see [classify_hybrids()]) and those
#'   are mostly F1s and backcrosses rather than F2-like animals, the ages
#'   follow directly from the classes: F1 = 1, first backcross = 1.5, second
#'   backcross = 1.75, F2 = 2 generations (hybrid age, both species), and
#'   likewise for the rare species' lines and the number of backcrosses.
#'   Intervals come from resampling animals, widened to at least the error
#'   of class-based ages in the reference at a similar share of early hybrids.
#' * **Population clock.** Otherwise the trained clock is used: the sample's
#'   summaries are compared with the simulated reference.
#'
#' F2-like animals are not trusted as second-generation hybrids on their own:
#' in a swarm that mates among itself, F3, F4 and later generations have the
#' same expected hybrid index and heterozygosity as F2s, so one animal cannot
#' show its generation. Only the population can, through linkage blocks that
#' shrink and drift that accumulates over generations, which is what the
#' clock reads.
#'
#' The share of hybrids in early classes at which the answer switches to the
#' classes was chosen when the clock was trained, as the point above which
#' the classes' large errors (90th percentile) on held-out simulated
#' histories were no worse than the clock's, or under about 30%. With few or
#' weak markers, later-generation animals still fit early classes, so the
#' classes are never used. Column
#' `method` says which tier gave the answer; the other tier's estimates are
#' kept in columns starting `clock_` and `classes_`.
#'
#' Populations sampled more heavily than the reference are summarized in
#' random subsets of the reference sample size and the results averaged (on
#' the model's scale), so that every summary is computed exactly as in
#' training.
#'
#' @param clock Output of [train_clock()].
#' @param scores Output of [hybrid_scores()] (columns `S`, `H`, `n_markers`).
#' @param population Population label of each animal.
#' @param diag Genotypes at the diagnostic markers ([orient_genotypes()]),
#'   needed if the reference used them (it always does).
#' @param shared Genotypes at shared loci, needed if the reference had
#'   shared loci.
#' @param positions Marker positions matching the columns of `diag`
#'   (`chromosome`, `position`), needed if the reference used a fixed map.
#' @param markers,effective_markers,het_dropout No longer used, kept so that
#'   older scripts run. The classes that decide between the two tiers are
#'   made exactly as in the reference (from S, H and the number of markers,
#'   with default settings), so that the rule learned on the reference
#'   applies; use [classify_hybrids()] with these settings for the
#'   individual calls.
#' @param repeats Number of random subsets for larger populations.
#' @param seed Random seed.
#' @return A data.frame with one row per population: `n`, `n_hybrids`
#'   (animals not classified as pure), `rare_species` (the species with less
#'   ancestry in the hybrids), `rare_share` (its mean share in the hybrids),
#'   `method`, `early_share` (share of hybrids in early classes), an estimate
#'   and interval (`_lower`, `_upper`) for `hybrid_age_both`,
#'   `hybrid_age_rare` and `backcrosses`, the arrival probability
#'   `arriving` and rate `arrival_rate` (always from the clock), and notes.
#'   `arriving_note` flags young swarms (hybrid age under 5 generations),
#'   where pure founders of the rarer species still form new hybrids without
#'   any arrivals and the arrival call was right only 63-67% of the time in
#'   validation. `information` is 1 minus the width of the hybrid-age
#'   interval relative to the spread of ages in the reference (log scale):
#'   near 1 the data pin the age down, near 0 the answer mostly reflects the
#'   reference's age range (few or weak markers), and `note` says so. `note`
#'   also flags populations with too few hybrids.
#' @examples
#' \donttest{
#' ex <- hybridclock_example()
#' mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
#' sc <- hybrid_scores(ex$genotypes, mk)
#' ref <- build_reference(n_sims = 150, generations = 12, n_markers = nrow(mk),
#'                        breeders = c(20, 100))
#' clock <- train_clock(ref, targets = "hybrid_age_both", committees = 5)
#' hyb <- ex$groups %in% c("site_A", "site_B")
#' estimate_history(clock, sc[hyb, ], ex$groups[hyb], diag = orient_genotypes(ex$genotypes, mk)[hyb, ])
#' }
#' @export
estimate_history <- function(clock, scores, population, diag = NULL, shared = NULL, markers = NULL,
                             positions = NULL, effective_markers = NULL, het_dropout = 0,
                             repeats = 30, seed = 1) {
  stopifnot(inherits(clock, "hybridclock_model"), nrow(scores) == length(population))
  need_diag <- "marker_spread" %in% clock$features
  need_shared <- "shared_fis" %in% clock$features
  need_n <- any(class_features %in% clock$features)
  need_pos <- "switch_rate" %in% clock$features
  if (need_pos && (is.null(positions) || is.null(diag)))
    stop("This clock uses genotype switches along chromosomes; supply `diag` and `positions`.", call. = FALSE)
  if (need_n && is.null(scores$n_markers)) stop("This clock uses `scores$n_markers`; use hybrid_scores().", call. = FALSE)
  if (need_diag && is.null(diag)) stop("This clock uses `diag` genotypes; supply them.", call. = FALSE)
  if (need_shared && is.null(shared)) stop("This clock was trained with shared loci; supply `shared`.", call. = FALSE)
  # animals without scores (no marker genotyped) are set aside
  ok <- is.finite(scores$S) & is.finite(scores$H)
  if (!all(ok)) {
    warning(sum(!ok), " animals without hybrid scores (no markers genotyped) were left out.", call. = FALSE)
    scores <- scores[ok, , drop = FALSE]; population <- population[ok]
    if (!is.null(diag)) diag <- diag[ok, , drop = FALSE]
    if (!is.null(shared)) shared <- shared[ok, , drop = FALSE]
  }
  n_ref <- clock$settings$sample_size
  cutoff <- if (is.null(clock$tier)) Inf else clock$tier$cutoff
  age_targets <- intersect(c("hybrid_age_both", "hybrid_age_rare", "backcrosses"), clock$targets)
  with_seed(seed, {
    # classes made exactly as in the reference summaries (from S, H and the
    # number of markers, default settings), so that the rule for choosing
    # between the classes and the clock, learned on the reference, applies
    cs <- clock$settings$class_settings
    cls_all <- classify_counts(class_counts(scores$S, scores$H, scores$n_markers, diag, cs), cs)
    rows <- lapply(split(seq_along(population), population), function(idx) {
      n <- length(idx)
      S <- scores$S[idx]; H <- scores$H[idx]
      hy <- !(cls_all[idx] %in% c("pure1", "pure2"))
      rare <- if (any(hy) && mean(S[hy]) > 0.5) "species1" else "species2"
      base <- data.frame(population = population[idx[1]], n = n, n_hybrids = sum(hy), rare_species = rare,
                         rare_share = if (any(hy)) mean(if (rare == "species1") 1 - S[hy] else S[hy]) else NA,
                         method = NA_character_, early_share = NA_real_, f2_share = NA_real_, note = "")
      # class-based tier
      k <- orient_classes(cls_all[idx], rare)
      ca <- class_ages(class_shares(k))
      base$early_share <- ca$early_share; base$f2_share <- ca$f2_share
      boot <- if (sum(hy) >= 3) do.call(rbind, lapply(seq_len(500), function(b)
        class_ages(class_shares(k[sample.int(n, n, TRUE)])))) else NULL
      # intervals: resampling animals, widened to at least the class-based
      # error seen in the reference at this early share
      f1_rule <- isTRUE(clock$tier$f1_ok) && isTRUE(ca$f1_share >= min_f1_share) &&
        isTRUE(ca$early_share >= min_early_share)
      q <- if (f1_rule) max(clock$tier$f1_spread, class_floor) else class_spread(clock$tier, ca$early_share)
      for (t in age_targets) {
        est <- ca[[t]]
        lo <- if (is.null(boot)) NA else stats::quantile(boot[[t]], 0.05, na.rm = TRUE, names = FALSE)
        hi <- if (is.null(boot)) NA else stats::quantile(boot[[t]], 0.95, na.rm = TRUE, names = FALSE)
        if (!is.null(boot) && is.finite(est)) {
          z <- to_scale(t, est)
          lo <- min(lo, from_scale(t, z - q)); hi <- max(hi, from_scale(t, z + q))
        }
        base[[paste0("classes_", t)]] <- est
        base[[paste0("classes_", t, "_lower")]] <- lo
        base[[paste0("classes_", t, "_upper")]] <- hi
      }
      # clock tier
      if (n >= n_ref && any(hy)) {
        draws <- if (n == n_ref) list(seq_len(n)) else replicate(repeats, sample.int(n, n_ref), simplify = FALSE)
        X <- do.call(rbind, lapply(draws, function(j) sample_summaries(
          S[j], H[j], diag = if (need_diag) diag[idx[j], , drop = FALSE],
          shared = if (need_shared) shared[idx[j], , drop = FALSE],
          n_markers = if (need_n) scores$n_markers[idx[j]], positions = if (need_pos) positions,
          class_settings = cs)))
        X <- prep_features(X, clock$features)
        for (target in clock$targets) {
          if (target == "arriving") { base$arriving <- mean(arriving_prob(clock$models$arriving, X)); next }
          z <- stats::predict(clock$models[[target]], X)
          iv <- apply_conformal(z, clock$calibration[[target]], sample_type(X), X)
          nm <- if (target %in% age_targets) paste0("clock_", target) else target
          base[[nm]] <- from_scale(target, mean(z))
          base[[paste0(nm, "_lower")]] <- from_scale(target, mean(iv$lower))
          base[[paste0(nm, "_upper")]] <- from_scale(target, mean(iv$upper))
        }
      } else if (n < n_ref) {
        warning("Population ", population[idx[1]], " has ", n, " animals, fewer than the ", n_ref,
                " the clock was trained on; build a reference with sample_size = ", n, ".", call. = FALSE)
      }
      # choose the tier
      use_classes <- sum(hy) >= 3 && (isTRUE(ca$early_share >= cutoff) && isTRUE(ca$f2_share <= max_f2_share) || f1_rule)
      use_clock <- !use_classes && !is.null(base[[paste0("clock_", age_targets[1])]])
      base$method <- if (sum(hy) < 3) NA_character_ else if (use_classes) "individual classes" else
        if (use_clock) "population clock" else NA_character_
      if (sum(hy) < 3) base$note <- "too few hybrids for a population estimate"
      for (t in age_targets) for (sfx in c("", "_lower", "_upper")) {
        src <- if (use_classes) paste0("classes_", t, sfx) else paste0("clock_", t, sfx)
        base[[paste0(t, sfx)]] <- if (!is.na(base$method) && !is.null(base[[src]])) base[[src]] else NA
      }
      # information: 1 - interval width / spread of ages in the reference (log
      # scale); near 0, the answer mostly reflects the reference's age range
      if (!is.null(clock$age_range) && isTRUE(base$hybrid_age_both_upper > 0)) {
        w <- log(base$hybrid_age_both_upper) - log(base$hybrid_age_both_lower)
        base$information <- max(0, min(1, 1 - w / diff(clock$age_range)))
        if (base$information < low_information && !nzchar(base$note))
          base$note <- "little information: the answer mostly reflects the reference's age range; narrow `generations` if the history allows"
      }
      if ("arriving" %in% clock$targets && !is.null(base$arriving))
        base$arriving_note <- if (isTRUE(base$hybrid_age_both < 5)) "uncertain: young swarm" else ""
      base
    })
    cols <- unique(unlist(lapply(rows, names)))
    out <- do.call(rbind, lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
    front <- c("population", "n", "n_hybrids", "rare_species", "rare_share", "method", "early_share", "f2_share",
               unlist(lapply(age_targets, function(t) paste0(t, c("", "_lower", "_upper")))))
    out <- out[, c(intersect(front, names(out)), setdiff(names(out), front))]
    rownames(out) <- NULL
    out
  })
}

# relabel classes so that species 2 is the rarer species in the population
orient_classes <- function(k, rare) {
  if (rare == "species2") {
    map <- c(pure1 = "pure_common", pure2 = "pure_rare", BC1_species1 = "BC1_common", BC1_species2 = "BC1_rare",
             BC2_species1 = "BC2_common", BC2_species2 = "BC2_rare")
  } else {
    map <- c(pure1 = "pure_rare", pure2 = "pure_common", BC1_species1 = "BC1_rare", BC1_species2 = "BC1_common",
             BC2_species1 = "BC2_rare", BC2_species2 = "BC2_common")
  }
  ifelse(k %in% names(map), map[k], k)
}

class_shares <- function(k) {
  lev <- sub("^class_", "", class_features)
  stats::setNames(as.list(vapply(lev, function(l) mean(k == l), 0)), class_features)
}

# below this, an answer is flagged as mostly reflecting the reference's age range
low_information <- 0.4
