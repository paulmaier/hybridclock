#' Train the hybrid clock on a reference set
#'
#' Fits one Cubist model (Quinlan 1993; Kuhn and Quinlan 2023) per target,
#' from the summaries of each simulated sample to its true history, and
#' gives each estimate an interval by split-conformal calibration (Vovk et
#' al. 2005; Lei et al. 2018). The simulated runs are divided into three sets
#' so that no run informs more than one step: one to fit the models, one to
#' calibrate the intervals, and one to test the finished estimator.
#'
#' New hybridization is estimated in two parts. `arriving` is the
#' probability that pure animals of the rarer species are still arriving
#' (over the three generations before sampling), from a probability forest
#' (Wright and Ziegler 2017) trained on equal numbers of histories with and
#' without arrivals, so the prior odds are even. Its probabilities are then
#' recalibrated on the calibration runs by logistic regression on the
#' forest's log-odds, the estimated hybrid age and their interaction (Platt
#' 1999), because the same summaries mean different things in young and old
#' swarms: in a young swarm, pure founders still form new hybrids without
#' any arrivals. `arrival_rate` is the rate
#' if they are, from a Cubist model trained only on histories with arrivals.
#'
#' Ages and rates are modeled on a log scale (`log(age)`,
#' `log(1 + backcrosses)`, `log(arrival_rate)`). Fitting and calibration sets are balanced so
#' that every range of the target is equally common; otherwise the most
#' common ages in the reference set pull estimates toward them. Intervals
#' are calibrated separately within ranges of the estimate (Mondrian
#' conformal prediction), so that young and old estimates each get
#' intervals of the right width. Coverage is guaranteed on average over
#' histories like those in the reference set, which is why its settings
#' should match the study system.
#'
#' @param reference Output of [build_reference()].
#' @param targets Which parts of the history to estimate. See
#'   [simulate_hybrid_zone()] for definitions. `arriving` and
#'   `arrival_rate` describe new hybridization: whether pure animals of the
#'   rarer species are still arriving, and if so what share of breeders
#'   they make up per generation.
#' @param level Interval coverage.
#' @param split Shares of runs used to fit, calibrate and test.
#' @param max_rows Most samples used to fit each model.
#' @param committees Number of Cubist committees.
#' @param bins Number of ranges of the estimate within which intervals are
#'   calibrated.
#' @param intervals How intervals are made: `"binned"` (split-conformal
#'   residuals within ranges of the estimate and by sample type) or `"cqr"`
#'   (conformalized quantile regression: a quantile forest gives each sample
#'   its own interval from its summaries, widened on calibration runs to the
#'   stated coverage; Romano et al. 2019), or `"normalized"` (intervals
#'   centred on the estimate, scaled by a forest's prediction of each
#'   sample's typical error, and calibrated on held-out runs; Papadopoulos
#'   et al. 2011).
#' @param cores Number of processes (targets are fitted in parallel).
#' @param seed Random seed.
#' @return An object of class `hybridclock_model`: the fitted models, the
#'   interval calibration, the features used, the reference settings and a
#'   `validation` table of accuracy on the test runs.
#' @references
#'   Kuhn M, Quinlan R (2023) Cubist: Rule- and instance-based regression
#'   modeling. R package.
#'
#'   Lei J, G'Sell M, Rinaldo A, Tibshirani RJ, Wasserman L (2018)
#'   Distribution-free predictive inference for regression. Journal of the
#'   American Statistical Association 113:1094-1111.
#'
#'   Quinlan JR (1993) Combining instance-based and model-based learning.
#'   Proceedings of the Tenth International Conference on Machine Learning,
#'   236-243.
#'
#'   Vovk V, Gammerman A, Shafer G (2005) Algorithmic Learning in a Random
#'   World. Springer.
#'
#'   Platt J (1999) Probabilistic outputs for support vector machines and
#'   comparisons to regularized likelihood methods. In: Advances in Large
#'   Margin Classifiers, MIT Press, 61-74.
#'
#'   Wright MN, Ziegler A (2017) ranger: A fast implementation of random
#'   forests for high dimensional data in C++ and R. Journal of Statistical
#'   Software 77:1-17.
#' @examples
#' \donttest{
#' ref <- build_reference(n_sims = 150, generations = 12, n_markers = 200, breeders = c(20, 100))
#' clock <- train_clock(ref, targets = "hybrid_age_both", committees = 5)
#' clock
#' }
#' @export
train_clock <- function(reference, targets = c("hybrid_age_both", "hybrid_age_rare", "backcrosses", "arriving",
                                               "arrival_rate"),
                        level = 0.9, split = c(fit = 0.6, calibrate = 0.2, test = 0.2), max_rows = 40000,
                        committees = 20, bins = 8, intervals = c("normalized", "binned", "cqr"), cores = 1, seed = 1) {
  targets <- match.arg(targets, all_targets, several.ok = TRUE)
  intervals <- match.arg(intervals)
  with_seed(seed, {
    ref <- reference[reference$n_hybrids > 0, , drop = FALSE]
    feats <- intersect(feature_names(), names(ref))
    feats <- feats[vapply(ref[feats], function(x) any(is.finite(x)), NA)]
    runs <- unique(ref$run)
    counts <- floor(split / sum(split) * length(runs))
    counts[1] <- counts[1] + length(runs) - sum(counts)
    part <- sample(rep(names(split), counts))
    sets <- lapply(stats::setNames(names(split), names(split)), function(s) ref$run %in% runs[part == s])
    fit_one <- function(target) with_seed(seed + match(target, targets), {
      z <- to_scale(target, ref[[target]]); ok <- is.finite(z)
      fit <- balance_rows(z, which(ok & sets$fit), max_rows)
      cal <- balance_rows(z, which(ok & sets$calibrate), max_rows)
      tst <- balance_rows(z, which(ok & sets$test), max_rows)
      m <- Cubist::cubist(prep_features(ref[fit, ], feats), z[fit], committees = committees)
      Xcal <- prep_features(ref[cal, ], feats)
      pc <- stats::predict(m, Xcal)
      calib <- switch(intervals,
        cqr = cqr_fit(prep_features(ref[fit, ], feats), z[fit], Xcal, z[cal], level, sample_type(ref[cal, ])),
        normalized = normalized_fit(Xcal, z[cal] - pc, level, sample_type(ref[cal, ])),
        binned = conformal_fit(pc, z[cal] - pc, level, bins, sample_type(ref[cal, ])))
      val <- NULL
      if (length(tst)) {
        Xtst <- prep_features(ref[tst, ], feats)
        pt <- stats::predict(m, Xtst)
        iv <- apply_conformal(pt, calib, sample_type(ref[tst, ]), Xtst)
        val <- validate(target, ref[[target]][tst], from_scale(target, pt),
                        from_scale(target, iv$lower), from_scale(target, iv$upper))
      }
      list(model = m, calibration = calib, validation = val)
    })
    constant <- vapply(targets, function(t) {
      z <- if (t == "arriving") ref$arriving else to_scale(t, ref[[t]])
      length(unique(z[is.finite(z)])) < 2 }, NA)
    if (any(constant)) {
      warning("Skipping targets that do not vary in the reference: ",
              paste(names(constant)[constant], collapse = ", "), call. = FALSE)
      targets <- setdiff(targets, names(constant)[constant])
    }
    main <- setdiff(targets, "arriving")
    fits <- if (cores > 1 && .Platform$OS.type != "windows") {
      parallel::mclapply(main, fit_one, mc.cores = min(cores, length(main)))
    } else lapply(main, fit_one)
    names(fits) <- main
    if ("arriving" %in% targets) {
      age_model <- if ("hybrid_age_both" %in% main) fits$hybrid_age_both$model else
        Cubist::cubist(prep_features(ref[sets$fit, ], feats), to_scale("hybrid_age_both", ref$hybrid_age_both[sets$fit]),
                       committees = 5)
      fits$arriving <- with_seed(seed, fit_arriving(ref, feats, sets, max_rows, age_model))
      fits <- fits[targets]
    }
    clock <- structure(list(models = lapply(fits, `[[`, "model"), calibration = lapply(fits, `[[`, "calibration"),
                            validation = do.call(rbind, lapply(fits, `[[`, "validation")),
                            features = feats, targets = targets, level = level,
                            settings = attr(reference, "settings")),
                       class = "hybridclock_model")
    # when to trust individual classes instead of the clock (calibration runs)
    if ("hybrid_age_both" %in% targets && all(paste0("class_", early_classes$class) %in% feats))
      clock$tier <- choose_tier(clock, ref[sets$calibrate, , drop = FALSE], level = level)
    # spread of hybrid ages the reference allows, to judge how much an
    # answer owes to the data rather than to the reference's range
    if ("hybrid_age_both" %in% targets)
      clock$age_range <- stats::quantile(log(pmax(ref$hybrid_age_both, 1)), c(0.05, 0.95), names = FALSE,
                                         na.rm = TRUE)
    clock
  })
}

#' @export
print.hybridclock_model <- function(x, ...) {
  s <- x$settings
  cat("hybridclock model:", paste(x$targets, collapse = ", "), "\n")
  if (!is.null(s)) cat("  reference:", s$n_markers, "markers,", s$n_shared, "shared loci, samples of",
                       s$sample_size, "animals, breeders", s$breeders[1], "to", s$breeders[2], "\n")
  v <- x$validation
  if (is.null(v)) return(invisible(x))
  cat("  ", round(100 * x$level), "% intervals; accuracy on held-out runs:\n", sep = "")
  overall <- do.call(rbind, lapply(split(v, v$target), function(d) data.frame(
    target = d$target[1], typical_error = stats::weighted.mean(d$typical_error, d$n),
    coverage = stats::weighted.mean(d$coverage, d$n))))
  print(overall, row.names = FALSE, digits = 2)
  if ("arriving" %in% x$targets) cat("  (for arriving, coverage is the share called correctly at p > 0.5)\n")
  if (!is.null(x$tier)) {
    if (is.finite(x$tier$cutoff)) cat("  individual classes are used when at least ", round(100 * x$tier$cutoff),
                                      "% of hybrids fit an early class\n", sep = "")
    else cat("  individual classes are not used by share of early hybrids: with these markers they could not date",
             "simulated samples\n")
    if (isTRUE(x$tier$f1_ok)) cat("  samples made up mostly (80% or more) of F1s are dated by their classes\n")
  }
  invisible(x)
}

#' Predict histories from sample summaries
#'
#' @param object Output of [train_clock()].
#' @param newdata A data.frame (or matrix) of summaries, one row per sample,
#'   as from [sample_summaries()] or [build_reference()].
#' @param ... Not used.
#' @return A data.frame with an estimate and interval (`_lower`, `_upper`)
#'   for each target.
#' @examples
#' \donttest{
#' ref <- build_reference(n_sims = 150, generations = 12, n_markers = 200, breeders = c(20, 100))
#' clock <- train_clock(ref, targets = "hybrid_age_both", committees = 5)
#' head(predict(clock, ref))
#' }
#' @export
predict.hybridclock_model <- function(object, newdata, ...) {
  X <- prep_features(newdata, object$features)
  out <- lapply(object$targets, function(target) {
    if (target == "arriving") return(data.frame(arriving = arriving_prob(object$models$arriving, X)))
    z <- stats::predict(object$models[[target]], X)
    iv <- apply_conformal(z, object$calibration[[target]], sample_type(newdata), X)
    stats::setNames(data.frame(from_scale(target, z), from_scale(target, iv$lower), from_scale(target, iv$upper)),
                    c(target, paste0(target, c("_lower", "_upper"))))
  })
  do.call(cbind, out)
}

arriving_inputs <- function(m, X) {
  p <- stats::predict(m$forest, X, num.threads = 1)$predictions[, "1"]
  data.frame(logit = stats::qlogis(pmin(pmax(p, 0.005), 0.995)),
             log_age = stats::predict(m$age_model, X))
}

arriving_prob <- function(m, X) {
  d <- arriving_inputs(m, X)
  b <- m$platt
  stats::plogis(b[1] + b[2] * d$logit + b[3] * d$log_age + b[4] * d$logit * d$log_age)
}

feature_names <- function() names(sample_summaries(c(0, 0.5), c(0, 1)))

prep_features <- function(x, feats) {
  x <- as.data.frame(x)[, feats, drop = FALSE]
  x[] <- lapply(x, function(v) { v[!is.finite(v)] <- -1; v })
  x
}

all_targets <- c("hybrid_age_both", "hybrid_age_rare", "backcrosses", "arriving", "arrival_rate")

to_scale <- function(target, y) switch(target,
  hybrid_age_both = , hybrid_age_rare = log(pmax(y, 1)), backcrosses = log1p(y),
  arrival_rate = ifelse(y > 0, log(y), NA))
from_scale <- function(target, z) switch(target,
  hybrid_age_both = , hybrid_age_rare = pmax(exp(z), 1), backcrosses = pmax(expm1(z), 0), arrival_rate = exp(z))

# probability that the rarer species is still arriving: a probability forest
# on equal numbers of histories with and without arrivals
fit_arriving <- function(ref, feats, sets, max_rows, age_model) {
  y <- ref$arriving
  even <- function(rows) {
    by <- split(rows, y[rows]); n <- min(lengths(by), floor(max_rows / 2))
    unlist(lapply(by, function(i) i[sample.int(length(i), n)]), use.names = FALSE)
  }
  fit <- even(which(sets$fit)); cal <- even(which(sets$calibrate)); tst <- even(which(sets$test))
  forest <- ranger::ranger(x = prep_features(ref[fit, ], feats), y = factor(y[fit], levels = c(0, 1)),
                           num.trees = 500, probability = TRUE, num.threads = 1, verbose = FALSE)
  m <- list(forest = forest, age_model = age_model, platt = NULL)
  d <- arriving_inputs(m, prep_features(ref[cal, ], feats))
  d$y <- y[cal]
  m$platt <- stats::coef(stats::glm(y ~ logit * log_age, family = stats::binomial(), data = d))
  val <- NULL
  if (length(tst)) {
    p <- arriving_prob(m, prep_features(ref[tst, ], feats))
    truth <- y[tst]
    val <- do.call(rbind, lapply(c(0, 1), function(k) data.frame(
      target = "arriving", true_range = c("not arriving", "arriving")[k + 1], n = sum(truth == k),
      median_estimate = stats::median(p[truth == k]), typical_error = NA,
      coverage = mean((p[truth == k] > 0.5) == (k == 1)), median_width = NA)))
  }
  list(model = m, calibration = NULL, validation = val)
}

balance_rows <- function(z, rows, max_rows, n_bins = 20) {
  if (length(rows) <= max_rows) return(rows)
  b <- cut(z[rows], seq(min(z[rows]), max(z[rows]), length.out = n_bins + 1), include.lowest = TRUE)
  per <- ceiling(max_rows / nlevels(b))
  unlist(lapply(split(rows, b), function(i) if (length(i) > per) i[sample.int(length(i), per)] else i),
         use.names = FALSE)
}

conformal_bins <- function(pred, resid, level, bins) {
  edges <- unique(stats::quantile(pred, seq(0, 1, length.out = bins + 1), names = FALSE))
  edges[1] <- -Inf; edges[length(edges)] <- Inf
  b <- cut(pred, edges, include.lowest = TRUE)
  a <- (1 - level) / 2
  q <- t(vapply(split(resid, b), function(r) {
    n <- length(r)
    if (n < 10) return(c(NA_real_, NA_real_))
    c(stats::quantile(r, max(0, floor((n + 1) * a) / n), names = FALSE, type = 1),
      stats::quantile(r, min(1, ceiling((n + 1) * (1 - a)) / n), names = FALSE, type = 1))
  }, numeric(2)))
  all_q <- stats::quantile(resid, c(a, 1 - a), names = FALSE)
  q[is.na(q[, 1]), ] <- rep(all_q, each = sum(is.na(q[, 1])))
  list(edges = edges, lower = q[, 1], upper = q[, 2])
}

apply_conformal <- function(pred, calib, type = NULL, X = NULL) {
  if (identical(calib$method, "cqr")) return(cqr_apply(pred, calib, type, X))
  if (identical(calib$method, "normalized")) return(normalized_apply(pred, calib, type, X))
  one <- function(cb, p) {
    b <- as.integer(cut(p, cb$edges, include.lowest = TRUE))
    list(lower = p + cb$lower[b], upper = p + cb$upper[b])
  }
  if (is.null(calib$base)) return(one(calib, pred))   # clocks trained before sample types
  out <- one(calib$base, pred)
  if (!is.null(type)) for (t in intersect(unique(type), names(calib$strata))) {
    i <- which(type == t); r <- one(calib$strata[[t]], pred[i])
    out$lower[i] <- r$lower; out$upper[i] <- r$upper
  }
  out
}

# intervals calibrated separately for sample types whose errors differ: samples
# made up mostly of early hybrids that are mostly F2-like (F2 families and
# later generations of closed swarms look alike) versus all others; a type
# needs enough calibration samples, otherwise the shared calibration is used
conformal_fit <- function(pred, resid, level, bins, type, min_n = 200) {
  out <- list(base = conformal_bins(pred, resid, level, bins), strata = list())
  for (t in unique(type)) {
    i <- which(type == t)
    if (length(i) >= min_n) out$strata[[t]] <- conformal_bins(pred[i], resid[i], level, bins)
  }
  out
}

sample_type <- function(x) {
  x <- as.data.frame(x)
  if (!all(paste0("class_", c(early_classes$class, "later")) %in% names(x))) return(rep("other", nrow(x)))
  ca <- class_ages(x)
  ifelse(ca$early_share >= min_early_share & ca$f2_share > max_f2_share & !is.na(ca$f2_share), "f2_like", "other")
}

validate <- function(target, y, est, lower, upper) {
  br <- switch(target, hybrid_age_both = , hybrid_age_rare = c(1, 2.5, 4, 6, 9, 14, 20, 30, Inf),
               backcrosses = c(-Inf, 0.5, 1.5, 2.5, 4, 6, 9, Inf),
               arrival_rate = c(0, 0.005, 0.01, 0.02, 0.04, 0.06, Inf))
  b <- cut(y, br, include.lowest = TRUE)
  err <- switch(target, backcrosses = abs(log((est + 1) / (y + 1))), abs(log(est / y)))
  do.call(rbind, lapply(split(seq_along(y), b), function(i) if (length(i)) data.frame(
    target = target, true_range = as.character(b[i[1]]), n = length(i),
    median_estimate = stats::median(est[i]), typical_error = stats::median(err[i]),
    coverage = mean(y[i] >= lower[i] & y[i] <= upper[i]),
    median_width = stats::median(upper[i] - lower[i]))))
}

# conformalized quantile regression (Romano, Patterson & Candes 2019): a
# quantile forest gives each sample its own interval from its summaries, and
# the calibration runs widen all intervals by the amount needed for the
# stated coverage (separately by sample type when there are enough runs)
cqr_fit <- function(Xfit, zfit, Xcal, zcal, level, type, min_n = 200, max_fit = 20000) {
  a <- (1 - level) / 2
  if (nrow(Xfit) > max_fit) { i <- sample.int(nrow(Xfit), max_fit); Xfit <- Xfit[i, , drop = FALSE]; zfit <- zfit[i] }
  qf <- ranger::ranger(x = Xfit, y = zfit, quantreg = TRUE, num.trees = 200, min.node.size = 10, num.threads = 1,
                       verbose = FALSE)
  q <- stats::predict(qf, Xcal, type = "quantiles", quantiles = c(a, 1 - a), num.threads = 1)$predictions
  score <- pmax(q[, 1] - zcal, zcal - q[, 2])
  corr <- function(e) { n <- length(e); stats::quantile(e, min(1, ceiling((n + 1) * level) / n), names = FALSE, type = 1) }
  out <- list(method = "cqr", forest = qf, a = a, base = corr(score), strata = list())
  for (t in unique(type)) { i <- which(type == t); if (length(i) >= min_n) out$strata[[t]] <- corr(score[i]) }
  out
}

cqr_apply <- function(pred, calib, type = NULL, X = NULL) {
  if (is.null(X)) stop("Quantile-forest intervals need the summaries.", call. = FALSE)
  q <- stats::predict(calib$forest, X, type = "quantiles", quantiles = c(calib$a, 1 - calib$a),
                      num.threads = 1)$predictions
  e <- rep(calib$base, length(pred))
  if (!is.null(type)) for (t in intersect(unique(type), names(calib$strata))) e[type == t] <- calib$strata[[t]]
  list(lower = pmin(q[, 1] - e, pred), upper = pmax(q[, 2] + e, pred))
}

# normalized split-conformal intervals: half of the calibration runs teach a
# forest how large each sample's error tends to be (from its summaries); the
# other half sets how many such units the interval needs on each side
normalized_fit <- function(Xcal, resid, level, type, min_n = 200) {
  a <- (1 - level) / 2
  h <- sample(rep(1:2, length.out = length(resid)))
  sf <- ranger::ranger(x = Xcal[h == 1, , drop = FALSE], y = abs(resid[h == 1]), num.trees = 200,
                       min.node.size = 10, num.threads = 1, verbose = FALSE)
  sig <- pmax(stats::predict(sf, Xcal[h == 2, , drop = FALSE], num.threads = 1)$predictions, 1e-3)
  sc <- resid[h == 2] / sig; ty <- type[h == 2]
  qq <- function(e) { n <- length(e)
    c(stats::quantile(e, max(0, floor((n + 1) * a) / n), names = FALSE, type = 1),
      stats::quantile(e, min(1, ceiling((n + 1) * (1 - a)) / n), names = FALSE, type = 1)) }
  out <- list(method = "normalized", forest = sf, base = qq(sc), strata = list())
  for (t in unique(ty)) { i <- which(ty == t); if (length(i) >= min_n / 2) out$strata[[t]] <- qq(sc[i]) }
  out
}

normalized_apply <- function(pred, calib, type = NULL, X = NULL) {
  if (is.null(X)) stop("Normalized intervals need the summaries.", call. = FALSE)
  sig <- pmax(stats::predict(calib$forest, X, num.threads = 1)$predictions, 1e-3)
  lo <- rep(calib$base[1], length(pred)); hi <- rep(calib$base[2], length(pred))
  if (!is.null(type)) for (t in intersect(unique(type), names(calib$strata))) {
    i <- type == t; lo[i] <- calib$strata[[t]][1]; hi[i] <- calib$strata[[t]][2] }
  list(lower = pmin(pred + lo * sig, pred), upper = pmax(pred + hi * sig, pred))
}
