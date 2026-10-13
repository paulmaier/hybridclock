# Validation study for hybridclock with generic settings: 1,000 diagnostic
# markers, 200 shared loci, samples of 20 animals, 10% missing data.
# A quarter of the simulated histories is set aside; clocks trained on the
# rest with different assumed population-size ranges are scored on it.
# Usage (from the package folder):
#   Rscript data-raw/validation/run_validation.R <output folder> <cores>
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "data-raw/validation/output"
cores <- if (length(args) > 1) as.integer(args[2]) else 4
dir.create(out, showWarnings = FALSE, recursive = TRUE)
devtools::load_all(".", quiet = TRUE)
set.seed(1)
shared_freq <- runif(200, 0.05, 0.95)
ref <- build_reference(n_sims = 16000, generations = 50, breeders = c(10, 500), n_markers = 1000,
                       shared_freq = shared_freq, sample_size = 20, missing = 0.1, cores = cores, seed = 3e6)
settings <- attr(ref, "settings")
runs <- unique(ref$run); hold <- sample(runs, round(0.25 * length(runs)))
pool <- ref[!(ref$run %in% hold), ]; test <- ref[ref$run %in% hold & ref$n_hybrids > 0, ]
targets <- c("hybrid_age_both", "hybrid_age_rare", "backcrosses", "arrival_rate")
train_range <- function(lo, hi) {
  p <- pool[pool$breeders >= lo & pool$breeders <= hi, ]
  attr(p, "settings") <- modifyList(settings, list(breeders = c(lo, hi)))
  train_clock(p, targets, split = c(fit = 0.75, calibrate = 0.25, test = 0), cores = min(cores, 4))
}
clocks <- list(wide_10_500 = train_range(10, 500), narrow_25_100 = train_range(25, 100))
pw <- predict(clocks$wide_10_500, test)
overall <- do.call(rbind, lapply(targets, function(t) validate(
  t, test[[t]], pw[[t]], pw[[paste0(t, "_lower")]], pw[[paste0(t, "_upper")]])))
write.csv(overall, file.path(out, "overall.csv"), row.names = FALSE)
band <- cut(test$breeders, c(9, 20, 24, 100, 199, 500), labels = c("10-20", "21-24", "25-100", "101-199", "200-500"))
kind <- ifelse(test$kind == "closed", "closed swarm", "arrivals")
sens <- list()
for (target in c("hybrid_age_both", "hybrid_age_rare")) for (b in c("10-20", "25-100", "200-500"))
  for (k in c("closed swarm", "arrivals")) for (cn in names(clocks)) {
    rows <- which(band == b & kind == k & is.finite(test[[target]]) & test[[target]] >= 5)
    p <- predict(clocks[[cn]], test[rows, ]); y <- test[[target]][rows]; e <- p[[target]]
    sens[[length(sens) + 1]] <- data.frame(target = target, true_breeders = b, history = k, clock = cn,
      n = length(y), ratio = median(e / y), typical_error = median(abs(log(e / y))),
      coverage = mean(y >= p[[paste0(target, "_lower")]] & y <= p[[paste0(target, "_upper")]]))
  }
write.csv(do.call(rbind, sens), file.path(out, "size_sensitivity.csv"), row.names = FALSE)
writeLines(c(paste("R", getRversion()), paste("hybridclock", packageVersion("hybridclock")),
             paste("Cubist", packageVersion("Cubist")), format(Sys.time())), file.path(out, "RUN_INFO.txt"))
