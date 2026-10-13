# Builds data/wolves.rda: Mexican wolves from crosses between the captive
# McBride and Ghost Ranch lineages.
#
# Source: Fitak RR, Rinkevich SE, Culver M (2018) Journal of Heredity,
# doi:10.1093/jhered/esy009. Data: Dryad doi:10.5061/dryad.g68k008 (CC0).
#
# The genotypes, groups and per-animal truth come from the truth-set loader
# (load.R and TRUTH.md in the folder below), which reads the Dryad files and
# checks each wolf's pedigree class against its genotypes. Only the class and
# hybrid age of each wolf are kept from the pedigree.
#
# Run from the package directory: Rscript data-raw/wolves.R

truth_dir <- path.expand(Sys.getenv("HYBRIDCLOCK_TRUTH_DATA", "~/scripts/hybridclock_analysis/truth_data"))
source(file.path(truth_dir, "wolf_mexican_fitak2018", "load.R"))
w <- load_wolf_mexican_fitak2018(file.path(truth_dir, "wolf_mexican_fitak2018"))

g <- w$genotypes
grp <- w$groups
n_markers <- 500

## ---- SNPs fixed between the two lineages, typed in every reference wolf ----
r1 <- grp == "ref1"; r2 <- grp == "ref2"
f1 <- colMeans(g[r1, ], na.rm = TRUE) / 2
f2 <- colMeans(g[r2, ], na.rm = TRUE) / 2
typed <- colSums(is.na(g[r1 | r2, ])) == 0
fixed <- which(typed & abs(f2 - f1) == 1)

## ---- spread them over the genome: each chromosome gets a share in
## proportion to its length, picked at evenly spaced positions ----
m <- w$map
chr_len <- tapply(m$position, m$chromosome, max)
pick <- integer(0)
for (ch in names(chr_len)) {
  cand <- fixed[m$chromosome[fixed] == as.integer(ch)]
  if (!length(cand)) next
  k <- min(length(cand), max(1, round(n_markers * chr_len[[ch]] / sum(chr_len))))
  target <- (seq_len(k) - 0.5) / k * chr_len[[ch]]
  for (t in target) {
    left <- setdiff(cand, pick)
    pick <- c(pick, left[which.min(abs(m$position[left] - t))])
  }
}
pick <- pick[order(m$chromosome[pick], m$position[pick])]

## ---- count Ghost Ranch (species 2) alleles ----
geno <- g[, pick]
flip <- f2[pick] == 0
geno[, flip] <- 2L - geno[, flip]
storage.mode(geno) <- "integer"
stopifnot(all(colMeans(geno[r1, ], na.rm = TRUE) == 0), all(colMeans(geno[r2, ], na.rm = TRUE) == 2))

map <- data.frame(snp = colnames(geno), chromosome = m$chromosome[pick],
                  position = round(m$cM_approx[pick], 3), bp = as.integer(m$position[pick]),
                  stringsAsFactors = FALSE)

## ---- per-animal truth: class and hybrid age only ----
meta <- w$meta
stopifnot(identical(meta$id, rownames(geno)))
cls <- c(ref1 = "pure1", ref2 = "pure2", F1 = "F1", BC1_MB = "BC1_species1", BC2_MB = "BC2_species1")
animals <- data.frame(id = meta$id, group = meta$group, class = unname(cls[meta$group]),
                      hybrid_age = ifelse(meta$group %in% c("ref1", "ref2"), 0, meta$truth_age_both),
                      stringsAsFactors = FALSE)
stopifnot(!anyNA(animals$class), !anyNA(animals$hybrid_age))

lev <- c("ref1", "ref2", "F1", "BC1_MB", "BC2_MB")
truth <- data.frame(
  group = lev,
  n = as.integer(table(factor(grp, levels = lev))),
  class = unname(cls[lev]),
  description = c("pure McBride lineage", "pure Ghost Ranch lineage",
                  "McBride x Ghost Ranch first-generation crosses",
                  "first backcrosses (F1 x pure McBride)",
                  "second backcrosses (first backcross x pure McBride)"),
  hybrid_age = c(0, 0, 1, 1.5, 1.75),
  stringsAsFactors = FALSE)
stopifnot(all(tapply(animals$hybrid_age, factor(animals$group, levels = lev), mean) == truth$hybrid_age))

wolves <- list(
  genotypes = geno,
  groups = grp,
  map = map,
  truth = truth,
  animals = animals,
  description = paste(
    "Captive Mexican wolves (Canis lupus baileyi): pure McBride (ref1) and Ghost Ranch (ref2)",
    "lineage wolves and their crosses (F1, first and second backcrosses to McBride), bred",
    "after the lineages were merged in 1995. Genotypes count Ghost Ranch alleles at",
    ncol(geno), "SNPs fixed between the lineages, spread over the 38 autosomes.",
    "Fitak et al. 2018, J Hered, doi:10.1093/jhered/esy009; Dryad doi:10.5061/dryad.g68k008 (CC0)."))

usethis::use_data(wolves, overwrite = TRUE, compress = "xz")
