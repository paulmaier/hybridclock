# Builds data/bettongs.rda: burrowing bettongs (boodies) from two subspecies
# released together into a fenced reserve in 2010, and their offspring
# sampled each year from 2010 to 2018.
#
# Source: Rick K, Ottewell K, Lohr C, Thavornkanlapachai R, Byrne M,
# Kennington WJ (2019) Genes 10:851, doi:10.3390/genes10110851.
# Data: Mendeley Data doi:10.17632/vg983472v3.1 (CC BY 4.0). Founder source and
# sampling year of each animal come from Table S6 of the paper (CC BY 4.0).
#
# The genotypes and groups come from the truth-set loader (load.R and TRUTH.md
# in the folder below), which reads the Stacks VCF and Table S6. No locality
# data are kept.
#
# Run from the package directory: Rscript data-raw/bettongs.R

truth_dir <- path.expand(Sys.getenv("HYBRIDCLOCK_TRUTH_DATA", "~/scripts/hybridclock_MEE/truth_data"))
source(file.path(truth_dir, "bettong_rick2019", "load.R"))
b <- load_bettong_rick2019(truth_dir)

g <- b$genotypes
grp <- b$groups
n_markers <- 500

## ---- SNPs fixed between the founder groups (typed in >= 5 of each) ----
r1 <- grp == "ref1"; r2 <- grp == "ref2"
f1 <- colMeans(g[r1, ], na.rm = TRUE) / 2
f2 <- colMeans(g[r2, ], na.rm = TRUE) / 2
typed <- colSums(!is.na(g[r1, ])) >= 5 & colSums(!is.na(g[r2, ])) >= 5
fixed <- which(typed & abs(f2 - f1) >= 0.99)
set.seed(2019)
pick <- sort(sample(fixed, n_markers))

## ---- count Dryandra / Shark Bay (species 2) alleles ----
geno <- g[, pick]
flip <- f2[pick] < 0.5
geno[, flip] <- 2L - geno[, flip]
storage.mode(geno) <- "integer"
stopifnot(all(colMeans(geno[r1, ], na.rm = TRUE) == 0), all(colMeans(geno[r2, ], na.rm = TRUE) == 2))

## ---- groups, per-animal labels and truth ----
lev <- c("ref1", "ref2", paste0("Matuwa", 2010:2018))
meta <- b$meta
stopifnot(identical(meta$id, rownames(geno)))
o <- order(match(grp, lev), meta$id)
geno <- geno[o, ]; grp <- grp[o]; meta <- meta[o, ]

animals <- data.frame(id = meta$id, group = meta$group, year_sampled = meta$year_sampled,
                      published_class = meta$newhybrids_class, stringsAsFactors = FALSE)
rownames(animals) <- NULL

yr <- c(NA, NA, 2010:2018)
truth <- data.frame(
  group = lev,
  n = as.integer(table(factor(grp, levels = lev))),
  year = as.integer(yr),
  years_since_mixing = as.integer(yr - 2010),
  max_hybrid_age = as.integer(yr - 2010),
  generations_at_3yr = round((yr - 2010) / 3, 2),
  stringsAsFactors = FALSE)
truth$max_hybrid_age[1:2] <- 0L

bettongs <- list(
  genotypes = geno,
  groups = grp,
  truth = truth,
  animals = animals,
  description = paste(
    "Burrowing bettongs (Bettongia lesueur): founders from Barrow Island (ref1) and from the",
    "Shark Bay subspecies via Dryandra (ref2), released together into a closed, fenced reserve",
    "in 2010, and animals born there, grouped by year of first capture (2010-2018). Genotypes",
    "count Shark Bay alleles at", ncol(geno), "unmapped RAD SNPs fixed between the founder groups;",
    "about a quarter of heterozygotes are called homozygous. Rick et al. 2019, Genes,",
    "doi:10.3390/genes10110851; Mendeley Data doi:10.17632/vg983472v3.1 (CC BY 4.0)."))

usethis::use_data(bettongs, overwrite = TRUE, compress = "xz")
