# Builds data/possums.rda: genetic rescue of the Mount Buller mountain pygmy
# possum population with males from the central region, at 24 microsatellites.
#
# Source: Weeks AR, Heinze D, Perrin L, et al. (2017) Nature Communications
# 8:1071, doi:10.1038/s41467-017-01182-3. Data: Dryad doi:10.5061/dryad.7g988 (CC0).
#
# The raw allele sizes and groups come from the truth-set loader (load.R and
# TRUTH.md in the folder below). All 24 loci are kept as allele sizes; no
# biallelic recoding is applied. No locality data are kept.
#
# Run from the package directory: Rscript data-raw/possums.R

truth_dir <- path.expand(Sys.getenv("HYBRIDCLOCK_TRUTH_DATA", "~/scripts/hybridclock_MEE/truth_data"))
source(file.path(truth_dir, "possum_weeks2017", "load.R"))
p <- load_possum_weeks2017(truth_dir)

lev <- c("ref1", "ref2", paste0("Buller", 2011:2015), "Higginbotham_males_2011", "TimmsSpur_males_2014")
o <- order(match(p$groups, lev), p$meta$year, p$meta$id)
alleles <- p$alleles[o, ]
storage.mode(alleles) <- "integer"
grp <- p$groups[o]
meta <- p$meta[o, ]
stopifnot(identical(meta$row_id, rownames(alleles)), length(unique(sub("\\.a[12]$", "", colnames(alleles)))) == 24)

animals <- data.frame(id = meta$row_id, tag = meta$id, group = grp, year = meta$year,
                      new_this_year = meta$new_this_year, stringsAsFactors = FALSE)
rownames(animals) <- NULL

truth <- data.frame(
  group = lev,
  n = as.integer(table(factor(grp, levels = lev))),
  year = c(2010L, 2012L, 2011:2015, 2011L, 2014L),
  years_since_rescue = c(NA, NA, 0:4, NA, NA),
  hybrid_age = c(0, 0, 0, 1, 1.5, NA, NA, 0, 0),
  max_hybrid_age = c(0L, 0L, 0L, 1L, 2L, 3L, 4L, 0L, 0L),
  note = c(
    "Mount Buller before the rescue",
    "Mount Higginbotham, central region: the source population",
    "sampled as the first males were released; no hybrid young yet",
    "hybrids are F1s from the four 2011 males (published: 13 F1)",
    "new adults are F2s and backcrosses, no new F1s (published: 6 F2, 23 backcrosses); returning F1s",
    "at most 3 generations; new F1s from a fifth 2011 male and young of captive-bred animals released in 2013",
    "at most 4 generations; new F1s from the Timms Spur males released in September 2014",
    "males released at Mount Buller in September 2011",
    "males released at Mount Buller in September 2014"),
  stringsAsFactors = FALSE)

possums <- list(
  alleles = alleles,
  groups = grp,
  truth = truth,
  animals = animals,
  description = paste(
    "Mountain pygmy possums (Burramys parvus): the small, inbred Mount Buller population before",
    "(2010, ref1) and after (2011-2015) males from Mount Higginbotham (ref2, central region)",
    "were released in September 2011, with the released males and a second group of males",
    "released in 2014. Raw allele sizes at 24 microsatellites; animals caught in several",
    "years appear once per year. Weeks et al. 2017, Nat Commun, doi:10.1038/s41467-017-01182-3;",
    "Dryad doi:10.5061/dryad.7g988 (CC0)."))

usethis::use_data(possums, overwrite = TRUE, compress = "xz")
