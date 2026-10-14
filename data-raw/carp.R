# Builds data/carp.rda: bighead x silver carp hybridizing in the Mississippi
# River Basin, with native-range reference fish from China.
#
# Sources: Lamer JT, Ruebush BC, Arbieva ZH, McClelland MA, Epifanio JM, Sass GG
# (2015) Molecular Ecology 24:3931, doi:10.1111/mec.13285; Dryad
# doi:10.5061/dryad.kp6j2 (CC0). Reference fish: Lamer JT, Sass GG, Boone JQ,
# Arbieva ZH, Green SJ, Epifanio JM (2013) Molecular Ecology Resources,
# doi:10.1111/1755-0998.12152; Dryad doi:10.5061/dryad.5j2v6 (CC0).
#
# The genotypes and groups come from the truth-set loader (load.R and TRUTH.md
# in the folder below), which reads both Dryad packages. Three of the nine
# river sites are kept, one from each of three river reaches.
#
# Run from the package directory: Rscript data-raw/carp.R

truth_dir <- path.expand(Sys.getenv("HYBRIDCLOCK_TRUTH_DATA", "~/scripts/hybridclock_MEE/truth_data"))
source(file.path(truth_dir, "carp_lamer2015", "load.R"))
cp <- load_carp_lamer2015(truth_dir, include_selected = FALSE)

lev <- c("ref1", "ref2", "IMAR", "PL26", "MKY")
keep <- cp$groups %in% lev
o <- which(keep)[order(match(cp$groups[keep], lev), which(keep))]
geno <- cp$genotypes[o, ]
storage.mode(geno) <- "integer"
grp <- cp$groups[o]
meta <- cp$meta[o, ]
stopifnot(identical(meta$id, rownames(geno)))
stopifnot(all(colMeans(geno[grp == "ref1", ], na.rm = TRUE) == 0),
          all(colMeans(geno[grp == "ref2", ], na.rm = TRUE) == 2))

cls <- ifelse(meta$published_class %in% "Hybrid", NA, meta$published_class)
animals <- data.frame(id = meta$id, group = grp, published_class = cls,
                      birth_year = as.integer(meta$birth_year), stringsAsFactors = FALSE)
rownames(animals) <- NULL

truth <- data.frame(
  group = lev,
  n = as.integer(table(factor(grp, levels = lev))),
  site = c("Pearl and Yangtze rivers, China (native range): bighead carp",
           "Pearl and Yangtze rivers, China (native range): silver carp",
           "Illinois River, Marseilles Reach: the invasion front",
           "Upper Mississippi River, Pool 26",
           "Lower Mississippi River, Kentucky"),
  sampled = c(NA, NA, "2009-2011", "2009-2011", "2009-2011"),
  max_hybrid_age = c(0L, 0L, 12L, 12L, 12L),
  stringsAsFactors = FALSE)

carp <- list(
  genotypes = geno,
  groups = grp,
  truth = truth,
  animals = animals,
  description = paste(
    "Bighead carp (Hypophthalmichthys nobilis, species 1) and silver carp (H. molitrix,",
    "species 2): native-range references from China and fish from three invaded sites in",
    "the Mississippi River Basin, sampled 2009-2011. Genotypes count silver carp alleles at",
    ncol(geno), "species-diagnostic SNPs. Both species were imported in the early 1970s, so",
    "hybrids are at most about 12 generations old. Lamer et al. 2015, Mol Ecol,",
    "doi:10.1111/mec.13285; Dryad doi:10.5061/dryad.kp6j2 and doi:10.5061/dryad.5j2v6 (CC0)."))

usethis::use_data(carp, overwrite = TRUE, compress = "xz")
