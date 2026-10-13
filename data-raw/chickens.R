# Builds data/chickens.rda: generation F8 of a Red Junglefowl x White Leghorn
# advanced intercross.
#
# Source: Johnsson M, Williams MJ, Jensen P, Wright D (2016) Genetics 202:327-340,
# doi:10.1534/genetics.115.179010. Data: figshare doi:10.6084/m9.figshare.1265060
# (Martin Johnsson, CC0).
#
# The genotypes and map come from the truth-set loader (load.R and TRUTH.md in
# the folder below). It counts the B allele, which it shows to be the White
# Leghorn allele from the Z chromosome, drops the Z, the markers with partial
# calls and (here) the few markers on an unplaced chromosome-8 scaffold.
#
# Run from the package directory: Rscript data-raw/chickens.R

truth_dir <- path.expand(Sys.getenv("HYBRIDCLOCK_TRUTH_DATA", "~/scripts/hybridclock_analysis/truth_data"))
source(file.path(truth_dir, "chicken_rjfwl_ail", "load.R"))
ch <- load_chicken_rjfwl_ail(file.path(truth_dir, "chicken_rjfwl_ail"), drop_random = TRUE)

geno <- ch$genotypes
storage.mode(geno) <- "integer"
rownames(geno) <- paste0("bird", rownames(geno))

map <- data.frame(snp = ch$map$snp, chromosome = as.integer(ch$map$chromosome),
                  position = round(ch$map$position, 3), position_f8 = ch$map$position_f8map,
                  stringsAsFactors = FALSE)
stopifnot(identical(map$snp, colnames(geno)), !anyNA(map$chromosome))

truth <- data.frame(
  group = "F8", n = nrow(geno), generation = 8L, hybrid_age = 8, backcrosses = 0L,
  breeders = "about 100 birds per generation from F3 to F7; 122 F7 parents",
  stringsAsFactors = FALSE)

chickens <- list(
  genotypes = geno,
  groups = ch$groups,
  map = map,
  truth = truth,
  description = paste(
    "Advanced intercross line between Red Junglefowl (species 1) and White Leghorn (species 2)",
    "chickens, generation F8: one Red Junglefowl male x three White Leghorn females, about 100",
    "birds per generation, no backcrossing. Genotypes count White Leghorn alleles at",
    ncol(geno), "autosomal SNPs. Johnsson et al. 2016, Genetics, doi:10.1534/genetics.115.179010;",
    "figshare doi:10.6084/m9.figshare.1265060 (CC0)."))

usethis::use_data(chickens, overwrite = TRUE, compress = "xz")
