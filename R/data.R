#' Mexican wolves: crosses between two captive lineages
#'
#' Captive Mexican wolves (*Canis lupus baileyi*) from two lineages that were
#' kept apart for 15-20 years, McBride and Ghost Ranch, and wolves from the
#' crosses made after the lineages were merged in 1995: first-generation
#' crosses (F1), first backcrosses to McBride and second backcrosses to
#' McBride. A young hybrid population with a known pedigree for every
#' animal.
#'
#' @format A list with:
#' \describe{
#'   \item{genotypes}{Integer matrix, 56 wolves x 486 SNPs, counting Ghost
#'     Ranch (species 2) alleles (0, 1, 2; `NA` missing). The SNPs are fixed
#'     between the two reference groups and typed in every reference wolf,
#'     and are spread over the 38 autosomes in proportion to their length.}
#'   \item{groups}{Group of each wolf: `"ref1"` (pure McBride, 31),
#'     `"ref2"` (pure Ghost Ranch, 6), `"F1"` (5), `"BC1_MB"` (11) and
#'     `"BC2_MB"` (3).}
#'   \item{map}{Data frame with `snp`, `chromosome`, `position` (cM,
#'     approximate: Mb x 0.97, the genome-wide rate of the dog linkage map)
#'     and `bp` (CanFam2).}
#'   \item{truth}{Data frame, one row per group: `group`, `n`, `class` (as
#'     named by [classify_hybrids()]), `description` and `hybrid_age`.}
#'   \item{animals}{Data frame, one row per wolf: `id`, `group`, `class` and
#'     `hybrid_age` from the pedigree.}
#'   \item{description}{A short description.}
#' }
#'
#' @section Truth:
#' Every wolf's class comes from the pedigree and was checked against its
#' genotypes. Hybrid age is the number of generations back to the nearest
#' pure ancestor, averaged over the pedigree lines: F1 = 1, first backcross
#' = 1.5, second backcross = 1.75. The sample labelled `"547"` is a
#' re-identified museum voucher and is an F1. The two lineages differ only
#' by drift from very few founders, the Ghost Ranch reference is only six
#' related wolves (including parents of the hybrids), and several hybrids are
#' full siblings. With at most two meioses since the F1, ancestry blocks are
#' long, so markers on one chromosome are far from independent (see
#' [effective_markers()]).
#'
#' @section Licence:
#' The genotypes are from Dryad and were released under CC0 1.0 (public
#' domain dedication). They were subset to the markers and wolves described
#' above and recoded.
#'
#' @source Data package of Fitak et al. (2018) on Dryad,
#'   \doi{10.5061/dryad.g68k008} (files `MW.ped` and `MW.map`). CC0 1.0.
#' @references Fitak RR, Rinkevich SE, Culver M (2018) Genome-wide analysis
#'   of SNPs is consistent with no domestic dog ancestry in the endangered
#'   Mexican wolf (*Canis lupus baileyi*). Journal of Heredity.
#'   \doi{10.1093/jhered/esy009}
#' @examples
#' mk <- diagnostic_markers(wolves$genotypes, wolves$groups, "ref1", "ref2")
#' sc <- hybrid_scores(wolves$genotypes, mk)
#' cl <- classify_hybrids(sc)
#' table(truth = wolves$animals$class, called = cl$class)
"wolves"

#' Chickens: an advanced intercross after eight generations
#'
#' Generation F8 of an advanced intercross line between Red Junglefowl and
#' White Leghorn chickens: one Red Junglefowl male crossed with three White
#' Leghorn females, then kept as a closed line of about 100 birds per
#' generation with no backcrossing. Every bird is eight generations from the
#' pure parents.
#'
#' @format A list with:
#' \describe{
#'   \item{genotypes}{Integer matrix, 572 birds x 514 autosomal SNPs,
#'     counting White Leghorn (species 2) alleles (0, 1, 2; `NA` missing).}
#'   \item{groups}{`"F8"` for every bird. No pure birds were genotyped.}
#'   \item{map}{Data frame with `snp`, `chromosome`, `position` (cM per
#'     meiosis, approximate: the authors' F8 map divided by 3.09, their
#'     genome-wide map expansion over eight generations) and `position_f8`
#'     (cM on the authors' F8 map).}
#'   \item{truth}{Data frame with `group`, `n`, `generation`, `hybrid_age`
#'     (8), `backcrosses` (0) and `breeders`.}
#'   \item{description}{A short description.}
#' }
#'
#' @section Truth:
#' All birds are F8, so the hybrid age is 8, with no backcrosses and no
#' arrivals. The line was kept at about 100 birds per generation from F3 to
#' F7 and produced the F8 from 122 F7 parents; the mating scheme is not
#' reported. The loss of heterozygosity (mean H about 0.47 instead of 0.5)
#' suggests about 50 effective breeders. The founders were outbred, so not
#' every marker was fixed between them, which also lowers H. Which allele
#' is White Leghorn is not stated by the authors; it was inferred from the Z
#' chromosome (excluded here), where White Leghorn alleles make up one third
#' of copies because of the direction of the founding cross.
#'
#' @section Licence:
#' The genotypes are from figshare and were released under CC0 1.0 (public
#' domain dedication). Markers with partial calls, on the Z chromosome or on
#' an unplaced scaffold were dropped, and genotypes were recoded.
#'
#' @source Johnsson M (2017) Open field phenotypes and genotypes (R/qtl csv
#'   format). figshare, \doi{10.6084/m9.figshare.1265060}. CC0 1.0.
#' @references Johnsson M, Williams MJ, Jensen P, Wright D (2016) Genetical
#'   genomics of behavior: a novel chicken genomic model for anxiety
#'   behavior. Genetics 202:327-340. \doi{10.1534/genetics.115.179010}
#'
#'   Johnsson M, Gustafson I, Rubin C-J, et al. (2012) A sexual ornament in
#'   chickens is affected by pleiotropic alleles at HAO1 and BMP2, selected
#'   during domestication. PLoS Genetics 8:e1002914.
#'   \doi{10.1371/journal.pgen.1002914}
#' @examples
#' # every SNP already counts White Leghorn alleles
#' g <- chickens$genotypes
#' mk <- data.frame(column = seq_len(ncol(g)), name = colnames(g), flip = FALSE,
#'                  freq1 = 0, freq2 = 1)
#' sc <- hybrid_scores(g, mk)
#' colMeans(sc[, c("S", "H")])
"chickens"

#' Burrowing bettongs: two subspecies mixed in a fenced reserve
#'
#' Burrowing bettongs, or boodies (*Bettongia lesueur*), from two subspecies
#' were released together into a predator-proof fenced reserve in Western
#' Australia in 2010: animals from Barrow Island and Shark Bay animals from
#' a breeding facility at Dryandra. The reserve is closed, so no animals
#' arrived later. The founders and animals born in the reserve, sampled each
#' year from 2010 to 2018, were genotyped by RAD sequencing.
#'
#' @format A list with:
#' \describe{
#'   \item{genotypes}{Integer matrix, 215 animals x 500 SNPs, counting Shark
#'     Bay (species 2) alleles (0, 1, 2; `NA` missing). The SNPs are a random
#'     500 of the 1,233 that are fixed between the two founder groups, each on
#'     a separate RAD locus. They are not mapped.}
#'   \item{groups}{`"ref1"` (Barrow Island founders, 22), `"ref2"` (Shark
#'     Bay founders, 24) and `"Matuwa2010"` to `"Matuwa2018"` (animals born
#'     in the reserve, by year of first capture; 4 to 37 per year).}
#'   \item{truth}{Data frame, one row per group: `group`, `n`, `year`,
#'     `years_since_mixing`, `max_hybrid_age` (an upper bound) and
#'     `generations_at_3yr`.}
#'   \item{animals}{Data frame, one row per animal: `id`, `group`,
#'     `year_sampled` and `published_class` (the NewHybrids class published
#'     by Rick et al. 2019, Table S6: `"BWI"`, `"DRY"`, `"F1"`, `"F2"`,
#'     `"F1BWI"`, `"F1DRY"`, or founder).}
#'   \item{description}{A short description.}
#' }
#'
#' @section Truth:
#' Mixing started in 2010. Bettongs can breed in their first year, so an
#' animal first caught in year Y is at most Y - 2010 generations from the
#' pure founders (`max_hybrid_age`); this upper bound is the only firm fact.
#' With the average generation length of three years, the expected age is
#' about (Y - 2010) / 3 (`generations_at_3yr`). Pure animals of both sources
#' kept breeding inside the fence, so new F1s and first backcrosses appear
#' every year, and most hybrids are much younger than the bound. Year is the
#' year of first capture, not of birth.
#'
#' About a quarter of heterozygotes are called homozygous (F1s show H of
#' about 0.75 instead of 1), so use `het_dropout = 0.25` in
#' [classify_hybrids()]. The published classes are model-based, not known
#' pedigrees; a few disagree with the genotypes (D6056 and D6057 are listed
#' as F1 but are pure Barrow Island; D5893 is listed as Barrow Island but
#' has a hybrid index of 0.38).
#'
#' @section Licence:
#' The genotypes are from Mendeley Data and were released under CC BY 4.0,
#' which requires attribution: cite Rick et al. (2019) and the dataset
#' (see Source). The founder source and sampling year of each animal, and
#' the published classes, come from Table S6 of Rick et al. (2019), an open
#' access article under CC BY 4.0. Changes made: 500 SNPs were selected,
#' genotypes were recoded to count Shark Bay alleles, and only the animal
#' labels listed above were kept. No locality data are included.
#'
#' @source Rick K (2019) SNP dataset and Capture-Mark-Recapture dataset used
#'   in the analysis for the paper Rick et al.: Population genomics of
#'   Bettongia lesueur: admixing increases genetic diversity with no evidence
#'   of outbreeding depression. Mendeley Data, version 1,
#'   \doi{10.17632/vg983472v3.1} (file "SNP dataset of Bettongia lesueur at
#'   Matuwa from STACKS pipeline.vcf"). Licensed under CC BY 4.0; a subset of
#'   SNPs was taken and genotypes recoded (see data-raw/bettongs.R).
#' @references Rick K, Ottewell K, Lohr C, Thavornkanlapachai R, Byrne M,
#'   Kennington WJ (2019) Population genomics of *Bettongia lesueur*:
#'   admixing increases genetic diversity with no evidence of outbreeding
#'   depression. Genes 10:851. \doi{10.3390/genes10110851}
#' @examples
#' mk <- diagnostic_markers(bettongs$genotypes, bettongs$groups, "ref1", "ref2")
#' sc <- hybrid_scores(bettongs$genotypes, mk)
#' cl <- classify_hybrids(sc, het_dropout = 0.25)
#' table(bettongs$groups, cl$class)
"bettongs"

#' Mountain pygmy possums: genetic rescue at 24 microsatellites
#'
#' The mountain pygmy possum (*Burramys parvus*) population at Mount Buller
#' was small and inbred. In September 2011 six males from Mount Higginbotham,
#' in the central part of the range (isolated from Mount Buller for at least
#' 20,000 years), were released there, and four of them sired young that
#' season. More males, from Timms Spur, were released in September 2014.
#' Possums were trapped every spring and genotyped at 24 microsatellites.
#'
#' @format A list with:
#' \describe{
#'   \item{alleles}{Integer matrix, 536 rows x 48 columns: allele sizes at
#'     24 microsatellites, two columns per locus (`locus.a1`, `locus.a2`;
#'     `NA` missing), as used by [microsat_markers()].}
#'   \item{groups}{`"ref1"` (Mount Buller 2010, before the rescue, 29),
#'     `"ref2"` (Mount Higginbotham 2012, the source population, 104),
#'     `"Buller2011"` to `"Buller2015"` (Mount Buller by year; 40, 43, 67,
#'     103 and 138), `"Higginbotham_males_2011"` (the six males released in
#'     2011) and `"TimmsSpur_males_2014"` (the six released in 2014).}
#'   \item{truth}{Data frame, one row per group: `group`, `n`, `year`,
#'     `years_since_rescue`, `hybrid_age` (expected mean hybrid age of the
#'     hybrids, where known), `max_hybrid_age` (upper bound) and `note`.}
#'   \item{animals}{Data frame, one row per row of `alleles`: `id`, `tag`,
#'     `group`, `year` and `new_this_year` (first year the tag was caught at
#'     Mount Buller).}
#'   \item{description}{A short description.}
#' }
#'
#' @section Truth:
#' Possums breed as one-year-olds, so the oldest possible hybrid generation
#' equals the number of years since the 2011 breeding season. Buller 2011
#' was trapped as the first males were released and is pure. In 2012 every
#' hybrid is an F1 (hybrid age 1; the authors found 13). In 2013 the new
#' adults were F2s and backcrosses with no new F1s (6 F2 and 23 backcrosses),
#' and surviving F1s returned, so the hybrids average about 1.5 generations.
#' In 2014 hybrids are at most 3 generations old, but new F1s also appeared
#' (one 2011 male first bred in 2013, and captive-bred animals with some
#' Higginbotham ancestry were released after the 2013 trapping); in 2015
#' they are at most 4 generations old, with new F1s from the Timms Spur
#' males. The mean generation time is 1.54 years. Animals caught in several
#' years appear once per year, so year groups overlap. Alleles from Timms
#' Spur that Mount Higginbotham lacks appear in the 2015 sample.
#'
#' @section Licence:
#' The genotypes are from Dryad and were released under CC0 1.0 (public
#' domain dedication). Rows were reordered and labelled; allele sizes are
#' unchanged. No locality data are included.
#'
#' @source Data package of Weeks et al. (2017) on Dryad,
#'   \doi{10.5061/dryad.7g988} (file
#'   `Burramys_parvus_genotypic_data_Mt_Buller.xlsx`). CC0 1.0.
#' @references Weeks AR, Heinze D, Perrin L, Stoklosa J, Hoffmann AA, van
#'   Rooyen A, Kelly T, Mansergh I (2017) Genetic rescue increases fitness
#'   and aids rapid recovery of an endangered marsupial population. Nature
#'   Communications 8:1071. \doi{10.1038/s41467-017-01182-3}
#' @examples
#' ms <- microsat_markers(possums$alleles, possums$groups, "ref1", "ref2")
#' y12 <- possums$groups == "Buller2012"
#' table(classify_hybrids_microsat(possums$alleles[y12, ], ms)$class)
#' sc <- hybrid_scores_microsat(possums$alleles, ms)
#' aggregate(cbind(S, H) ~ group, data.frame(sc, group = possums$groups), mean)
"possums"

#' Bighead and silver carp: an invasive hybrid zone
#'
#' Bighead carp (*Hypophthalmichthys nobilis*) and silver carp (*H.
#' molitrix*) were imported from China to the United States in the early
#' 1970s, escaped, and now hybridize widely in the Mississippi River Basin.
#' Native-range fish from China serve as references, with fish from three
#' invaded sites sampled in 2009-2011, typed at 57 species-diagnostic SNPs.
#'
#' @format A list with:
#' \describe{
#'   \item{genotypes}{Integer matrix, 1,354 fish x 57 SNPs, counting silver
#'     carp (species 2) alleles (0, 1, 2; `NA` missing). The SNPs are on
#'     separate assays and are not mapped.}
#'   \item{groups}{`"ref1"` (bighead carp from the Pearl and Yangtze rivers,
#'     16), `"ref2"` (silver carp from the same rivers, 20), `"IMAR"`
#'     (Illinois River, Marseilles Reach, the invasion front, 545), `"PL26"`
#'     (upper Mississippi River, Pool 26, 365) and `"MKY"` (lower Mississippi
#'     River, Kentucky, 408).}
#'   \item{truth}{Data frame, one row per group: `group`, `n`, `site`,
#'     `sampled` and `max_hybrid_age` (an upper bound).}
#'   \item{animals}{Data frame, one row per fish: `id`, `group`,
#'     `published_class` (the class published by Lamer et al. 2015 from
#'     NewHybrids and manual checks: `"BH"`, `"SV"`, `"F1"`, `"F2"`,
#'     backcrosses `"BxBH"` to `"Bx4SV"`, and later `"FxBH"`, `"FxSV"`) and
#'     `birth_year` (capture year minus the age read from the postcleithrum bone).}
#'   \item{description}{A short description.}
#' }
#'
#' @section Truth:
#' Both species were imported in 1973 and were in open waters by the early
#' 1980s. With a generation time of at least 3 years, hybrids sampled in
#' 2009-2011 are at most about 12 generations old (37 years / 3), and
#' probably 6-12; each site was colonized later than the basin as a whole,
#' so its own hybrids are likely younger. Hybridization is ongoing: F1s were
#' born as late as 2006. The published classes are inferred, not known
#' pedigrees. The references are partly circular, because some Yangtze fish
#' were used to choose the SNPs, so they look perfectly fixed. Nets targeted
#' fish over 400 mm, so recent cohorts are under-represented.
#'
#' @section Licence:
#' The genotypes are from two Dryad packages, released under CC0 1.0
#' (public domain dedication). Three of the nine invaded sites were kept,
#' and genotypes were recoded to count silver carp alleles.
#'
#' @source Data package of Lamer et al. (2015) on Dryad,
#'   \doi{10.5061/dryad.kp6j2} (invaded sites). CC0 1.0.
#'
#'   Data package of Lamer et al. (2013) on Dryad,
#'   \doi{10.5061/dryad.5j2v6} (native-range reference fish). CC0 1.0.
#' @references Lamer JT, Ruebush BC, Arbieva ZH, McClelland MA, Epifanio JM,
#'   Sass GG (2015) Diagnostic SNPs reveal widespread introgressive
#'   hybridization between introduced bighead and silver carp in the
#'   Mississippi River Basin. Molecular Ecology 24:3931.
#'   \doi{10.1111/mec.13285}
#'
#'   Lamer JT, Sass GG, Boone JQ, Arbieva ZH, Green SJ, Epifanio JM (2013)
#'   Restriction site-associated DNA sequencing generates high-quality single
#'   nucleotide polymorphisms for assessing hybridization between bighead and
#'   silver carp in the United States and China. Molecular Ecology
#'   Resources. \doi{10.1111/1755-0998.12152}
#' @examples
#' mk <- diagnostic_markers(carp$genotypes, carp$groups, "ref1", "ref2")
#' sc <- hybrid_scores(carp$genotypes, mk)
#' aggregate(cbind(S, H) ~ group, data.frame(sc, group = carp$groups), mean)
#' table(carp$groups, classify_hybrids(sc)$class)
"carp"
