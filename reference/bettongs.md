# Burrowing bettongs: two subspecies mixed in a fenced reserve

Burrowing bettongs, or boodies (*Bettongia lesueur*), from two
subspecies were released together into a predator-proof fenced reserve
in Western Australia in 2010: animals from Barrow Island and Shark Bay
animals from a breeding facility at Dryandra. The reserve is closed, so
no animals arrived later. The founders and animals born in the reserve,
sampled each year from 2010 to 2018, were genotyped by RAD sequencing.

## Usage

``` r
bettongs
```

## Format

A list with:

- genotypes:

  Integer matrix, 215 animals x 500 SNPs, counting Shark Bay (species 2)
  alleles (0, 1, 2; `NA` missing). The SNPs are a random 500 of the
  1,233 that are fixed between the two founder groups, each on a
  separate RAD locus. They are not mapped.

- groups:

  `"ref1"` (Barrow Island founders, 22), `"ref2"` (Shark Bay
  founders, 24) and `"Matuwa2010"` to `"Matuwa2018"` (animals born in
  the reserve, by year of first capture; 4 to 37 per year).

- truth:

  Data frame, one row per group: `group`, `n`, `year`,
  `years_since_mixing`, `max_hybrid_age` (an upper bound) and
  `generations_at_3yr`.

- animals:

  Data frame, one row per animal: `id`, `group`, `year_sampled` and
  `published_class` (the NewHybrids class published by Rick et al. 2019,
  Table S6: `"BWI"`, `"DRY"`, `"F1"`, `"F2"`, `"F1BWI"`, `"F1DRY"`, or
  founder).

- description:

  A short description.

## Source

Rick K (2019) SNP dataset and Capture-Mark-Recapture dataset used in the
analysis for the paper Rick et al.: Population genomics of Bettongia
lesueur: admixing increases genetic diversity with no evidence of
outbreeding depression. Mendeley Data, version 1,
[doi:10.17632/vg983472v3.1](https://doi.org/10.17632/vg983472v3.1) (file
"SNP dataset of Bettongia lesueur at Matuwa from STACKS pipeline.vcf").
Licensed under CC BY 4.0; a subset of SNPs was taken and genotypes
recoded (see data-raw/bettongs.R).

## Truth

Mixing started in 2010. Bettongs can breed in their first year, so an
animal first caught in year Y is at most Y - 2010 generations from the
pure founders (`max_hybrid_age`); this upper bound is the only firm
fact. With the average generation length of three years, the expected
age is about (Y - 2010) / 3 (`generations_at_3yr`). Pure animals of both
sources kept breeding inside the fence, so new F1s and first backcrosses
appear every year, and most hybrids are much younger than the bound.
Year is the year of first capture, not of birth.

About a quarter of heterozygotes are called homozygous (F1s show H of
about 0.75 instead of 1), so use `het_dropout = 0.25` in
[`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md).
The published classes are model-based, not known pedigrees; a few
disagree with the genotypes (D6056 and D6057 are listed as F1 but are
pure Barrow Island; D5893 is listed as Barrow Island but has a hybrid
index of 0.38).

## Licence

The genotypes are from Mendeley Data and were released under CC BY 4.0,
which requires attribution: cite Rick et al. (2019) and the dataset (see
Source). The founder source and sampling year of each animal, and the
published classes, come from Table S6 of Rick et al. (2019), an open
access article under CC BY 4.0. Changes made: 500 SNPs were selected,
genotypes were recoded to count Shark Bay alleles, and only the animal
labels listed above were kept. No locality data are included.

## References

Rick K, Ottewell K, Lohr C, Thavornkanlapachai R, Byrne M, Kennington WJ
(2019) Population genomics of *Bettongia lesueur*: admixing increases
genetic diversity with no evidence of outbreeding depression. Genes
10:851.
[doi:10.3390/genes10110851](https://doi.org/10.3390/genes10110851)

## Examples

``` r
mk <- diagnostic_markers(bettongs$genotypes, bettongs$groups, "ref1", "ref2")
sc <- hybrid_scores(bettongs$genotypes, mk)
cl <- classify_hybrids(sc, het_dropout = 0.25)
table(bettongs$groups, cl$class)
#>             
#>              BC1_species1 BC1_species2 BC2_species1 F1 F2 later pure1 pure2
#>   Matuwa2010            0            0            0  0  0     0     4     0
#>   Matuwa2011            0            0            0  2  0     1     9     3
#>   Matuwa2012            0            1            0  7  0     3     5     3
#>   Matuwa2013            2            1            0  7  1    10     9     7
#>   Matuwa2014            5            1            1  4  1    14     6     1
#>   Matuwa2015            1            1            1  3  1     7     2     0
#>   Matuwa2016            0            0            0  1  0    12     2     0
#>   Matuwa2017            4            1            0  1  0     5     4     0
#>   Matuwa2018            3            0            2  1  1     6     2     0
#>   ref1                  0            0            0  0  0     0    22     0
#>   ref2                  0            0            0  0  0     0     0    24
```
