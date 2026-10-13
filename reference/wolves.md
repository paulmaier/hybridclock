# Mexican wolves: crosses between two captive lineages

Captive Mexican wolves (*Canis lupus baileyi*) from two lineages that
were kept apart for 15-20 years, McBride and Ghost Ranch, and wolves
from the crosses made after the lineages were merged in 1995:
first-generation crosses (F1), first backcrosses to McBride and second
backcrosses to McBride. A young hybrid population with a known pedigree
for every animal.

## Usage

``` r
wolves
```

## Format

A list with:

- genotypes:

  Integer matrix, 56 wolves x 486 SNPs, counting Ghost Ranch (species 2)
  alleles (0, 1, 2; `NA` missing). The SNPs are fixed between the two
  reference groups and typed in every reference wolf, and are spread
  over the 38 autosomes in proportion to their length.

- groups:

  Group of each wolf: `"ref1"` (pure McBride, 31), `"ref2"` (pure Ghost
  Ranch, 6), `"F1"` (5), `"BC1_MB"` (11) and `"BC2_MB"` (3).

- map:

  Data frame with `snp`, `chromosome`, `position` (cM, approximate: Mb x
  0.97, the genome-wide rate of the dog linkage map) and `bp` (CanFam2).

- truth:

  Data frame, one row per group: `group`, `n`, `class` (as named by
  [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md)),
  `description` and `hybrid_age`.

- animals:

  Data frame, one row per wolf: `id`, `group`, `class` and `hybrid_age`
  from the pedigree.

- description:

  A short description.

## Source

Data package of Fitak et al. (2018) on Dryad,
[doi:10.5061/dryad.g68k008](https://doi.org/10.5061/dryad.g68k008)
(files `MW.ped` and `MW.map`). CC0 1.0.

## Truth

Every wolf's class comes from the pedigree and was checked against its
genotypes. Hybrid age is the number of generations back to the nearest
pure ancestor, averaged over the pedigree lines: F1 = 1, first backcross
= 1.5, second backcross = 1.75. The sample labelled `"547"` is a
re-identified museum voucher and is an F1. The two lineages differ only
by drift from very few founders, the Ghost Ranch reference is only six
related wolves (including parents of the hybrids), and several hybrids
are full siblings. With at most two meioses since the F1, ancestry
blocks are long, so markers on one chromosome are far from independent
(see
[`effective_markers()`](https://paulmaier.github.io/hybridclock/reference/effective_markers.md)).

## Licence

The genotypes are from Dryad and were released under CC0 1.0 (public
domain dedication). They were subset to the markers and wolves described
above and recoded.

## References

Fitak RR, Rinkevich SE, Culver M (2018) Genome-wide analysis of SNPs is
consistent with no domestic dog ancestry in the endangered Mexican wolf
(*Canis lupus baileyi*). Journal of Heredity.
[doi:10.1093/jhered/esy009](https://doi.org/10.1093/jhered/esy009)

## Examples

``` r
mk <- diagnostic_markers(wolves$genotypes, wolves$groups, "ref1", "ref2")
sc <- hybrid_scores(wolves$genotypes, mk)
cl <- classify_hybrids(sc)
table(truth = wolves$animals$class, called = cl$class)
#>               called
#> truth          BC1_species1 BC2_species1 F1 pure1 pure2
#>   BC1_species1           11            0  0     0     0
#>   BC2_species1            0            3  0     0     0
#>   F1                      0            0  5     0     0
#>   pure1                   0            0  0    31     0
#>   pure2                   0            0  0     0     6
```
