# Bighead and silver carp: an invasive hybrid zone

Bighead carp (*Hypophthalmichthys nobilis*) and silver carp (*H.
molitrix*) were imported from China to the United States in the early
1970s, escaped, and now hybridize widely in the Mississippi River Basin.
Native-range fish from China serve as references, with fish from three
invaded sites sampled in 2009-2011, typed at 57 species-diagnostic SNPs.

## Usage

``` r
carp
```

## Format

A list with:

- genotypes:

  Integer matrix, 1,354 fish x 57 SNPs, counting silver carp (species 2)
  alleles (0, 1, 2; `NA` missing). The SNPs are on separate assays and
  are not mapped.

- groups:

  `"ref1"` (bighead carp from the Pearl and Yangtze rivers, 16),
  `"ref2"` (silver carp from the same rivers, 20), `"IMAR"` (Illinois
  River, Marseilles Reach, the invasion front, 545), `"PL26"` (upper
  Mississippi River, Pool 26, 365) and `"MKY"` (lower Mississippi River,
  Kentucky, 408).

- truth:

  Data frame, one row per group: `group`, `n`, `site`, `sampled` and
  `max_hybrid_age` (an upper bound).

- animals:

  Data frame, one row per fish: `id`, `group`, `published_class` (the
  class published by Lamer et al. 2015 from NewHybrids and manual
  checks: `"BH"`, `"SV"`, `"F1"`, `"F2"`, backcrosses `"BxBH"` to
  `"Bx4SV"`, and later `"FxBH"`, `"FxSV"`) and `birth_year` (capture
  year minus the age read from the postcleithrum bone).

- description:

  A short description.

## Source

Data package of Lamer et al. (2015) on Dryad,
[doi:10.5061/dryad.kp6j2](https://doi.org/10.5061/dryad.kp6j2) (invaded
sites). CC0 1.0.

Data package of Lamer et al. (2013) on Dryad,
[doi:10.5061/dryad.5j2v6](https://doi.org/10.5061/dryad.5j2v6)
(native-range reference fish). CC0 1.0.

## Truth

Both species were imported in 1973 and were in open waters by the early
1980s. With a generation time of at least 3 years, hybrids sampled in
2009-2011 are at most about 12 generations old (37 years / 3), and
probably 6-12; each site was colonized later than the basin as a whole,
so its own hybrids are likely younger. Hybridization is ongoing: F1s
were born as late as 2006. The published classes are inferred, not known
pedigrees. The references are partly circular, because some Yangtze fish
were used to choose the SNPs, so they look perfectly fixed. Nets
targeted fish over 400 mm, so recent cohorts are under-represented.

## Licence

The genotypes are from two Dryad packages, released under CC0 1.0
(public domain dedication). Three of the nine invaded sites were kept,
and genotypes were recoded to count silver carp alleles.

## References

Lamer JT, Ruebush BC, Arbieva ZH, McClelland MA, Epifanio JM, Sass GG
(2015) Diagnostic SNPs reveal widespread introgressive hybridization
between introduced bighead and silver carp in the Mississippi River
Basin. Molecular Ecology 24:3931.
[doi:10.1111/mec.13285](https://doi.org/10.1111/mec.13285)

Lamer JT, Sass GG, Boone JQ, Arbieva ZH, Green SJ, Epifanio JM (2013)
Restriction site-associated DNA sequencing generates high-quality single
nucleotide polymorphisms for assessing hybridization between bighead and
silver carp in the United States and China. Molecular Ecology Resources.
[doi:10.1111/1755-0998.12152](https://doi.org/10.1111/1755-0998.12152)

## Examples

``` r
mk <- diagnostic_markers(carp$genotypes, carp$groups, "ref1", "ref2")
sc <- hybrid_scores(carp$genotypes, mk)
aggregate(cbind(S, H) ~ group, data.frame(sc, group = carp$groups), mean)
#>   group         S          H
#> 1  IMAR 0.6115810 0.10918435
#> 2   MKY 0.9752853 0.01700759
#> 3  PL26 0.4701654 0.05260538
#> 4  ref1 0.0000000 0.00000000
#> 5  ref2 1.0000000 0.00000000
table(carp$groups, classify_hybrids(sc)$class)
#>       
#>        BC1_species1 BC1_species2 BC2_species1 BC2_species2  F1  F2 later pure1
#>   IMAR           14            8           23           86  14   1    60   150
#>   MKY             0            0            1           11   0   0    25     6
#>   PL26            2            0            6           23   4   3    37   172
#>   ref1            0            0            0            0   0   0     0    16
#>   ref2            0            0            0            0   0   0     0     0
#>       
#>        pure2
#>   IMAR   189
#>   MKY    365
#>   PL26   118
#>   ref1     0
#>   ref2    20
```
