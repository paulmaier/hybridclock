# Mountain pygmy possums: genetic rescue at 24 microsatellites

The mountain pygmy possum (*Burramys parvus*) population at Mount Buller
was small and inbred. In September 2011 six males from Mount
Higginbotham, in the central part of the range (isolated from Mount
Buller for at least 20,000 years), were released there, and four of them
sired young that season. More males, from Timms Spur, were released in
September 2014. Possums were trapped every spring and genotyped at 24
microsatellites.

## Usage

``` r
possums
```

## Format

A list with:

- alleles:

  Integer matrix, 536 rows x 48 columns: allele sizes at 24
  microsatellites, two columns per locus (`locus.a1`, `locus.a2`; `NA`
  missing), as used by
  [`microsat_markers()`](https://paulmaier.github.io/hybridclock/reference/microsat_markers.md).

- groups:

  `"ref1"` (Mount Buller 2010, before the rescue, 29), `"ref2"` (Mount
  Higginbotham 2012, the source population, 104), `"Buller2011"` to
  `"Buller2015"` (Mount Buller by year; 40, 43, 67, 103 and 138),
  `"Higginbotham_males_2011"` (the six males released in 2011) and
  `"TimmsSpur_males_2014"` (the six released in 2014).

- truth:

  Data frame, one row per group: `group`, `n`, `year`,
  `years_since_rescue`, `hybrid_age` (expected mean hybrid age of the
  hybrids, where known), `max_hybrid_age` (upper bound) and `note`.

- animals:

  Data frame, one row per row of `alleles`: `id`, `tag`, `group`, `year`
  and `new_this_year` (first year the tag was caught at Mount Buller).

- description:

  A short description.

## Source

Data package of Weeks et al. (2017) on Dryad,
[doi:10.5061/dryad.7g988](https://doi.org/10.5061/dryad.7g988) (file
`Burramys_parvus_genotypic_data_Mt_Buller.xlsx`). CC0 1.0.

## Truth

Possums breed as one-year-olds, so the oldest possible hybrid generation
equals the number of years since the 2011 breeding season. Buller 2011
was trapped as the first males were released and is pure. In 2012 every
hybrid is an F1 (hybrid age 1; the authors found 13). In 2013 the new
adults were F2s and backcrosses with no new F1s (6 F2 and 23
backcrosses), and surviving F1s returned, so the hybrids average about
1.5 generations. In 2014 hybrids are at most 3 generations old, but new
F1s also appeared (one 2011 male first bred in 2013, and captive-bred
animals with some Higginbotham ancestry were released after the 2013
trapping); in 2015 they are at most 4 generations old, with new F1s from
the Timms Spur males. The mean generation time is 1.54 years. Animals
caught in several years appear once per year, so year groups overlap.
Alleles from Timms Spur that Mount Higginbotham lacks appear in the 2015
sample.

## Licence

The genotypes are from Dryad and were released under CC0 1.0 (public
domain dedication). Rows were reordered and labelled; allele sizes are
unchanged. No locality data are included.

## References

Weeks AR, Heinze D, Perrin L, Stoklosa J, Hoffmann AA, van Rooyen A,
Kelly T, Mansergh I (2017) Genetic rescue increases fitness and aids
rapid recovery of an endangered marsupial population. Nature
Communications 8:1071.
[doi:10.1038/s41467-017-01182-3](https://doi.org/10.1038/s41467-017-01182-3)

## Examples

``` r
ms <- microsat_markers(possums$alleles, possums$groups, "ref1", "ref2")
y12 <- possums$groups == "Buller2012"
table(classify_hybrids_microsat(possums$alleles[y12, ], ms)$class)
#> 
#>    F1 pure1 
#>    13    30 
sc <- hybrid_scores_microsat(possums$alleles, ms)
aggregate(cbind(S, H) ~ group, data.frame(sc, group = possums$groups), mean)
#>                     group         S           H
#> 1              Buller2011 0.0005000 0.001000000
#> 2              Buller2012 0.1451163 0.284651163
#> 3              Buller2013 0.1692537 0.300000000
#> 4              Buller2014 0.1951456 0.322330097
#> 5              Buller2015 0.2617391 0.347246377
#> 6 Higginbotham_males_2011 0.9766667 0.020000000
#> 7    TimmsSpur_males_2014 0.9233333 0.063333333
#> 8                    ref1 0.0000000 0.000000000
#> 9                    ref2 0.9909615 0.009423077
```
