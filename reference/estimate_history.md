# Estimate the hybrid history of sampled populations

Gives one answer per population, from one of two tiers:

- **Individual classes.** When most hybrids confidently fit an early
  class (F1, F2, first or second backcross; see
  [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md))
  and those are mostly F1s and backcrosses rather than F2-like animals,
  the ages follow directly from the classes: F1 = 1, first backcross =
  1.5, second backcross = 1.75, F2 = 2 generations (hybrid age, both
  species), and likewise for the rare species' lines and the number of
  backcrosses. Intervals come from resampling animals, widened to at
  least the error of class-based ages in the reference at a similar
  share of early hybrids.

- **Population clock.** Otherwise the trained clock is used: the
  sample's summaries are compared with the simulated reference.

## Usage

``` r
estimate_history(
  clock,
  scores,
  population,
  diag = NULL,
  shared = NULL,
  markers = NULL,
  positions = NULL,
  effective_markers = NULL,
  het_dropout = 0,
  repeats = 30,
  seed = 1
)
```

## Arguments

- clock:

  Output of
  [`train_clock()`](https://paulmaier.github.io/hybridclock/reference/train_clock.md).

- scores:

  Output of
  [`hybrid_scores()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores.md)
  (columns `S`, `H`, `n_markers`).

- population:

  Population label of each animal.

- diag:

  Genotypes at the diagnostic markers
  ([`orient_genotypes()`](https://paulmaier.github.io/hybridclock/reference/orient_genotypes.md)),
  needed if the reference used them (it always does).

- shared:

  Genotypes at shared loci, needed if the reference had shared loci.

- markers, effective_markers, het_dropout:

  No longer used, kept so that older scripts run. The classes that
  decide between the two tiers are made exactly as in the reference
  (from S, H and the number of markers, with default settings), so that
  the rule learned on the reference applies; use
  [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md)
  with these settings for the individual calls.

- positions:

  Marker positions matching the columns of `diag` (`chromosome`,
  `position`), needed if the reference used a fixed map.

- repeats:

  Number of random subsets for larger populations.

- seed:

  Random seed.

## Value

A data.frame with one row per population: `n`, `n_hybrids` (animals not
classified as pure), `rare_species` (the species with less ancestry in
the hybrids), `rare_share` (its mean share in the hybrids), `method`,
`early_share` (share of hybrids in early classes), an estimate and
interval (`_lower`, `_upper`) for `hybrid_age_both`, `hybrid_age_rare`
and `backcrosses`, the arrival probability `arriving` and rate
`arrival_rate` (always from the clock), and notes. `arriving_note` flags
young swarms (hybrid age under 5 generations), where pure founders of
the rarer species still form new hybrids without any arrivals and the
arrival call was right only 63-67% of the time in validation.
`information` is 1 minus the width of the hybrid-age interval relative
to the spread of ages in the reference (log scale): near 1 the data pin
the age down, near 0 the answer mostly reflects the reference's age
range (few or weak markers), and `note` says so. `note` also flags
populations with too few hybrids.

## Details

F2-like animals are not trusted as second-generation hybrids on their
own: in a swarm that mates among itself, F3, F4 and later generations
have the same expected hybrid index and heterozygosity as F2s, so one
animal cannot show its generation. Only the population can, through
linkage blocks that shrink and drift that accumulates over generations,
which is what the clock reads.

The share of hybrids in early classes at which the answer switches to
the classes was chosen when the clock was trained, as the point above
which the classes' large errors (90th percentile) on held-out simulated
histories were no worse than the clock's, or under about 30%. With few
or weak markers, later-generation animals still fit early classes, so
the classes are never used. Column `method` says which tier gave the
answer; the other tier's estimates are kept in columns starting `clock_`
and `classes_`.

Populations sampled more heavily than the reference are summarized in
random subsets of the reference sample size and the results averaged (on
the model's scale), so that every summary is computed exactly as in
training.

## Examples

``` r
# \donttest{
ex <- hybridclock_example()
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
sc <- hybrid_scores(ex$genotypes, mk)
ref <- build_reference(n_sims = 150, generations = 12, n_markers = nrow(mk),
                       breeders = c(20, 100))
clock <- train_clock(ref, targets = "hybrid_age_both", committees = 5)
hyb <- ex$groups %in% c("site_A", "site_B")
estimate_history(clock, sc[hyb, ], ex$groups[hyb], diag = orient_genotypes(ex$genotypes, mk)[hyb, ])
#>   population  n n_hybrids rare_species rare_share           method early_share
#> 1     site_A 30        30     species2  0.4427813 population clock   0.5666667
#> 2     site_B 30        30     species2  0.1213907 population clock   0.1000000
#>   f2_share hybrid_age_both hybrid_age_both_lower hybrid_age_both_upper note
#> 1        1        4.066974              2.667928              5.859943     
#> 2        0        5.867665              3.774258              8.599688     
#>   classes_hybrid_age_both classes_hybrid_age_both_lower
#> 1                    2.00                      1.809675
#> 2                    1.75                      1.583465
#>   classes_hybrid_age_both_upper clock_hybrid_age_both
#> 1                      2.210342              4.066974
#> 2                      1.934049              5.867665
#>   clock_hybrid_age_both_lower clock_hybrid_age_both_upper information
#> 1                    2.667928                    5.859943   0.6529643
#> 2                    3.774258                    8.599688   0.6367845
# }
```
