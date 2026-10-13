# Simulate a hybrid population forward in time

Simulates a population founded by two species that interbreed for a
number of generations, with optional ongoing arrival of pure animals of
either species, mating between siblings, and separate breeding pools.
Each generation, a sample of offspring (for example tadpoles) is drawn
from that generation's matings, scored at unlinked diagnostic markers,
and summarized with
[`sample_summaries()`](https://paulmaier.github.io/hybridclock/reference/sample_summaries.md).
The true history of every sampled animal is tracked alongside, so the
output can be used to train and test estimators.

## Usage

``` r
simulate_hybrid_zone(
  generations = 20,
  breeders = 100,
  start_mix = 0.5,
  arrivals = c(0, 0),
  sib_mating = 0,
  pools = 1,
  pool_crossing = 0.1,
  n_markers = 1000,
  shared_freq = NULL,
  sample_size = 20,
  missing = 0,
  sample_at = NULL,
  seed = NULL,
  founders = c("pure", "F1", "BC1_1", "BC1_2"),
  map = NULL,
  assort = 0,
  het_dropout = 0,
  junctions = FALSE,
  freq = NULL,
  microsat = NULL,
  founder_impurity = 0,
  sample_founders = FALSE,
  class_settings = NULL
)
```

## Arguments

- generations:

  Number of generations to simulate.

- breeders:

  Number of breeding adults per generation.

- start_mix:

  Share of species 2 among the founders (pure founders only).

- arrivals:

  Share of breeders replaced each generation by newly arrived pure
  animals of species 1 and species 2: a length-2 vector for constant
  rates, or a matrix with one row per generation and two columns for
  rates that change over time (arrivals that start or stop).

- sib_mating:

  Probability that an animal mates with a sibling (shared mother) when
  one is available.

- pools:

  Number of breeding pools. Animals mate within their pool except with
  probability `pool_crossing`.

- pool_crossing:

  Probability that a mating crosses pools.

- n_markers:

  Number of diagnostic markers.

- shared_freq:

  Allele frequencies of shared loci, or `NULL` for none.

- sample_size:

  Number of offspring sampled each generation.

- missing:

  Per-animal missing-data rates to draw from (a single value or a vector
  of observed rates).

- sample_at:

  Generations at which to return the sampled animals' genotypes
  (attribute `samples`).

- seed:

  Optional random seed. The global random state is restored afterwards.

- founders:

  Who founds the population: `"pure"` animals of both species (in
  proportions `start_mix`), or all `"F1"`, or all first backcrosses to
  species 1 (`"BC1_1"`) or species 2 (`"BC1_2"`), as in experimental and
  hatchery populations.

- map:

  Optional marker map for linkage: a data.frame with columns
  `chromosome` and `position` (cM), one row per diagnostic marker, as
  from
  [`random_map()`](https://paulmaier.github.io/hybridclock/reference/random_map.md).
  Markers are simulated in map order, and recombination between
  neighbours on a chromosome follows Haldane's map function. `NULL` (the
  default) treats all markers as unlinked. Shared loci are always
  unlinked.

- assort:

  Probability that a mating is assortative by ancestry: the mate is the
  closest in ancestry among five random candidates (for example through
  mate choice or different breeding times).

- het_dropout:

  Probability that a heterozygous genotype is called as homozygous
  (either allele), as with allele dropout in reduced- representation
  sequencing. Applied to the observed genotypes.

- junctions:

  Also summarize genotype switches along chromosomes (requires `map`;
  see
  [`sample_summaries()`](https://paulmaier.github.io/hybridclock/reference/sample_summaries.md)).
  Use this only when the observed data have marker positions too.

- freq:

  Optional allele frequencies for markers that are not fixed: a list or
  data.frame with `freq1` and `freq2`, the frequency of the counted
  allele at each diagnostic marker in species 1 and species 2 (as from
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md)).
  Pure animals then draw their alleles from these frequencies, and S and
  H are maximum-likelihood estimates (see
  [`hybrid_scores()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores.md)).

- microsat:

  Optional output of
  [`microsat_markers()`](https://paulmaier.github.io/hybridclock/reference/microsat_markers.md):
  the diagnostic markers are then microsatellites. Each copy's species
  of origin is simulated, observed alleles are drawn from that species'
  frequencies, and S and H are estimated with
  [`hybrid_scores_microsat()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores_microsat.md).
  `n_markers` is set to the number of loci.

- founder_impurity:

  Share of diagnostic alleles that founders carry from the other
  species, for founders that were not quite pure (for example, cross
  parents with some old introgression). Arrivals are always pure.

- sample_founders:

  With F1 or backcross founders, also sample the founders themselves
  (generation 0), so that the output covers families of F1s or first
  backcrosses as well as their descendants.

- class_settings:

  Settings of the classes in the summaries (see
  [`class_settings()`](https://paulmaier.github.io/hybridclock/reference/class_settings.md)).

## Value

A data.frame with one row per generation: `generation`, the true history
of the sampled hybrids (`hybrid_age_both`, `hybrid_age_rare`,
`backcrosses`, `recent_share` for the share that are F1, F2 or first
backcrosses, `rare_species`, the current `arrival_rate` of that species
(its mean over the three generations before sampling) and `arriving`
(whether that rate is above zero)), and the summaries from
[`sample_summaries()`](https://paulmaier.github.io/hybridclock/reference/sample_summaries.md).

## Details

Genotypes count species-2 alleles: 0 is a pure species-1 animal and 2 a
pure species-2 animal at every diagnostic marker. Shared loci (optional)
are polymorphic in both species at the same frequencies and are used to
measure inbreeding and linkage disequilibrium.

Two hybrid ages and a backcross count are tracked for each hybrid:

- `hybrid_age_both` (hybrid age, both species): every line of the
  pedigree is followed back to its first pure ancestor, of either
  species, and the number of generations is averaged over lines. An F1
  is 1, a first backcross 1.5, an F2 2. Recent backcrossing keeps it
  low.

- `hybrid_age_rare` (hybrid age, rare species): only the lines leading
  back to pure animals of the rarer species in the sample, weighted by
  how much ancestry each contributes: how long ago, on average, the
  rarer species' genes entered.

- `backcrosses`: the number of matings with a pure animal since the
  first cross (averaged over the two parents when both are hybrids).

## Examples

``` r
sim <- simulate_hybrid_zone(generations = 8, breeders = 60, start_mix = 0.2,
                            arrivals = c(0.05, 0), n_markers = 300, seed = 1)
sim[, c("generation", "hybrid_age_both", "hybrid_age_rare", "backcrosses", "mean_S", "mean_H")]
#>   generation hybrid_age_both hybrid_age_rare backcrosses     mean_S    mean_H
#> 1          1        1.000000        1.000000    0.000000 0.22500000 0.4500000
#> 2          2        1.500000        2.000000    1.000000 0.08858333 0.1771667
#> 3          3        2.171875        2.875000    1.125000 0.21000000 0.3160000
#> 4          4        2.930556        4.000000    1.430556 0.19941667 0.3135000
#> 5          5        3.503125        4.988889    1.956250 0.17000000 0.2770000
#> 6          6        4.231250        5.985000    1.984375 0.16558333 0.2768333
#> 7          7        5.079688        6.958701    2.151562 0.16358333 0.2538333
#> 8          8        6.100391        7.948456    1.944531 0.17441667 0.2815000
```
