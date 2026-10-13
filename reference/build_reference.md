# Simulate a reference set of hybrid histories

Runs
[`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)
many times with settings drawn at random from the ranges given, and
stacks every generation of every run. The ranges act as the prior of the
estimator trained on the result, so they should cover what is plausible
for the study system. The breeding population size matters most: past
about five generations of a closed hybrid swarm, drift is the only clock
left, and its speed depends on population size.

## Usage

``` r
build_reference(
  n_sims = 4000,
  generations = 50,
  breeders = c(10, 500),
  start_mix = c(0.01, 0.99),
  arrival_max = 0.1,
  kinds = c(closed = 0.2, one = 0.4, both = 0.4),
  sib_mating = c(0, 0.9),
  pools = 4,
  n_markers = 1000,
  shared_freq = NULL,
  sample_size = 20,
  missing = 0,
  founders = "pure",
  map = NULL,
  map_scale = c(1, 1),
  assort = c(0, 0.8),
  het_dropout = c(0, 0.3),
  freq = NULL,
  microsat = NULL,
  founder_impurity = c(0, 0),
  founder_samples = FALSE,
  class_settings = NULL,
  cores = 1,
  seed = 1
)
```

## Arguments

- n_sims:

  Number of simulated histories.

- generations:

  Generations per history.

- breeders:

  Range of breeding population sizes (drawn on a log scale). An
  effective size (Ne) range is usually the most relevant figure.

- start_mix:

  Range of the founding share of species 2.

- arrival_max:

  Largest share of breeders replaced by arrivals per generation.

- kinds:

  Probabilities of closed, one-species and two-species arrival
  histories.

- sib_mating:

  Range of sibling-mating probability; half the runs use 0.

- pools:

  Largest number of breeding pools; half the runs use 1.

- n_markers, shared_freq, sample_size, missing:

  Passed to
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md).
  Match these to the observed data.

- founders:

  Founding designs to draw from (see
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)):
  `"pure"` for natural contact, or `"F1"`, `"BC1_1"`, `"BC1_2"` for
  experimental or hatchery populations.

- map:

  Linkage: `NULL` for unlinked markers; a marker map (see
  [`random_map()`](https://paulmaier.github.io/hybridclock/reference/random_map.md))
  used by every history; or a list with ranges `chromosomes` and
  `length_cM`, from which each history draws its own random map, so that
  an unknown genome enters the reference as uncertainty. A fixed map
  should hold the observed markers' real positions, in the order of the
  diagnostic markers; then genotype switches along chromosomes are
  summarized too (see
  [`sample_summaries()`](https://paulmaier.github.io/hybridclock/reference/sample_summaries.md)).

- map_scale:

  With a fixed map, a range by which positions are multiplied in each
  history, for an uncertain recombination rate (for example, positions
  in Mb with `map_scale = c(0.5, 3)` for 0.5-3 cM per Mb).

- assort:

  Range of the probability of mating by ancestry (see
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md));
  half the runs use 0.

- het_dropout:

  Range of the probability that a heterozygote is called homozygous;
  half the runs use 0 unless the range excludes 0. Give `c(0, 0)` for
  genotyping methods without allele dropout.

- freq:

  Optional allele frequencies for markers that are not fixed (see
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md));
  observed scores must then come from
  `hybrid_scores(method = "likelihood")`.

- microsat:

  Optional microsatellite markers (see
  [`microsat_markers()`](https://paulmaier.github.io/hybridclock/reference/microsat_markers.md));
  observed scores must then come from
  [`hybrid_scores_microsat()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores_microsat.md),
  and `n_markers` is set to the number of loci.

- founder_impurity:

  Range of the share of diagnostic alleles that founders carry from the
  other species (see
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)),
  drawn for each history.

- founder_samples:

  With F1 or backcross founders, also include samples of the founders
  themselves (generation 0), for experimental families of F1s or first
  backcrosses (see
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)).

- class_settings:

  Settings of the classes behind the class shares in the summaries (see
  [`class_settings()`](https://paulmaier.github.io/hybridclock/reference/class_settings.md)):
  genotyping error, heterozygote dropout, linked markers, and reference
  allele frequencies for markers that are not fixed. They are stored
  with the reference and used again by
  [`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md),
  so that the data are classified the same way.

- cores:

  Number of processes (forked; 1 on Windows).

- seed:

  Base random seed; run `i` uses `seed + i`.

## Value

A data.frame of all generations of all runs, with the run's settings
(`run`, `breeders`, `start_mix`, `arrival1`, `arrival2` and their
`timing1`, `timing2` and `switch1`, `switch2` generations, `kind`,
`sib_mating`, `pools`), the true history and the summaries. Attribute
`settings` records the arguments.

## Details

Each run is one of three kinds of history, with probabilities `kinds`: a
closed swarm (no arrivals after founding), arrivals of one species only
(repeated backcrossing), or arrivals of both species (ongoing
hybridization). Arrivals of each species are, with equal probability,
constant, stop at a random generation, or start at a random generation,
so that the reference includes hybridization that has ended and
hybridization that began late.

## Examples

``` r
# \donttest{
ref <- build_reference(n_sims = 20, generations = 10, n_markers = 200)
dim(ref)
#> [1] 200 127
# }
```
