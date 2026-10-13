# Summarize a population sample of hybrid scores

Turns the hybrid index (S) and interspecies heterozygosity (H) of a
sample of animals, and optionally their genotypes, into the summary
statistics that
[`train_clock()`](https://paulmaier.github.io/hybridclock/reference/train_clock.md)
learns from. The same function is used for simulated and observed
samples, so the two are always comparable.

## Usage

``` r
sample_summaries(
  S,
  H,
  diag = NULL,
  shared = NULL,
  n_markers = NULL,
  positions = NULL,
  class_settings = NULL
)
```

## Arguments

- S, H:

  Hybrid index and interspecies heterozygosity of each animal, as from
  [`hybrid_scores()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores.md).

- diag:

  Optional genotype matrix (animals x diagnostic markers) counting
  species-2 alleles, as from
  [`orient_genotypes()`](https://paulmaier.github.io/hybridclock/reference/orient_genotypes.md).
  Missing values allowed.

- shared:

  Optional genotype matrix (animals x shared loci) of allele counts.
  Missing values allowed.

- n_markers:

  Optional number of markers genotyped in each animal, as from
  [`hybrid_scores()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores.md).

- positions:

  Optional marker positions matching the columns of `diag`: a data.frame
  with `chromosome` and `position`. With them, the rate at which each
  hybrid's genotype switches between neighbouring markers along a
  chromosome is summarized (`switch_rate`, `switch_rate_sd`). Switches
  are ancestry junctions plus genotyping noise; they accumulate with
  every generation of recombination (Janzen et al. 2018), so they date
  early generations even when animals look alike in S and H.

- class_settings:

  Settings of the classes behind the class shares (see
  [`class_settings()`](https://paulmaier.github.io/hybridclock/reference/class_settings.md));
  the same settings must be used for the reference and the data.

## Value

A named numeric vector.

## Details

The sample is first oriented so that the rarer species among its hybrids
counts as species 2 (S below 0.5). The summaries then describe:

- the spread of animals over the S-H triangle (means, percentiles, and
  the share of animals in each cell of a 10 x 6 grid);

- heterozygosity relative to its maximum for each animal's ancestry
  (`2 * min(S, 1 - S)`) and to the minimum possible for unrelated
  parents (`2 * S * (1 - S)`). Values below the latter point to
  inbreeding;

- with `diag`, how unevenly the rarer species' alleles are spread across
  markers (drift);

- with `shared`, inbreeding (`shared_fis`) and excess linkage
  disequilibrium (`shared_ld`) at loci that do not differ between
  species;

- with `n_markers`, the share of animals in each class of
  [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md)
  (pure, F1, F2, first and second backcrosses toward the common or the
  rare species, and later generations).

## Examples

``` r
set.seed(1)
S <- c(0, 0, 0.12, 0.25, 0.3, 0.06, 0.5)
H <- c(0, 0, 0.22, 0.48, 0.41, 0.12, 1)
round(sample_summaries(S, H)[1:12], 3)
#>         n_hybrids share_pure_common   share_pure_rare            mean_S 
#>             5.000             0.286             0.000             0.176 
#>              sd_S             q10_S             q50_S             q90_S 
#>             0.185             0.000             0.120             0.380 
#>            mean_H              sd_H             q10_H             q50_H 
#>             0.319             0.354             0.000             0.220 
```
