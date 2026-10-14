# Classify recent hybrids individually

Tests each animal against the pure, first- and second-generation hybrid
classes and second backcrosses. Each class is a pair of parents, each
passing on species-2 alleles with a set probability: 0 (pure species 1),
1 (pure species 2), 0.5 (F1 or F2 parent), 0.25 or 0.75 (first backcross
parent). For each animal the likelihood of its genotypes under every
class is compared with the best fit under any pair of parental
ancestries. If no class fits (likelihood-ratio test at `alpha`), the
animal is called a later-generation hybrid.

## Usage

``` r
classify_hybrids(
  scores,
  alpha = 0.01,
  error = 0.005,
  genotypes = NULL,
  markers = NULL,
  effective_markers = NULL,
  het_dropout = 0
)
```

## Arguments

- scores:

  Output of
  [`hybrid_scores()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores.md)
  (`S`, `H`, `n_markers`).

- alpha:

  Significance level of the fit test.

- error:

  Genotyping error rate, or `"estimate"` (with `genotypes`). An estimate
  that reaches the top of its range (0.1) gives a warning: with many
  linked markers, extra variation in heterozygosity can be taken for
  genotyping error, and a fixed rate is then safer.

- genotypes:

  Optional genotypes at the diagnostic markers, counting species-2
  alleles, as from
  [`orient_genotypes()`](https://paulmaier.github.io/hybridclock/reference/orient_genotypes.md),
  one row per row of `scores`.

- markers:

  Output of
  [`diagnostic_markers()`](https://paulmaier.github.io/hybridclock/reference/diagnostic_markers.md)
  matching the columns of `genotypes` (columns `freq1` and `freq2` are
  used).

- effective_markers:

  Optional number of effectively independent markers. Markers on the
  same chromosome are inherited in blocks, so they carry less
  independent information than their number suggests. With this,
  log-likelihoods are scaled by `effective_markers / n_markers` (a
  composite-likelihood adjustment). The best class is unchanged, but
  animals are no longer rejected from every class just because their
  markers come in blocks. See
  [`effective_markers()`](https://paulmaier.github.io/hybridclock/reference/effective_markers.md).

- het_dropout:

  Probability that a heterozygote is called homozygous (allele dropout,
  common in reduced-representation sequencing). It is not estimated,
  because heterozygosity lost to dropout cannot be told apart from
  heterozygosity lost to inbreeding one animal at a time; set it from
  replicate genotyping or known error rates.

## Value

`scores` with columns `class` (best class or `"later"`), `probability`
(of the best class among the tested classes) and `fit_p` (p-value of the
fit test). With `error = "estimate"`, attribute `error` holds the
estimated rate.

## Details

With `genotypes` and `markers`, the likelihood uses each marker's allele
frequencies in the two reference groups, so markers that are nearly but
not perfectly fixed are handled correctly, as in NewHybrids (Anderson
and Thompson 2002) and HIest (Fitzpatrick 2012). A parent with species-2
ancestry `a` passes on the species-2 allele at marker `l` with
probability `a * freq2[l] + (1 - a) * freq1[l]`. The genotyping error
rate can be estimated from the data (`error = "estimate"`). Without
them, markers are treated as fixed and only the genotype counts implied
by `S` and `H` are used.

Several third-generation classes are indistinguishable from earlier ones
at unlinked markers whatever the number of markers: F1 x F2 and F2 x F2
look like F2, and F2 x pure looks like a first backcross. Inbreeding
(related parents) lowers heterozygosity and moves animals off every
class, so they are called later-generation hybrids. Use
[`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md)
for populations of advanced hybrids. Markers must be diagnostic in the
populations that actually hybridized, not only in the reference samples:
a marker at which one source population carries the "wrong" allele makes
its hybrids look like no class.

## References

Anderson EC, Thompson EA (2002) A model-based method for identifying
species hybrids using multilocus genetic data. Genetics 160:1217-1229.

## Examples

``` r
ex <- hybridclock_example()
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
sc <- hybrid_scores(ex$genotypes, mk)
table(ex$groups, classify_hybrids(sc)$class)
#>           
#>            BC2_species1 F2 later pure1 pure2
#>   site_A              0 17    13     0     0
#>   site_B              3  0    27     0     0
#>   species1            0  0     0    20     0
#>   species2            0  0     0     0    20
cl <- classify_hybrids(sc, error = "estimate", genotypes = orient_genotypes(ex$genotypes, mk),
                       markers = mk)
attr(cl, "error")
#> [1] 5e-04
```
