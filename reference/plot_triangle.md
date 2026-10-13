# Plot animals on the hybrid triangle

Plots interspecies heterozygosity (H) against hybrid index (S). Possible
values fill a triangle with pure animals at the bottom corners and F1s
at the top. The expected positions of first- and second-generation
classes and second backcrosses are marked. Animals below the dashed
curve `H = 2S(1 - S)` are less heterozygous than any offspring of
unrelated parents can be on average, a sign of inbreeding.

## Usage

``` r
plot_triangle(scores, group = NULL, col = NULL, classes = TRUE, ...)
```

## Arguments

- scores:

  Output of
  [`hybrid_scores()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores.md).

- group:

  Optional grouping for colors.

- col:

  Colors, one per group.

- classes:

  Mark the expected class positions.

- ...:

  Passed to
  [`graphics::points()`](https://rdrr.io/r/graphics/points.html).

## Value

`NULL`, invisibly.

## Examples

``` r
ex <- hybridclock_example()
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
plot_triangle(hybrid_scores(ex$genotypes, mk), ex$groups)
```
