# Random marker map

Places markers at random positions on chromosomes, for simulating
linkage when no genetic map is available. Use the species' chromosome
number and, if known, its map length; otherwise give
[`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md)
a range so that the uncertainty is part of the reference.

## Usage

``` r
random_map(n_markers, chromosomes, length_cM)
```

## Arguments

- n_markers:

  Number of markers.

- chromosomes:

  Number of chromosomes.

- length_cM:

  Total map length in centimorgans, shared equally among chromosomes.

## Value

A data.frame with `chromosome` and `position` (cM), sorted.

## Examples

``` r
head(random_map(500, chromosomes = 11, length_cM = 1500))
#>    chromosome  position
#> 15          1  5.816880
#> 9           1  6.923341
#> 5           1  7.035726
#> 14          1 11.781188
#> 24          1 16.564295
#> 29          1 17.531547
```
