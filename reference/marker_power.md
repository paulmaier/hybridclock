# Chance that a backcrossed hybrid looks pure

After `k` generations of backcrossing to one species, each fixed
diagnostic marker is heterozygous with probability `0.5^k`, so the
chance that none of `n` markers is heterozygous (and the hybrid passes
as pure) is `(1 - 0.5^k)^n` (Boecklen and Howard 1997).

## Usage

``` r
marker_power(n_markers, generations = 1:10)

markers_needed(generations = 1:10, error = 0.05)
```

## Arguments

- n_markers:

  Number of diagnostic markers (vector allowed).

- generations:

  Generations of backcrossing (1 is a first backcross).

- error:

  Acceptable chance of missing a hybrid.

## Value

A data.frame of `generations`, `n_markers` and `p_missed`.

For `markers_needed()`, the smallest number of markers that keeps the
chance of missing a hybrid of each backcross generation below `error`.

## References

Boecklen WJ, Howard DJ (1997) Genetic analysis of hybrid zones: numbers
of markers and power of resolution. Ecology 78:2611-2616.

## Examples

``` r
marker_power(c(100, 1000), generations = 1:8)
#>    generations n_markers      p_missed
#> 1            1       100  7.888609e-31
#> 2            2       100  3.207202e-13
#> 3            3       100  1.587835e-06
#> 4            4       100  1.574446e-03
#> 5            5       100  4.179954e-02
#> 6            6       100  2.070416e-01
#> 7            7       100  4.564310e-01
#> 8            8       100  6.761165e-01
#> 9            1      1000 9.332636e-302
#> 10           2      1000 1.151499e-125
#> 11           3      1000  1.018716e-58
#> 12           4      1000  9.360012e-29
#> 13           5      1000  1.628229e-14
#> 14           6      1000  1.447356e-07
#> 15           7      1000  3.924202e-04
#> 16           8      1000  1.996251e-02
markers_needed(1:10)
#>    generations n_markers
#> 1            1         5
#> 2            2        11
#> 3            3        23
#> 4            4        47
#> 5            5        95
#> 6            6       191
#> 7            7       382
#> 8            8       766
#> 9            9      1533
#> 10          10      3067
```
