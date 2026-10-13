# Effective number of independent markers

Markers on one chromosome travel together through a meiosis, broken only
by crossovers. An early-generation hybrid's genotype is shaped by the
recombinant gamete(s) of its hybrid parent(s): one for a backcross, two
for an F2. The number of independently inherited pieces is about
`meioses * (chromosomes + map length in Morgans)`. The default of one
meiosis suits backcrosses, the classes most affected; in the bear
pedigree truth set it let every known backcross be classified correctly
while animals with old introgression were still called later-generation.

## Usage

``` r
effective_markers(chromosomes, length_cM = 100 * chromosomes, meioses = 1)
```

## Arguments

- chromosomes:

  Number of chromosomes (haploid).

- length_cM:

  Total map length in centimorgans. If unknown, about 50-150 cM per
  chromosome is typical.

- meioses:

  Number of recombinant gametes shaping the genotype (1 for a backcross,
  2 for an F2).

## Value

A number to pass to
[`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md)
as `effective_markers`.

## Examples

``` r
effective_markers(chromosomes = 37, length_cM = 2500)
#> [1] 62
```
