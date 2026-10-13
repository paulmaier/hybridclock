# Find loci shared by both species

Keeps SNPs that are common in both reference groups and differ little
between them. Inbreeding and linkage disequilibrium measured at these
loci help separate inbreeding from hybrid age (see
[`sample_summaries()`](https://paulmaier.github.io/hybridclock/reference/sample_summaries.md)).

## Usage

``` r
shared_loci(
  geno,
  groups,
  ref1,
  ref2,
  min_maf = 0.05,
  max_diff = 0.1,
  locus = NULL
)
```

## Arguments

- geno:

  Genotype matrix, animals in rows and SNPs in columns, counting copies
  of one allele (0, 1, 2; `NA` for missing).

- groups:

  Group label of each animal.

- ref1, ref2:

  Labels of the two reference groups. `ref2` defines species 2, whose
  alleles the hybrid index counts.

- min_maf:

  Smallest minor allele frequency required in each group.

- max_diff:

  Largest allele frequency difference allowed between groups.

- locus:

  Optional locus name of each SNP.

## Value

A data.frame with `column`, `name` and the pooled allele frequency
`freq`.

## Examples

``` r
ex <- hybridclock_example()
nrow(shared_loci(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus))
#> [1] 106
```
