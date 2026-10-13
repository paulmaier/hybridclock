# Find diagnostic markers between two reference groups

Keeps SNPs that are nearly fixed for different alleles in two reference
groups of (presumed) pure animals: at most `threshold` frequency in one
group and at least `1 - threshold` in the other. With `locus`, only the
first such SNP per locus is kept, so markers are effectively unlinked.

## Usage

``` r
diagnostic_markers(geno, groups, ref1, ref2, threshold = 0.01, locus = NULL)
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

- threshold:

  Largest frequency allowed for the minor allele in each reference group
  (0.01 keeps markers within 1% of fixation).

- locus:

  Optional locus name of each SNP.

## Value

A data.frame with one row per marker: `column` (index in `geno`),
`name`, `flip` (`TRUE` when the counted allele belongs to species 1),
and the counted allele's frequency in each reference group (`freq1`,
`freq2`).

## Examples

``` r
ex <- hybridclock_example()
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
nrow(mk)
#> [1] 314
```
