# Chickens: an advanced intercross after eight generations

Generation F8 of an advanced intercross line between Red Junglefowl and
White Leghorn chickens: one Red Junglefowl male crossed with three White
Leghorn females, then kept as a closed line of about 100 birds per
generation with no backcrossing. Every bird is eight generations from
the pure parents.

## Usage

``` r
chickens
```

## Format

A list with:

- genotypes:

  Integer matrix, 572 birds x 514 autosomal SNPs, counting White Leghorn
  (species 2) alleles (0, 1, 2; `NA` missing).

- groups:

  `"F8"` for every bird. No pure birds were genotyped.

- map:

  Data frame with `snp`, `chromosome`, `position` (cM per meiosis,
  approximate: the authors' F8 map divided by 3.09, their genome-wide
  map expansion over eight generations) and `position_f8` (cM on the
  authors' F8 map).

- truth:

  Data frame with `group`, `n`, `generation`, `hybrid_age` (8),
  `backcrosses` (0) and `breeders`.

- description:

  A short description.

## Source

Johnsson M (2017) Open field phenotypes and genotypes (R/qtl csv
format). figshare,
[doi:10.6084/m9.figshare.1265060](https://doi.org/10.6084/m9.figshare.1265060)
. CC0 1.0.

## Truth

All birds are F8, so the hybrid age is 8, with no backcrosses and no
arrivals. The line was kept at about 100 birds per generation from F3 to
F7 and produced the F8 from 122 F7 parents; the mating scheme is not
reported. The loss of heterozygosity (mean H about 0.47 instead of 0.5)
suggests about 50 effective breeders. The founders were outbred, so not
every marker was fixed between them, which also lowers H. Which allele
is White Leghorn is not stated by the authors; it was inferred from the
Z chromosome (excluded here), where White Leghorn alleles make up one
third of copies because of the direction of the founding cross.

## Licence

The genotypes are from figshare and were released under CC0 1.0 (public
domain dedication). Markers with partial calls, on the Z chromosome or
on an unplaced scaffold were dropped, and genotypes were recoded.

## References

Johnsson M, Williams MJ, Jensen P, Wright D (2016) Genetical genomics of
behavior: a novel chicken genomic model for anxiety behavior. Genetics
202:327-340.
[doi:10.1534/genetics.115.179010](https://doi.org/10.1534/genetics.115.179010)

Johnsson M, Gustafson I, Rubin C-J, et al. (2012) A sexual ornament in
chickens is affected by pleiotropic alleles at HAO1 and BMP2, selected
during domestication. PLoS Genetics 8:e1002914.
[doi:10.1371/journal.pgen.1002914](https://doi.org/10.1371/journal.pgen.1002914)

## Examples

``` r
# every SNP already counts White Leghorn alleles
g <- chickens$genotypes
mk <- data.frame(column = seq_len(ncol(g)), name = colnames(g), flip = FALSE,
                 freq1 = 0, freq2 = 1)
sc <- hybrid_scores(g, mk)
colMeans(sc[, c("S", "H")])
#>         S         H 
#> 0.4774355 0.4719910 
```
