# Example data set

Simulates a small data set: 20 animals from each of two pure reference
groups and two hybrid sites. Site A is a young hybrid swarm (4
generations, no arrivals); site B has been backcrossing into species 1
for 15 generations. Animals carry 300 diagnostic SNPs, 40 SNPs that are
not quite fixed, and 150 shared SNPs, on 490 separate loci. The data are
simulated afresh with a fixed seed, so every call returns the same data.

## Usage

``` r
hybridclock_example()
```

## Value

A list with `genotypes` (animals x SNPs, counts of the reference
allele), `groups` (`"species1"`, `"species2"`, `"site_A"`, `"site_B"`),
`locus` (locus of each SNP) and `truth` (the simulated histories of the
two sites at the time of sampling).

## Examples

``` r
ex <- hybridclock_example()
table(ex$groups)
#> 
#>   site_A   site_B species1 species2 
#>       30       30       20       20 
ex$truth
#>      site generation hybrid_age_both hybrid_age_rare backcrosses arrival_rate
#> 4  site_A          4        3.641667        3.754127   0.5416667            0
#> 15 site_B         15        9.592318       14.787010   2.7037638            0
```
