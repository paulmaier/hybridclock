# Changelog

## hybridclock 0.1.1

- `classify_hybrids(error = "estimate")` warns when the estimated
  genotyping error rate reaches the top of its range (0.1). With many
  linked markers, as in F2 families, extra variation in heterozygosity
  can be taken for genotyping error; a fixed rate (the default, 0.005)
  is then safer.
- The getting-started vignette links to the full tutorial.
- The package description quotes ‘Cubist’ and cites Fitzpatrick (2012)
  in the form CRAN asks for.

## hybridclock 0.1.0

- First release: diagnostic markers and hybrid scores, individual
  classification of recent hybrids, forward simulation of hybrid
  populations, and estimates of hybrid age (both species and rare
  species), backcrosses and arrival rate with split-conformal intervals.

- New hybridization is estimated in two parts: the probability that the
  rarer species is still arriving (`arriving`), and the rate if it is
  (`arrival_rate`). Reference histories include arrivals that stop or
  start partway, and
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)
  accepts arrival schedules.

- The two ages are named `hybrid_age_both` and `hybrid_age_rare`.

- [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)
  and
  [`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md)
  gain `founders` (pure, F1 or first-backcross founders, for
  experimental populations) and `map` (linkage, from a marker map or a
  range of chromosome numbers and map lengths; see
  [`random_map()`](https://paulmaier.github.io/hybridclock/reference/random_map.md)).

- [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md)
  can use the reference allele frequencies of each marker and estimate
  the genotyping error rate;
  [`classification_power()`](https://paulmaier.github.io/hybridclock/reference/classification_power.md)
  gives the expected accuracy of individual classes.

- Class shares from
  [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md)
  are part of the population summaries, and the arrival probability is
  recalibrated by estimated age; young swarms are flagged.

- [`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md)
  gives one answer per population from two tiers: individual classes
  when most hybrids are F1s and backcrosses, otherwise the population
  clock (column `method`).

- [`effective_markers()`](https://paulmaier.github.io/hybridclock/reference/effective_markers.md)
  and the `effective_markers` argument of
  [`classify_hybrids()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids.md)
  and
  [`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md)
  correct the classifier for linked markers.

- Heterozygote dropout (`het_dropout`) and mating by ancestry (`assort`)
  in the simulator and reference; sibling mating now ranges up to 0.9.

- With a fixed marker map, genotype switches along chromosomes
  (junctions) are summarized; `map_scale` allows for an uncertain
  recombination rate.

- Markers that are not fixed: `hybrid_scores(method = "likelihood")`
  estimates S and H from reference allele frequencies (Fitzpatrick
  2012),
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)
  and
  [`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md)
  gain `freq`, and
  [`marker_information()`](https://paulmaier.github.io/hybridclock/reference/marker_information.md)
  summarizes how informative a marker set is.

- Microsatellites:
  [`microsat_markers()`](https://paulmaier.github.io/hybridclock/reference/microsat_markers.md),
  [`hybrid_scores_microsat()`](https://paulmaier.github.io/hybridclock/reference/hybrid_scores_microsat.md),
  [`classify_hybrids_microsat()`](https://paulmaier.github.io/hybridclock/reference/classify_hybrids_microsat.md),
  and `microsat` in the simulator and reference.

- The share of early-class hybrids needed for the individual-classes
  answer is learned from the reference without a fixed cap; with weak
  markers the classes are never used.

- `founder_impurity` in
  [`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)
  and
  [`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md):
  founders that carry some of the other species’ alleles at the
  diagnostic markers, such as cross parents with old introgression.

- Intervals of class-based answers include the error of class-based ages
  seen in the reference, not only the resampling of animals.

- The individual-classes answer is chosen by the classes’ large errors
  (90th percentile) on held-out reference histories, needs at least 80%
  of hybrids in early classes, and its intervals are at least +/-10%
  wide.

- [`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md)
  reports `information` (how much narrower the hybrid-age interval is
  than the reference’s age range) and flags answers that mostly reflect
  the reference’s age range; it gains `het_dropout` for the individual
  classes.

- `founder_samples`
  ([`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md))
  and `sample_founders`
  ([`simulate_hybrid_zone()`](https://paulmaier.github.io/hybridclock/reference/simulate_hybrid_zone.md)):
  references for lab families of F1s or first backcrosses.

- [`check_markers()`](https://paulmaier.github.io/hybridclock/reference/check_markers.md)
  flags markers that are not diagnostic in the animals that actually
  bred, using likely F1s.

- Samples made up mostly of F1s may be dated by their classes even with
  weak markers, if the reference shows this is accurate (F1 rule).

- Intervals are calibrated separately for samples made up mostly of
  F2-like animals and for all others.

- Tutorial: keeping reference ranges generous, breeding design, map
  length uncertainty, checking markers, and the number of markers needed
  to date F2 families.

- Intervals are now normalized split-conformal by default
  (`train_clock(intervals = "normalized")`): centred on the estimate,
  scaled by a forest’s prediction of each sample’s typical error, and
  calibrated on held-out runs. Coverage is even across ages and sample
  types (F2 families included). `"binned"` (earlier behavior, now also
  by sample type) and `"cqr"` (conformalized quantile regression) remain
  available.

- [`class_settings()`](https://paulmaier.github.io/hybridclock/reference/class_settings.md):
  settings of the classes behind the class shares (error, heterozygote
  dropout, linked markers, markers that are not fixed), given to
  [`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md),
  stored with the reference and reused by
  [`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md),
  so that the reference and the data are classified the same way. The
  `markers`, `effective_markers`, `het_dropout` and `classes` arguments
  of
  [`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md)
  no longer change its answer.

- [`estimate_history()`](https://paulmaier.github.io/hybridclock/reference/estimate_history.md)
  sets aside animals without scores, with a warning.

- [`build_reference()`](https://paulmaier.github.io/hybridclock/reference/build_reference.md)
  draws heterozygote dropout in every run when the given range excludes
  zero.

- Example datasets: `wolves`, `chickens`, `bettongs`, `possums`, `carp`.
