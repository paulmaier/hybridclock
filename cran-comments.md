## Resubmission

This version replaces my pending submission of hybridclock 0.1.0 (submitted
2026-10-13), which has not been reviewed yet. In this version I have:

* Put the package name 'Cubist' in single quotes in the Description, and
  written the reference as Fitzpatrick (2012) <doi:...>.
* Added a warning in `classify_hybrids()` when the estimated genotyping error
  rate reaches the top of its range.
* Linked the full online tutorial from the vignette.

## Test environments

* local macOS (aarch64), R 4.5.3
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

Examples that build and train a reference take more than 5 s and are in
`\donttest{}`. Examples, tests and the vignette use a single core.
