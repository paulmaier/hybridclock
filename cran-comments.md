## Test environments

* local macOS (aarch64), R 4.5.3
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

Examples that build and train a reference take more than 5 s and are in
`\donttest{}`. Examples, tests and the vignette use a single core.
