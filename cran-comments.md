## Comments

## 2017-8-2

Changed the numbers around to execute the example and in testing! Thank you for the suggestion.

(Initial errors were caused when using 4 cores.  It does work with 2 cores!)

- Fulya


## 2017-8-2

Thanks, we see the examples do not conatin a single lmmpar call outside of \sontrun{} and the tests are not executed on CRAN machines, hence your package does not get tested appropriately.

Why not add an lmmpar example running on 2 cores?

Best,
Uwe Ligges


## Test environments
* local OS X install, R 3.4.0
* ubuntu 14.04.5 (on travis-ci), R 3.4.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
