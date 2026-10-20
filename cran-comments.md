## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* local: macOS 26.6 (aarch64), R 4.6.1
* win-builder: R-release, R-devel
* macOS builder (mac.r-project.org): R-release

## Notes

* The package depends only on base R (stats, utils, graphics, grDevices).
  The packages in Suggests are used only in the vignettes and tests, behind
  `requireNamespace()` guards.
* The two DOIs in the Description field point to the methodological
  references (Chan 2022, Quantitative Economics; Rubio-Ramirez, Waggoner and
  Zha 2010, Review of Economic Studies).
