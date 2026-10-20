## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* "Possibly misspelled words in DESCRIPTION": Autoregressions, VARs and
  Rubio (as in Rubio-Ramirez, Waggoner and Zha 2010) are correct.

## Test environments

* local: macOS 26.6 (aarch64), R 4.6.1 — 1 NOTE (new submission)
* win-builder: R-release 4.6.1 and R-devel (x86_64-w64-mingw32) — 1 NOTE
  (new submission)
* macOS builder (mac.r-project.org): R 4.6.1 Patched, arm64 — OK

## Notes

* The package depends only on base R (stats, utils, graphics, grDevices).
  The packages in Suggests are used only in the vignettes and tests, behind
  `requireNamespace()` guards.
* The two DOIs in the Description field point to the methodological
  references (Chan 2022, Quantitative Economics; Rubio-Ramirez, Waggoner and
  Zha 2010, Review of Economic Studies).
