# ksamplesLTRC 0.1.1

## Bug fixes

* Fixed an incorrect boundary condition in `ltrc.cvm()` and `lt.cvm()`
  that could bias the computed Cramér–von Mises statistic, particularly 
  under unequal group sizes or non-proportional weights.

* Fixed a sign and continuity error in `lt.truncation.estimator()` that
  caused both the estimated truncation distribution and the truncation
  probability estimate (`gamma.hat`) to be incorrect.

## Other changes

* Added a public GitHub repository and issue tracker:
  https://github.com/adrian-lago/ksamplesLTRC

* Substantially expanded the package's unit test suite.

# ksamplesLTRC 0.1.0

* Initial CRAN submission.
