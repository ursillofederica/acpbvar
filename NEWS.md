# acpbvar 0.1.0

Initial release (package formerly named `samplerChan`).

* Bayesian structural VAR in the A-form of Chan (2022), with the asymmetric
  conjugate (Minnesota-type) prior; own-lag and cross-lag tightness chosen by
  Empirical Bayes on the closed-form marginal likelihood.
* Exact equation-by-equation posterior simulation (independent
  Normal-Inverse-Gamma draws), analytical recovery of the reduced form.
* Sign-restriction identification on the impact matrix following
  Rubio-Ramirez, Waggoner and Zha (2010); impulse responses with posterior
  bands.

Fixes relative to the development version of May 2026:

* `check_sign()` now checks every variable under partial identification
  (it previously looped over the number of shocks), and `sign_restrict()`
  defaults `n_id` to `nrow(S)`.
* `sign_restrict()` records the posterior draws for which no rotation was
  found within `max_tries` (attribute `failed_draws`) and warns when any.
* `optim_hyper()` no longer imposes a discontinuous `kappa < 1` penalty;
  `post_sim()` stores the optimiser convergence code in `spec$convergence`
  and warns on non-convergence.
* Added testthat regression tests and a DGP recovery test.
* All exported functions validate their arguments and stop with an
  informative message (`R/validate.R`).
