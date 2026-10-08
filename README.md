# acpbvar

Bayesian estimation of asymmetric structural VARs in A-form, with
Empirical-Bayes hyperparameter selection and sign-restriction identification.

## Overview

`acpbvar` implements the asymmetric A-form representation of structural
Vector Autoregressions following Chan (2022). The package combines

- an asymmetric Minnesota prior, with own-lag, cross-lag, and intercept
  tightness selected by Empirical Bayes through the closed-form log marginal
  likelihood;
- exact equation-by-equation posterior simulation in the structural form,
  with independent Normal-Inverse-Gamma draws;
- analytical recovery of the reduced-form quantities (impact matrix, lag
  coefficients, innovation covariance) from each posterior draw;
- structural identification via sign restrictions, following the rotation
  algorithm of Rubio-Ramirez, Waggoner, and Zha (2010);
- posterior impulse responses with credible bands.

The package has no external dependencies beyond base R.

## Installation

```r
# install.packages("devtools")
devtools::install_github("ursillofederica/acpbvar", build_vignettes = TRUE)
```

## Quick start

```r
library(acpbvar)

# Simulate a small structural VAR with known sign pattern
sim <- simulate_var_example(T = 200, seed = 42)
Y   <- sim$Y

# Posterior simulation under the asymmetric A-form prior
mcmc <- post_sim(Y, p = 1L, R = 1000L)

# Sign restrictions: supply (y+, pi-), demand (y+, pi+), monetary (y-, pi-, r+)
S <- rbind(
    supply   = c( 1, -1,  0),
    demand   = c( 1,  1,  0),
    monetary = c(-1, -1,  1)
)

accepted <- sign_restrict(mcmc, S, n_id = 3L, max_tries = 1000L)

# Impulse responses
irfs <- collect_irfs(accepted, n = 3L, p = 1L, H = 12L)
plot_irfs(irfs,
          var_names   = sim$var_names,
          shock_names = sim$shock_names)
```

The vignette `acpbvar_workflow` documents the full workflow, including
diagnostics and recovery of the true structural impact from the simulated
data-generating process.

## References

- Chan, J. C. C. (2022). Asymmetric conjugate priors for large Bayesian VARs.
  *Quantitative Economics*, 13(3), 1145-1169.
- Rubio-Ramirez, J. F., Waggoner, D. F., & Zha, T. (2010). Structural Vector
  Autoregressions: Theory of Identification and Algorithms for Inference.
  *Review of Economic Studies*, 77(2), 665-696.

## License

MIT &copy; 2026 Federica Ursillo.
