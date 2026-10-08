# DGP recovery test (CLAUDE.md 10a): on the package's own simulated VAR(1),
# the posterior must recover the reduced-form persistence, the reduced-form
# innovation covariance and the sign pattern of the structural impact matrix.
# Thresholds are deliberately loose: they guard against a broken pipeline
# (wrong sign convention, wrong lag alignment, wrong recovery of A), not
# against ordinary sampling noise.

skip_on_cran()

sim  <- simulate_var_example(T = 400, seed = 20261008)
mcmc <- post_sim(sim$Y, p = 1L, R = 1500L)

rf <- lapply(seq_len(mcmc$spec$R), function(r)
    recover_rf(mcmc$samples$theta[[r]], mcmc$samples$sigma2[r, ]))
B_mean     <- Reduce(`+`, lapply(rf, `[[`, "B_rf"))     / length(rf)
Sigma_mean <- Reduce(`+`, lapply(rf, `[[`, "Sigma_rf")) / length(rf)

test_that("reduced-form lag matrix is recovered", {
    B1 <- B_mean[, -1]                      # drop intercept column
    expect_lt(max(abs(B1 - sim$B_true)), 0.15)
    expect_lt(max(abs(B_mean[, 1] - sim$intercept_true)), 0.15)
})

test_that("reduced-form innovation covariance is recovered", {
    Sigma_true <- sim$impact_true %*% t(sim$impact_true)
    expect_lt(max(abs(Sigma_mean - Sigma_true)), 0.15)
})

test_that("sign-identified impact matrix recovers the true signs", {
    S <- rbind(supply = c(1, -1, 0), demand = c(1, 1, 0), monetary = c(-1, -1, 1))
    set.seed(1)
    acc <- sign_restrict(mcmc, S, max_tries = 2000L, verbose = FALSE)
    expect_gt(length(acc), 0.9 * mcmc$spec$R)   # the identified set is rarely empty here
    impact_mean <- Reduce(`+`, lapply(acc, `[[`, "impact")) / length(acc)
    restricted  <- S != 0
    expect_true(all(sign(t(impact_mean))[restricted] == S[restricted]))
})
