# Argument validation: every exported function must fail early with a message
# that names the offending argument.

sim <- simulate_var_example(T = 80, seed = 7)
Y   <- sim$Y

test_that("post_sim validates its inputs", {
    expect_error(post_sim(as.data.frame(Y), p = 1L, R = 5L), "numeric matrix")
    expect_error(post_sim(Y[, 1, drop = FALSE], p = 1L, R = 5L), "at least 2 columns")
    Yna <- Y; Yna[3, 2] <- NA
    expect_error(post_sim(Yna, p = 1L, R = 5L), "missing")
    expect_error(post_sim(Y, p = 0, R = 5L), "`p`")
    expect_error(post_sim(Y, p = 1.5, R = 5L), "`p`")
    expect_error(post_sim(Y, p = nrow(Y), R = 5L), "`p`")
    expect_error(post_sim(Y, p = 1L, R = 0), "`R`")
    expect_error(post_sim(Y, p = 1L, R = 5L, kappa3 = -1), "`kappa3`")
    expect_error(post_sim(Y, p = 1L, R = 5L, var_names = c("a", "b")), "`var_names`")
})

test_that("build_yX, build_s2 and optim_hyper validate their inputs", {
    expect_error(build_yX(Y, p = -1), "`p`")
    expect_error(build_yX(c(1, 2, 3), p = 1), "numeric matrix")
    expect_error(build_s2(list("a")), "numeric vectors")
    dat <- build_yX(Y, 1L); s2 <- build_s2(dat$y_list)$s2
    expect_error(optim_hyper(dat$y_list, dat$X_list[1:2], 1L, s2), "same length")
    expect_error(optim_hyper(dat$y_list, dat$X_list, 1L, s2[1:2]), "`s2`")
    expect_error(optim_hyper(dat$y_list, dat$X_list, 1L, -s2), "`s2`")
    expect_error(optim_hyper(dat$y_list, dat$X_list, 1L, s2, init = c(0.1, -0.1)), "`init`")
})

test_that("sign_restrict validates mcmc, S and the integer arguments", {
    mcmc <- post_sim(Y, p = 1L, R = 5L)
    S <- rbind(c(1, -1, 0), c(1, 1, 0), c(-1, -1, 1))
    expect_error(sign_restrict(list(), S), "post_sim")
    expect_error(sign_restrict(mcmc, S[, 1:2]), "columns")
    S2 <- S; S2[1, 1] <- 2
    expect_error(sign_restrict(mcmc, S2), "\\{-1, 0, 1\\}")
    S0 <- S; S0[3, ] <- 0
    expect_error(sign_restrict(mcmc, S0), "all zero")
    expect_error(sign_restrict(mcmc, S, n_id = 4), "`n_id`")
    expect_error(sign_restrict(mcmc, S, max_tries = 0), "`max_tries`")
    expect_error(sign_restrict(mcmc, S, verbose = NA), "`verbose`")
})

test_that("collect_irfs, compute_irf and plot_irfs validate their inputs", {
    mcmc <- post_sim(Y, p = 1L, R = 5L)
    S <- rbind(c(1, -1, 0), c(1, 1, 0), c(-1, -1, 1))
    acc <- sign_restrict(mcmc, S, max_tries = 500L, verbose = FALSE)
    expect_error(collect_irfs(list(), n = 3L, p = 1L, H = 4L), "empty")
    expect_error(collect_irfs(acc, n = 2L, p = 1L, H = 4L), "`n` = 2")
    expect_error(collect_irfs(acc, n = 3L, p = 2L, H = 4L), "`p` = 2")
    expect_error(collect_irfs(acc, n = 3L, p = 1L, H = -1), "`H`")
    expect_error(compute_irf(acc[[1]]$impact, acc[[1]]$B_rf[, 1:2], 3L, 1L, 4L), "`B_rf`")
    irfs <- collect_irfs(acc, n = 3L, p = 1L, H = 4L)
    expect_error(plot_irfs(irfs[, , , 1]), "4-dimensional")
    expect_error(plot_irfs(irfs, probs = c(0.9, 0.1)), "`probs`")
    expect_error(plot_irfs(irfs, var_names = "x"), "`var_names`")
})

test_that("recover_rf and simulate_var_example validate their inputs", {
    expect_error(recover_rf(list(1, 2), c(1)), "same length")
    expect_error(recover_rf(list(1, 2), c(1, -1)), "positive variances")
    expect_error(simulate_var_example(T = 3), "`T`")
})

test_that("check_S accepts a matrix or a 3-d array and rejects the rest", {
    S <- rbind(c(1, -1, 0), c(1, 1, 0))
    a <- acpbvar:::check_S(S, 3, horizons = 2L)
    expect_equal(dim(a), c(2, 3, 3))
    expect_true(all(a[, , 2] == S))
    A <- array(0, c(2, 3, 3)); A[1, 1, ] <- 1; A[2, 2, 2] <- -1
    expect_equal(dim(acpbvar:::check_S(A, 3, horizons = 2L)), c(2, 3, 3))
    expect_error(acpbvar:::check_S(A, 3, horizons = 1L), "slices")
    A0 <- A; A0[2, , ] <- 0
    expect_error(acpbvar:::check_S(A0, 3, horizons = 2L), "all zero")
    expect_error(acpbvar:::check_S(1:3, 3), "3-d array")
    sim  <- simulate_var_example(T = 80, seed = 7)
    mcmc <- post_sim(sim$Y, p = 1L, R = 5L)
    expect_error(sign_restrict(mcmc, S, horizons = -1), "`horizons`")
    expect_error(sign_restrict(mcmc, S, horizons = 1.5), "`horizons`")
})
