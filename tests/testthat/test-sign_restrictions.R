# Regression tests for the orientation of S (n_id x n_var): rows index
# identified shocks, columns index variables. Before the fix, check_sign()
# looped over nrow(S) variables and sign_restrict() defaulted n_id to ncol(S),
# so partial identification (n_id < n_var) silently dropped restrictions on
# the last variables, or failed with a subscript error.

test_that("check_sign reads restrictions on every variable under partial identification", {
    S_part <- rbind(supply = c(1, -1, -1),
                    demand = c(1,  1, -1))          # 2 shocks, 3 variables
    impact_ok  <- matrix(c(1, -1, -1,   1, 1, -1,   0, 0, 0), 3, 3)
    impact_bad <- matrix(c(1, -1,  5,   1, 1,  5,   0, 0, 0), 3, 3)  # violates variable 3
    expect_true(acpbvar:::check_sign(impact_ok,  S_part, n_id = 2))
    expect_false(acpbvar:::check_sign(impact_bad, S_part, n_id = 2))
})

test_that("sign_restrict defaults n_id to the number of rows of S and runs under partial identification", {
    sim  <- simulate_var_example(T = 120, seed = 1)
    mcmc <- post_sim(sim$Y, p = 1L, R = 20L, seed = 1)
    S_part <- rbind(supply = c(1, -1, 0),
                    demand = c(1,  1, 0))
    acc <- sign_restrict(mcmc, S_part, max_tries = 500L, verbose = FALSE)
    expect_gt(length(acc), 0)
    for (a in acc) {
        expect_true(all(a$impact[1, 1:2] >= -1e-12))   # output responds positively to both
        expect_true(a$impact[2, 1] <= 1e-12)            # inflation falls after supply
        expect_true(a$impact[2, 2] >= -1e-12)           # inflation rises after demand
    }
})

test_that("full identification (square S) behaves as before", {
    sim  <- simulate_var_example(T = 120, seed = 2)
    mcmc <- post_sim(sim$Y, p = 1L, R = 20L, seed = 2)
    S <- rbind(supply = c(1, -1, 0), demand = c(1, 1, 0), monetary = c(-1, -1, 1))
    acc <- sign_restrict(mcmc, S, n_id = 3L, max_tries = 500L, verbose = FALSE)
    expect_gt(length(acc), 0)
    expect_true(all(sapply(acc, function(a) acpbvar:::check_sign(a$impact, S, 3L))))
})

# --- restrictions beyond impact -------------------------------------------

test_that("a matrix S with horizons > 0 is imposed at every horizon 0..horizons", {
    sim  <- simulate_var_example(T = 200, seed = 42)
    mcmc <- post_sim(sim$Y, p = 1L, R = 60L, seed = 1)
    S <- rbind(supply = c(1, -1, 0), demand = c(1, 1, 0), monetary = c(-1, -1, 1))
    set.seed(1)
    acc <- sign_restrict(mcmc, S, horizons = 3L, max_tries = 5000L, verbose = FALSE)
    expect_gt(length(acc), 0)
    for (a in acc) {
        ir <- compute_irf(a$impact, a$B_rf, n = 3L, p = 1L, H = 3L)
        for (h in 0:3) expect_true(acpbvar:::check_sign(ir[, h + 1, ], S, 3L))
    }
})

test_that("horizons > 0 is at least as demanding as impact only", {
    sim  <- simulate_var_example(T = 200, seed = 42)
    mcmc <- post_sim(sim$Y, p = 1L, R = 60L, seed = 1)
    S <- rbind(c(1, -1, 0), c(1, 1, 0), c(-1, -1, 1))
    set.seed(1); a0 <- sign_restrict(mcmc, S, verbose = FALSE)
    set.seed(1); a3 <- sign_restrict(mcmc, S, horizons = 3L, verbose = FALSE)
    expect_gte(mean(attr(a3, "n_tries")), mean(attr(a0, "n_tries")))
})

test_that("a 3-d array S imposes different signs at different horizons", {
    sim  <- simulate_var_example(T = 200, seed = 42)
    mcmc <- post_sim(sim$Y, p = 1L, R = 60L, seed = 1)
    S <- rbind(c(1, -1, 0), c(1, 1, 0), c(-1, -1, 1))
    A <- array(0, c(3, 3, 3))
    A[, , 1] <- S          # full pattern on impact
    A[3, 3, 2] <- 1        # h = 1: rate still up after a monetary shock
    A[1, 1, 3] <- 1        # h = 2: output still up after a supply shock
    set.seed(1)
    acc <- sign_restrict(mcmc, A, horizons = 2L, max_tries = 5000L, verbose = FALSE)
    expect_gt(length(acc), 0)
    for (a in acc) {
        ir <- compute_irf(a$impact, a$B_rf, n = 3L, p = 1L, H = 2L)
        expect_true(acpbvar:::check_sign(a$impact, S, 3L))
        expect_gte(ir[3, 2, 3], -1e-12)
        expect_gte(ir[1, 3, 1], -1e-12)
    }
})

test_that("partial identification works with horizons > 0 and leaves the free column free", {
    sim  <- simulate_var_example(T = 200, seed = 42)
    mcmc <- post_sim(sim$Y, p = 1L, R = 60L, seed = 1)
    S_part <- rbind(supply = c(1, -1, 0), demand = c(1, 1, 0))
    set.seed(1)
    acc <- sign_restrict(mcmc, S_part, horizons = 2L, max_tries = 5000L, verbose = FALSE)
    expect_gt(length(acc), 0)
    signs_third <- sapply(acc, function(a) sign(a$impact[1, 3]))
    expect_true(any(signs_third > 0) && any(signs_third < 0))   # third column unrestricted
    for (a in acc) {
        ir <- compute_irf(a$impact, a$B_rf, n = 3L, p = 1L, H = 2L)
        for (h in 0:2) expect_true(acpbvar:::check_sign(ir[, h + 1, ], S_part, 2L))
    }
})
