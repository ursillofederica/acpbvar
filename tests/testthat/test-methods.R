# S3 classes and methods: the objects stay the documented list/array; the
# methods only read them. Small model so that the file runs in seconds.

sim  <- simulate_var_example(T = 120, seed = 3)
mcmc <- post_sim(sim$Y, p = 1L, R = 150L, var_names = sim$var_names, seed = 3)
S    <- rbind(supply = c(1, -1, 0), demand = c(1, 1, 0), monetary = c(-1, -1, 1))
acc  <- sign_restrict(mcmc, S, max_tries = 500L, verbose = FALSE)
irfs <- collect_irfs(acc, n = 3L, p = 1L, H = 6L)

test_that("returned objects carry their class and remain the documented structures", {
    expect_s3_class(mcmc, "acpbvar_post")
    expect_true(is.list(mcmc) && all(c("samples", "spec", "dat") %in% names(mcmc)))
    expect_s3_class(acc, "acpbvar_signs")
    expect_true(is.list(acc) && all(c("impact", "B_rf") %in% names(acc[[1]])))
    expect_s3_class(irfs, "acpbvar_irf")
    expect_equal(dim(irfs), c(3, 7, 3, length(acc)))
    expect_equal(dimnames(irfs)[[1]], sim$var_names)
    expect_equal(dimnames(irfs)[[3]], rownames(S))
    expect_equal(dimnames(irfs)[[2]], paste0("h", 0:6))
    # subsetting gives a plain array, as before
    expect_false(inherits(irfs[1, , 1, ], "acpbvar_irf"))
    expect_equal(dim(irfs[, , 1, ]), c(3, 7, length(acc)))
})

test_that("names are padded under partial identification", {
    acc2  <- sign_restrict(mcmc, S[1:2, ], max_tries = 500L, verbose = FALSE)
    irfs2 <- collect_irfs(acc2, n = 3L, p = 1L, H = 2L)
    expect_equal(dimnames(irfs2)[[3]], c("supply", "demand", "Shock 3"))
})

test_that("print methods run and mention the key quantities", {
    expect_output(print(mcmc), "3 variables, 1 lags")
    expect_output(print(mcmc), "kappa1")
    expect_output(print(acc), "accepted")
    expect_output(print(acc), "rotations per accepted draw")
    expect_output(print(irfs), "3 variables x 6 horizons")
    expect_invisible(print(irfs))
})

test_that("summary returns the same quantiles as a manual computation", {
    s <- summary(irfs, vars = "output", shocks = "monetary", horizons = c(0, 3, 6))
    expect_s3_class(s, "summary.acpbvar_irf")
    manual <- apply(irfs[1, c(1, 4, 7), 3, ], 1, quantile, probs = c(0.16, 0.5, 0.84))
    expect_equal(unname(s$monetary$output), unname(manual))
    expect_equal(colnames(s$monetary$output), c("h0", "h3", "h6"))
    expect_output(print(s), "Response of output to monetary")
    # integer selection and default horizons
    s2 <- summary(irfs, vars = 2L, shocks = 1L)
    expect_equal(dim(s2$supply$inflation), c(3, 7))
    expect_error(summary(irfs, horizons = 9), "between 0 and 6")
    expect_error(summary(irfs, vars = "price"), "not found")
})

test_that("plot method and plot_irfs selection and overlay run without error", {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off())
    expect_invisible(plot(irfs))
    expect_invisible(plot(irfs, shocks = "monetary", mfrow = c(1, 3), col_band = "grey85"))
    expect_invisible(plot_irfs(irfs, vars = c(1, 3), shocks = 2:3))
    # 4-d reference: true responses as a one-draw array
    irf_true <- compute_irf(sim$impact_true, cbind(sim$intercept_true, sim$B_true), 3L, 1L, 6L)
    expect_invisible(plot(irfs, ref = array(irf_true, c(dim(irf_true), 1))))
    # 3-d reference with a single shock
    expect_invisible(plot(irfs, shocks = 1, ref = irfs[, , 1, ]))
    expect_error(plot(irfs, ref = irfs[, , 1, ]), "single shock")
    expect_error(plot(irfs, vars = 7), "between 1 and 3")
    expect_error(plot(irfs, probs = c(0.9, 0.1)), "increasing")
})
