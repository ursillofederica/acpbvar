#' Posterior simulation for the asymmetric A-form BVAR
#'
#' End-to-end wrapper: prepares the equation-by-equation responses and
#' regressors, fits univariate AR variance estimates, optimises the Minnesota
#' tightness via Empirical Bayes, constructs the asymmetric prior, and runs
#' the conjugate Gibbs sampler equation by equation for `R` iterations.
#'
#' @param Y A `T x n` numeric matrix of observations on the endogenous
#'   variables.
#' @param p VAR lag order.
#' @param R Number of posterior draws to retain.
#' @param kappa3 Intercept prior variance, passed to [build_prior_asym()].
#' @param unit_root_mean Logical; whether to centre the first own-lag prior
#'   at one. Default `FALSE`.
#' @param var_names Optional character vector of variable names of length
#'   `n`. Stored in the returned spec.
#' @param seed Integer seed for reproducibility. Default `123456`.
#'
#' @return A list with `samples` (the posterior draws of `theta` and
#'   `sigma2`), `spec` (sampler configuration, optimised tightness, log
#'   marginal likelihood, the `convergence` code returned by [stats::optim()]
#'   for the hyperparameter search, residual variances, variable names), and `dat`
#'   (the prepared response/regressor objects).
#'
#' @export
post_sim <- function(Y, p, R,
                     kappa3 = 100, unit_root_mean = FALSE,
                     var_names = NULL, seed = 123456) {

    n     <- ncol(Y)
    set.seed(seed)

    dat   <- build_yX(Y, p)
    s2out <- build_s2(dat$y_list)

    opt <- optim_hyper(dat$y_list, dat$X_list, p, s2out$s2,
                       kappa3, unit_root_mean)
    if (opt$convergence != 0)
        warning(sprintf("optim_hyper did not converge (code %d); kappa estimates may be unreliable.",
                        opt$convergence), call. = FALSE)

    prior <- build_prior_asym(dat$y_list, p, s2out$s2,
                              opt$kappa1, opt$kappa2, kappa3, unit_root_mean)


    samples <- list(
        theta  = vector("list", R),
        sigma2 = matrix(NA, R, n)
    )

    show_pb <- isTRUE(getOption("mcmc.pb", TRUE)) && interactive()
    if (show_pb) {
        pb <- utils::txtProgressBar(min = 0, max = R, style = 3)
        on.exit(try(close(pb), silent = TRUE), add = TRUE)
    }

    for (i in seq_len(R)) {

        post <- draw_all(dat$y_list, dat$X_list,
                         prior$nu, prior$m, prior$V, prior$S)

        samples$theta[[i]]  <- post$theta
        samples$sigma2[i, ] <- post$sigma2

        if (show_pb && ((i %% 10L) == 0L || i == R)) {
            utils::setTxtProgressBar(pb, i)
            utils::flush.console()
        }
    }
    list(
        samples = samples,
        spec    = list(
            n         = n,
            p         = p,
            TT        = dat$TT,
            R         = R,
            kappa     = c(kappa1 = opt$kappa1, kappa2 = opt$kappa2),
            log_ml    = opt$log_ml,
            convergence = opt$convergence,
            s2        = s2out$s2,
            var_names = var_names
        ),
        dat = dat
    )
}
