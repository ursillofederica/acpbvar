#' Closed-form log marginal likelihood under the asymmetric Minnesota prior
#'
#' Returns the log marginal likelihood of the data under the asymmetric A-form
#' prior, exploiting Normal-Inverse-Gamma conjugacy equation by equation. The
#' computation reuses the Cholesky factor of the posterior precision matrix
#' to avoid forming explicit matrix inverses, and accumulates the per-equation
#' contributions into the joint log marginal likelihood.
#'
#' @inheritParams build_prior_asym
#' @param X_list Per-equation regressor matrices as produced by [build_yX()].
#'
#' @return Numeric scalar, the joint log marginal likelihood.
#'
#' @keywords internal
log_mlik_asym <- function(y_list, X_list, p, s2, kappa1, kappa2,
                          kappa3 = 100, unit_root_mean = FALSE) {

    prior <- build_prior_asym(y_list, p, s2, kappa1, kappa2, kappa3, unit_root_mean)
    m_i   <- prior$m
    V_i   <- prior$V
    nu_i  <- prior$nu
    S_i   <- prior$S

    n       <- length(y_list)
    log_ml  <- 0

    for (i in seq_len(n)) {

        y  <- y_list[[i]]
        X  <- X_list[[i]]
        TT <- nrow(X)

        Vinv <- 1 / V_i[[i]]
        m    <- m_i[[i]]

        # K_theta = V^{-1} + X'X
        K  <- crossprod(X) + diag(Vinv)
        UK <- chol(K)

        q         <- Vinv * m + drop(crossprod(X, y))
        z         <- forwardsolve(t(UK), q)
        theta_hat <- backsolve(UK, z)

        yy         <- drop(crossprod(y))
        m_Vinv_m   <- sum(m * Vinv * m)
        correction <- drop(crossprod(theta_hat, K %*% theta_hat))
        S_hat      <- S_i[i] + 0.5 * (yy + m_Vinv_m - correction)

        nu_post <- nu_i[i] + TT / 2

        log_det_Vinv <- sum(log(Vinv))
        log_det_K    <- 2 * sum(log(diag(UK)))

        log_ml_i <- -0.5 * TT * log(2 * pi) +
                     0.5 * log_det_Vinv -
                     0.5 * log_det_K +
                     lgamma(nu_post) - lgamma(nu_i[i]) +
                     nu_i[i] * log(S_i[i]) - nu_post * log(S_hat)

        log_ml <- log_ml + log_ml_i
    }

    log_ml
}


#' Empirical Bayes selection of Minnesota prior tightness
#'
#' Selects the own-lag and cross-lag tightness parameters
#' \eqn{(\kappa_1, \kappa_2)} by maximising the closed-form log marginal
#' likelihood. Optimisation is performed by Nelder-Mead in log-space, with the
#' transformation `kappa = exp(log_kappa)` ensuring positivity; values outside
#' the unit interval are penalised with a large constant.
#'
#' @inheritParams log_mlik_asym
#' @param init Numeric vector of length 2; starting values for
#'   `(kappa1, kappa2)`. Default `c(0.04, 0.004)`.
#'
#' @return A list with `kappa1`, `kappa2`, the optimised `log_ml`, and the
#'   `convergence` code returned by [stats::optim()].
#'
#' @export
optim_hyper <- function(y_list, X_list, p, s2,
                        kappa3 = 100, unit_root_mean = FALSE,
                        init = c(0.04, 0.004)) {

    neg_log_ml <- function(log_kappa) {
        k1 <- exp(log_kappa[1])
        k2 <- exp(log_kappa[2])
        if (k1 >= 1 || k2 >= 1 || k1 <= 0 || k2 <= 0) return(1e10)
        -log_mlik_asym(y_list, X_list, p, s2, k1, k2, kappa3, unit_root_mean)
    }

    opt <- stats::optim(log(init), neg_log_ml, method = "Nelder-Mead",
                        control = list(maxit = 5000))

    list(
        kappa1  = exp(opt$par[1]),
        kappa2  = exp(opt$par[2]),
        log_ml  = -opt$value,
        convergence = opt$convergence
    )
}
