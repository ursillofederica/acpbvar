#' Posterior quantities of one equation
#'
#' Computes, once per equation, everything the Normal-Inverse-Gamma posterior
#' needs: the Cholesky factor of the posterior precision, the posterior mean
#' of the coefficient vector, and the parameters of the inverse-gamma
#' posterior of the variance. None of these depends on the draw, so they are
#' computed once and reused by [draw()] for every posterior draw.
#'
#' @param y Response vector for the equation.
#' @param X Regressor matrix for the equation.
#' @param nu Shape hyperparameter of the prior on `sigma^2`.
#' @param m Prior mean of the structural coefficient vector.
#' @param V Diagonal of the prior covariance of the structural coefficient
#'   vector (passed as a numeric vector of prior variances).
#' @param S Scale hyperparameter of the prior on `sigma^2`.
#'
#' @return A list with `UK` (upper Cholesky factor of the posterior
#'   precision), `theta_hat` (posterior mean), `S_post` and `nu_post`
#'   (inverse-gamma posterior parameters).
#'
#' @keywords internal
posterior_eq <- function(y, X, nu, m, V, S) {

    Vinv <- 1 / V

    # K = V^{-1} + X'X ; theta_hat = K^{-1} (V^{-1} m + X'y)
    K  <- crossprod(X) + diag(Vinv)
    UK <- chol(K)
    q  <- Vinv * m + drop(crossprod(X, y))
    theta_hat <- backsolve(UK, forwardsolve(t(UK), q))

    S_post  <- S + 0.5 * (drop(crossprod(y)) + sum(m * Vinv * m) -
                          drop(crossprod(theta_hat, K %*% theta_hat)))
    nu_post <- nu + nrow(X) / 2

    list(UK = UK, theta_hat = theta_hat, S_post = S_post, nu_post = nu_post)
}


#' Single-equation Normal-Inverse-Gamma posterior draw
#'
#' Draws one value of the structural variance and of the coefficient vector
#' of a single equation from their joint Normal-Inverse-Gamma posterior,
#' using the quantities precomputed by [posterior_eq()]: `sigma^2` from its
#' marginal inverse-gamma posterior, then `theta | sigma^2` from the
#' conditional normal through the stored Cholesky factor (no inversion).
#'
#' @param post The list returned by [posterior_eq()] for the equation.
#'
#' @return A list with `sigma2` and `theta`.
#'
#' @keywords internal
draw <- function(post) {

    sigma_draw <- 1 / stats::rgamma(1, shape = post$nu_post, rate = post$S_post)

    z <- stats::rnorm(length(post$theta_hat))
    theta_draw <- post$theta_hat + sqrt(sigma_draw) * backsolve(post$UK, z)

    list(sigma2 = sigma_draw, theta = theta_draw)
}


#' Equation-by-equation posterior draw
#'
#' Loops [draw()] across all equations and collects the per-equation posterior
#' draws into list/vector containers.
#'
#' @param post_list List with one element per equation, each produced by
#'   [posterior_eq()].
#'
#' @return A list with `sigma2` (numeric vector of length `n`) and `theta`
#'   (list of length `n` of coefficient vectors).
#'
#' @keywords internal
draw_all <- function(post_list) {
    n <- length(post_list)

    sigma_draws <- numeric(n)
    theta_draws <- vector("list", n)

    for (i in seq_len(n)) {
        out <- draw(post_list[[i]])
        sigma_draws[i]   <- out$sigma2
        theta_draws[[i]] <- out$theta
    }

    list(sigma2 = sigma_draws, theta = theta_draws)
}
