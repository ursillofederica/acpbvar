#' Single-equation conjugate Gibbs draw
#'
#' Draws posterior values for the structural variance and coefficient vector
#' of a single equation in the A-form representation, exploiting Normal-
#' Inverse-Gamma conjugacy. The draw uses the Cholesky factor of the posterior
#' precision to avoid explicit inversion and to deliver an exact joint draw.
#'
#' @param y Response vector for the equation.
#' @param X Regressor matrix for the equation.
#' @param nu Shape hyperparameter of the prior on `sigma^2`.
#' @param m Prior mean of the structural coefficient vector.
#' @param V Diagonal of the prior covariance of the structural coefficient
#'   vector (passed as a numeric vector of prior variances).
#' @param S Scale hyperparameter of the prior on `sigma^2`.
#'
#' @return A list with `sigma2`, `theta` (the draw), and the posterior
#'   summaries `theta_hat`, `S_post`, `nu_post`.
#'
#' @keywords internal
draw <- function(y, X, nu, m, V, S) {

    TT <- nrow(X)
    yy <- drop(crossprod(y))
    XX <- crossprod(X)

    Vinv <- 1 / V

    # q = V^{-1} m + X'y
    q <- Vinv * m + drop(crossprod(X, y))

    # K = V^{-1} + X'X
    K <- XX + diag(Vinv)

    # theta_hat = K^{-1} q
    UK        <- chol(K)
    z         <- forwardsolve(t(UK), q)
    theta_hat <- backsolve(UK, z)

    m_Vinv_m   <- sum(m * Vinv * m)
    correction <- drop(crossprod(theta_hat, K %*% theta_hat))

    S_post  <- S + 0.5 * (yy + m_Vinv_m - correction)
    nu_post <- nu + TT / 2

    sigma_draw <- 1 / stats::rgamma(1, shape = nu_post, rate = S_post)


    z <- stats::rnorm(length(theta_hat))
    theta_draw <- theta_hat + sqrt(sigma_draw) * backsolve(UK, z)

    list(
        sigma2    = sigma_draw,
        theta     = theta_draw,
        theta_hat = theta_hat,
        S_post    = S_post,
        nu_post   = nu_post
    )
}


#' Equation-by-equation conjugate Gibbs sweep
#'
#' Loops [draw()] across all equations and collects the per-equation posterior
#' draws into list/vector containers.
#'
#' @param y_list,X_list Per-equation responses and regressors.
#' @param nu,m_list,V_list,S Prior components produced by
#'   [build_prior_asym()].
#'
#' @return A list with vectors and lists of length `n`: `sigma2`, `theta`,
#'   `theta_hat`, `S_post`, `nu_post`.
#'
#' @keywords internal
draw_all <- function(y_list, X_list, nu, m_list, V_list, S) {
    n <- length(y_list)

    sigma_draws <- numeric(n)
    theta_draws <- vector("list", n)
    theta_hats  <- vector("list", n)
    S_posts     <- numeric(n)
    nu_posts    <- numeric(n)

    for (i in seq_len(n)) {
        out <- draw(
            y  = y_list[[i]],
            X  = X_list[[i]],
            nu = nu[i],
            m  = m_list[[i]],
            V  = V_list[[i]],
            S  = S[i]
        )

        sigma_draws[i]   <- out$sigma2
        theta_draws[[i]] <- out$theta
        theta_hats[[i]]  <- out$theta_hat
        S_posts[i]       <- out$S_post
        nu_posts[i]      <- out$nu_post
    }

    list(
        sigma2    = sigma_draws,
        theta     = theta_draws,
        theta_hat = theta_hats,
        S_post    = S_posts,
        nu_post   = nu_posts
    )
}
