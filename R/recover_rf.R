#' Recover reduced-form quantities from a structural A-form draw
#'
#' Given a posterior draw of the structural coefficients and variances,
#' reconstructs the contemporaneous impact matrix `A`, its inverse, the
#' reduced-form coefficient matrix `B_rf = A^{-1} B_struct`, and the reduced-
#' form innovation covariance `Sigma_rf = A^{-1} diag(sigma2) A^{-1}'`.
#'
#' @param theta_draw A list of length `n` containing the equation-wise
#'   structural coefficient draws produced by [draw_all()].
#' @param sigma2_draw Numeric vector of length `n` of structural variance
#'   draws.
#'
#' @return A list with `A`, `A_inv`, `B_struct`, `B_rf`, and `Sigma_rf`.
#'
#' @export
recover_rf <- function(theta_draw, sigma2_draw) {

    n      <- length(theta_draw)
    k_beta <- length(theta_draw[[1]])

    A        <- diag(n)
    B_struct <- matrix(0, n, k_beta)

    for (i in seq_len(n)) {
        th      <- theta_draw[[i]]
        n_alpha <- i - 1

        if (n_alpha > 0) {
            A[i, 1:n_alpha] <- th[1:n_alpha]
        }
        B_struct[i, ] <- th[(n_alpha + 1):length(th)]
    }

    A_inv    <- solve(A)
    B_rf     <- A_inv %*% B_struct
    Sigma_rf <- A_inv %*% diag(sigma2_draw) %*% t(A_inv)

    list(A = A, A_inv = A_inv, B_struct = B_struct,
         B_rf = B_rf, Sigma_rf = Sigma_rf)
}
