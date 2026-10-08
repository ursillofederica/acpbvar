#' Build the VAR companion-form matrix
#'
#' Constructs the `np x np` companion matrix from the reduced-form coefficient
#' matrix, dropping the intercept column. Used internally by [compute_irf()].
#'
#' @param B_rf Reduced-form coefficient matrix as returned by [recover_rf()],
#'   with the first column the intercept and the remaining `n * p` columns
#'   the lag coefficients in vector form.
#' @param n Number of endogenous variables.
#' @param p VAR lag order.
#'
#' @return The `(n * p) x (n * p)` companion matrix.
#'
#' @keywords internal
build_companion <- function(B_rf, n, p) {

    B <- B_rf[, -1, drop = FALSE]
    np <- n * p
    F_comp <- matrix(0, np, np)
    F_comp[1:n, ] <- B
    if (p > 1) {
        F_comp[(n + 1):np, 1:(np - n)] <- diag(n * (p - 1))
    }
    F_comp
}


#' Compute structural impulse responses for a single draw
#'
#' Iterates the companion form for `H + 1` periods to deliver the structural
#' impulse responses for the given impact matrix and reduced-form coefficient
#' matrix.
#'
#' @param impact The `n x n` impact matrix; columns index shocks.
#' @param B_rf Reduced-form coefficient matrix as returned by [recover_rf()].
#' @param n,p Number of variables and VAR lag order.
#' @param H Maximum impulse-response horizon.
#'
#' @return A three-dimensional array of dimension `n x (H + 1) x n`, indexed
#'   by `[variable, horizon, shock]`.
#'
#' @export
compute_irf <- function(impact, B_rf, n, p, H) {
    if (!is.matrix(impact) || !all(dim(impact) == c(n, n)))
        stop("`impact` must be an n x n matrix with n = ", n, ".", call. = FALSE)
    if (!is.matrix(B_rf) || nrow(B_rf) != n || ncol(B_rf) != 1 + n * p)
        stop("`B_rf` must be an n x (1 + n * p) matrix (intercept first); got ",
             nrow(B_rf), " x ", ncol(B_rf), ".", call. = FALSE)
    check_count(H, "H", min = 0L)
    # impact: n x n (columns = shocks)
    # Returns: array n x (H+1) x n  [variable, horizon, shock]

    F_comp <- build_companion(B_rf, n, p)
    np     <- n * p
    J      <- matrix(0, np, n)
    J[1:n, ] <- diag(n)

    irf <- array(0, dim = c(n, H + 1, n))

    F_power <- diag(np)
    for (h in 0:H) {
        Phi_h <- t(J) %*% F_power %*% J
        for (k in seq_len(n)) {
            irf[, h + 1, k] <- Phi_h %*% impact[, k]
        }
        F_power <- F_power %*% F_comp
    }

    irf
}


#' Collect impulse responses across accepted draws
#'
#' Loops [compute_irf()] over the list of accepted rotations returned by
#' [sign_restrict()] and stacks the resulting impulse-response arrays along a
#' fourth dimension indexing the draw.
#'
#' @param accepted The list of accepted draws returned by [sign_restrict()].
#' @param n,p,H As in [compute_irf()].
#'
#' @return A four-dimensional array of dimension
#'   `n x (H + 1) x n x n_accepted`.
#'
#' @export
collect_irfs <- function(accepted, n, p, H) {

    check_accepted(accepted)
    check_count(n, "n", min = 1L)
    check_count(p, "p", min = 1L)
    check_count(H, "H", min = 0L)
    if (!all(dim(accepted[[1]]$impact) == c(n, n)))
        stop("`n` = ", n, " does not match the ", nrow(accepted[[1]]$impact),
             " x ", ncol(accepted[[1]]$impact), " impact matrices in `accepted`.", call. = FALSE)
    if (ncol(accepted[[1]]$B_rf) != 1 + n * p)
        stop("`p` = ", p, " does not match the lag structure in `accepted` (B_rf has ",
             ncol(accepted[[1]]$B_rf), " columns, expected 1 + n * p = ", 1 + n * p, ").", call. = FALSE)
    # Returns: array n x (H+1) x n x n_accepted

    n_acc <- length(accepted)
    irfs  <- array(NA, dim = c(n, H + 1, n, n_acc))

    for (a in seq_len(n_acc)) {
        irfs[, , , a] <- compute_irf(
            accepted[[a]]$impact,
            accepted[[a]]$B_rf,
            n, p, H
        )
    }

    irfs
}
