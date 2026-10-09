#' Compute structural impulse responses for a single draw
#'
#' Computes the structural impulse responses at horizons `0, ..., H` for the
#' given impact matrix and reduced-form coefficient matrix, through the VMA
#' recursion on the `n x n` lag matrices.
#'
#' @param impact The `n x n` impact matrix; columns index shocks.
#' @param B_rf Reduced-form coefficient matrix as returned by [recover_rf()].
#' @param n,p Number of variables and VAR lag order.
#' @param H Maximum impulse-response horizon.
#'
#' @return A three-dimensional array of dimension `n x (H + 1) x n`, indexed
#'   by `[variable, horizon, shock]`.
#'
#' @examples
#' sim  <- simulate_var_example(T = 120, seed = 1)
#' mcmc <- post_sim(sim$Y, p = 1L, R = 50L)
#' S    <- rbind(c(1, -1, 0), c(1, 1, 0), c(-1, -1, 1))
#' acc  <- sign_restrict(mcmc, S, max_tries = 500L, verbose = FALSE)
#'
#' irf1 <- compute_irf(acc[[1]]$impact, acc[[1]]$B_rf, n = 3L, p = 1L, H = 8L)
#' dim(irf1)                   # variable x horizon (0..H) x shock
#' irf1[1, , 3]                # response of variable 1 to shock 3 over horizons
#' @export
compute_irf <- function(impact, B_rf, n, p, H) {
    if (!is.matrix(impact) || !all(dim(impact) == c(n, n)))
        stop("`impact` must be an n x n matrix with n = ", n, ".", call. = FALSE)
    if (!is.matrix(B_rf) || nrow(B_rf) != n || ncol(B_rf) != 1 + n * p)
        stop("`B_rf` must be an n x (1 + n * p) matrix (intercept first); got ",
             nrow(B_rf), " x ", ncol(B_rf), ".", call. = FALSE)
    check_count(H, "H", min = 0L)

    # VMA recursion on the n x n lag matrices: Phi_0 = I,
    # Phi_h = sum_{j=1}^{min(h,p)} Phi_{h-j} B_j. Same numbers as powering the
    # companion matrix, at a fraction of the cost when n * p is large.
    B   <- B_rf[, -1, drop = FALSE]
    Phi <- vector("list", H + 1)
    Phi[[1]] <- diag(n)
    for (h in seq_len(H)) {
        acc <- matrix(0, n, n)
        for (j in seq_len(min(h, p))) {
            acc <- acc + Phi[[h - j + 1]] %*% B[, ((j - 1) * n + 1):(j * n), drop = FALSE]
        }
        Phi[[h + 1]] <- acc
    }

    irf <- array(0, dim = c(n, H + 1, n))
    for (h in 0:H) irf[, h + 1, ] <- Phi[[h + 1]] %*% impact
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
#' @examples
#' sim  <- simulate_var_example(T = 120, seed = 1)
#' mcmc <- post_sim(sim$Y, p = 1L, R = 200L)
#' S    <- rbind(c(1, -1, 0), c(1, 1, 0), c(-1, -1, 1))
#' acc  <- sign_restrict(mcmc, S, max_tries = 500L, verbose = FALSE)
#'
#' irfs <- collect_irfs(acc, n = 3L, p = 1L, H = 8L)
#' dim(irfs)                   # variable x horizon x shock x accepted draw
#'
#' # Posterior median and 68% band of variable 1 to shock 3, horizons 0..8
#' apply(irfs[1, , 3, ], 1, quantile, probs = c(0.16, 0.5, 0.84))
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
