#' Verify a candidate impact matrix against a sign-restriction pattern
#'
#' Checks whether the column `k` of the candidate impact matrix satisfies the
#' sign pattern `S[k, ]`, for `k = 1, ..., n_id`. Entries of `S` equal to `+1`,
#' `-1`, or `0` denote respectively a positive, negative, or unrestricted sign
#' for the corresponding variable response.
#'
#' @param impact The candidate `n_var x n_shock` impact matrix.
#' @param S An `n_id x n_var` matrix of sign restrictions.
#' @param n_id Integer; number of identified shocks (defaults to `nrow(S)`).
#' @param tol Numeric tolerance below which violations are ignored.
#'
#' @return `TRUE` if the candidate satisfies the restrictions, otherwise
#'   `FALSE`.
#'
#' @keywords internal
check_sign <- function(impact, S, n_id, tol = 1e-12) {
    # impact: n_var x n_shock
    # S:      n_id  x n_var
    # S[k, j] = required sign for variable j under shock k
    # impact[j, k] = response of variable j to shock k
    for (k in seq_len(n_id)) {
        for (j in seq_len(ncol(S))) {
            if (S[k, j] ==  1 && impact[j, k] < -tol) return(FALSE)
            if (S[k, j] == -1 && impact[j, k] >  tol) return(FALSE)
        }
    }
    TRUE
}


#' Structural identification via sign restrictions (RRWZ algorithm)
#'
#' Implements the Rubio-Ramirez, Waggoner and Zha (2010, Theorem 9)
#' rotation algorithm for sign-restriction identification. For each retained
#' posterior draw of the reduced form, candidate rotations `Q` are drawn from
#' the Haar measure (via QR decomposition of a standard normal matrix, with
#' sign normalisation on the diagonal of `R`); for each candidate, all
#' `2^n_id` sign flips on the columns are checked against the restrictions,
#' and the first accepted rotation is stored.
#'
#' @param mcmc A posterior object produced by [post_sim()].
#' @param S Sign restrictions, with entries `+1`, `-1` or `0` (unrestricted).
#'   Either an `n_id x n_var` matrix (rows: identified shocks, columns:
#'   variables), imposed on the impulse responses at every horizon
#'   `0, ..., horizons`; or an `n_id x n_var x (horizons + 1)` array with one
#'   matrix per horizon, so that different signs can be imposed at different
#'   horizons (a slice of zeros leaves that horizon unrestricted). Every
#'   identified shock must carry at least one non-zero restriction at some
#'   horizon. A matrix with `horizons = 0` restricts the impact only.
#' @param n_id Number of identified shocks (defaults to `nrow(S)`). Columns of
#'   the rotation beyond `n_id` are left unrestricted and unlabelled.
#' @param horizons Non-negative integer. Sign restrictions are imposed on the
#'   impulse responses at horizons `0, 1, ..., horizons`; the default `0L`
#'   restricts the impact matrix only. Responses at later horizons are
#'   computed with [compute_irf()] from the reduced form of each draw, so the
#'   cost per candidate rotation grows with `horizons`, and draws whose
#'   identified set is empty cost `max_tries` candidates each.
#' @param max_tries Maximum number of `Q` candidates per posterior draw before
#'   moving on. Default `10000`. Draws for which no rotation
#'   is found are dropped from the output and listed in the `failed_draws`
#'   attribute. This is the intended behaviour when the identified set of a
#'   draw is empty, but it also removes draws whose identified set is small
#'   relative to `max_tries`. Check that `max(attr(x, "n_tries"))` is well
#'   below `max_tries`: if it is close to it, the limit is binding and
#'   `max_tries` should be increased.
#' @param verbose Logical; print a progress bar and acceptance summary.
#'
#' @return A list of length equal to the number of accepted draws; each entry
#'   contains `draw_id`, `impact`, `Q`, `B_rf`, `Sigma_rf`.
#'   The result carries two attributes: `"n_tries"`, the number of `Q`
#'   candidates needed for each accepted draw, and `"failed_draws"`, the
#'   indices of the posterior draws dropped after `max_tries`.
#'
#' @examples
#' sim  <- simulate_var_example(T = 120, seed = 1)
#' mcmc <- post_sim(sim$Y, p = 1L, R = 200L)
#'
#' # Rows: identified shocks. Columns: variables (output, inflation, rate).
#' # +1 / -1 restrict the sign of the impact response, 0 leaves it free.
#' S <- rbind(supply   = c( 1, -1,  0),
#'            demand   = c( 1,  1,  0),
#'            monetary = c(-1, -1,  1))
#'
#' acc <- sign_restrict(mcmc, S, max_tries = 500L, verbose = FALSE)
#' length(acc)                           # accepted draws (one rotation each)
#' acc[[1]]$impact                       # impact matrix of the first draw
#' summary(attr(acc, "n_tries"))         # rotations screened per acceptance;
#'                                       # compare the max with max_tries
#' attr(acc, "failed_draws")             # draws dropped after max_tries
#'
#' # Partial identification: two shocks, the third column is left free
#' S2  <- S[1:2, ]
#' acc2 <- sign_restrict(mcmc, S2, max_tries = 500L, verbose = FALSE)
#'
#' # Same signs imposed on impact and on the next three horizons (as in
#' # Uhlig 2005): more demanding, so more rotations are screened per draw
#' acc3 <- sign_restrict(mcmc, S, horizons = 3L, max_tries = 2000L, verbose = FALSE)
#' summary(attr(acc3, "n_tries"))
#'
#' # Different signs at different horizons: one matrix per horizon
#' A <- array(0, c(3, 3, 3))
#' A[, , 1] <- S            # impact: full pattern
#' A[3, 3, 2] <- 1          # h = 1: interest rate still up after a monetary shock
#' A[1, 1, 3] <- 1          # h = 2: output still up after a supply shock
#' accA <- sign_restrict(mcmc, A, horizons = 2L, max_tries = 2000L, verbose = FALSE)
#' @export
sign_restrict <- function(mcmc, S, n_id = nrow(S), horizons = 0L,
                          max_tries = 10000, verbose = TRUE){

    check_mcmc(mcmc)
    check_count(horizons, "horizons", min = 0L)
    S <- check_S(S, mcmc$spec$n, horizons)
    check_count(n_id, "n_id", min = 1L, max = nrow(S))
    check_count(max_tries, "max_tries", min = 1L)
    if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose))
        stop("`verbose` must be TRUE or FALSE.", call. = FALSE)

    R_tot <- mcmc$spec$R
    n     <- mcmc$spec$n
    p     <- mcmc$spec$p

    n_flips <- expand.grid(rep(list(c(-1, 1)), n_id))
    flip_list <- lapply(seq_len(nrow(n_flips)), function(f){
        d <- rep(1L, n)
        d[seq_len(n_id)] <- as.integer(n_flips[f, ])
        d
    })
    # one sign matrix per horizon 0, ..., horizons (kept as matrices even if n_id = 1)
    S_h <- lapply(seq_len(horizons + 1), function(h) matrix(S[, , h], nrow = n_id))

    accepted  <- list()
    n_tries_v <- integer(0) # Q draws needed per accepted draw
    failed    <- integer(0)

    if (verbose) {
        pb <- utils::txtProgressBar(min = 0, max = R_tot, style = 3)
        on.exit(close(pb))
    }

    for (s in seq_len(R_tot)) {
        # recover reduced form
        rf       <- recover_rf(mcmc$samples$theta[[s]], mcmc$samples$sigma2[s, ])
        Sigma_rf <- rf$Sigma_rf
        B_rf     <- rf$B_rf

        # Cholesky factor (P0)
        P0 <- t(chol(Sigma_rf))

        found <- FALSE

        # orthogonal Q (Theorem 9 RRWZ)
        for (i in seq_len(max_tries)) {
            X <- matrix(stats::rnorm(n * n), n, n)
            qrX <- qr(X)
            Q <- qr.Q(qrX)
            R_qr <- qr.R(qrX)
            d <- sign(diag(R_qr))
            d[d == 0] <- 1
            Q <- Q %*% diag(d)

            # impact matrix
            P <- P0 %*% Q

            # check all the 2^n sign flips
            for (f in seq_along(flip_list)) {
                impact_f <- P %*% diag(flip_list[[f]])

                # impact first; responses at later horizons only if impact passes
                ok <- check_sign(impact_f, S_h[[1]], n_id)
                if (ok && horizons > 0L) {
                    irf_f <- compute_irf(impact_f, B_rf, n, p, H = horizons)
                    for (h in seq_len(horizons)) {
                        if (!check_sign(irf_f[, h + 1, ], S_h[[h + 1]], n_id)) {
                            ok <- FALSE
                            break
                        }
                    }
                }
                if (ok) {
                    accepted[[length(accepted) + 1]] <- list(
                        draw_id  = s,
                        impact   = impact_f,
                        Q        = Q %*% diag(flip_list[[f]]),
                        B_rf     = B_rf,
                        Sigma_rf = Sigma_rf
                    )
                    n_tries_v <- c(n_tries_v, i)
                    found <- TRUE
                    break
                }
            }
            if (found) break
        }
        if (!found) failed <- c(failed, s)
        if (verbose) utils::setTxtProgressBar(pb, s)
    }

    if (verbose) {
        cat(sprintf("\nAccepted %d / %d draws (%.1f%%); %d draws dropped after %d tries\n",
                    length(accepted), R_tot, 100 * length(accepted) / R_tot,
                    length(failed), max_tries))
        if (length(n_tries_v) > 0) {
            cat(sprintf("Q draws per acceptance: mean=%.1f, median=%d, max=%d\n",
                        mean(n_tries_v), stats::median(n_tries_v), max(n_tries_v)))
            # Implied per-Q acceptance rate (accounting for sign flips)
            p_per_flip <- 1 - (1 - 1/mean(n_tries_v))^(1/length(flip_list))
            cat(sprintf("Implied per-flip prob: %.6f  (per-Q prob: %.4f)\n",
                        p_per_flip, 1/mean(n_tries_v)))
        }


    }

    if (length(failed) > 0){
        warning(sprintf(paste0("%d of %d posterior draws dropped: no rotation found within max_tries = %d. ",
                               "If max(attr(., 'n_tries')) is close to max_tries, increase it."),
                        length(failed), R_tot, max_tries), call. = FALSE)
    }

    attr(accepted, "n_tries") <- n_tries_v
    attr(accepted, "failed_draws") <- failed
    accepted
}
