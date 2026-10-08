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
#' @param S An `n_id x n_var` matrix of sign restrictions.
#' @param n_id Number of identified shocks (defaults to `nrow(S)`).
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
#' @export
sign_restrict <- function(mcmc, S, n_id = nrow(S), max_tries = 10000, verbose = TRUE){

    check_mcmc(mcmc)
    check_S(S, mcmc$spec$n)
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
                if (check_sign(impact_f, S, n_id)) {
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
