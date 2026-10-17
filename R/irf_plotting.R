#' Plot posterior impulse responses with credible bands
#'
#' Produces a grid of base-R plots displaying the posterior median impulse
#' response of each selected variable to each selected shock, together with a
#' shaded credible band defined by the quantiles in `probs`. A second set of
#' responses (`ref`), for example from another package or the true responses
#' of a simulated process, can be overlaid for comparison.
#'
#' @param irfs A four-dimensional array of impulse responses as returned by
#'   [collect_irfs()] (variable x horizon x shock x draw).
#' @param vars,shocks Which variables and shocks to plot: integer indices or
#'   names matching the dimnames of `irfs`. Default `NULL` plots all.
#' @param var_names Optional character vector of variable names of length
#'   `n`; defaults to the dimnames of `irfs`, or `"Var i"`.
#' @param shock_names Optional character vector of shock names of length
#'   `n_shocks`; defaults to the dimnames of `irfs`, or `"Shock k"`.
#' @param probs Two-element numeric vector of quantiles for the credible
#'   band. Default `c(0.16, 0.84)`.
#' @param ref Optional responses to overlay (dashed median, dotted band at
#'   the same `probs`): either a four-dimensional array with the same
#'   variables and horizons as `irfs` (subset with the same `vars` and
#'   `shocks`), or a three-dimensional array (variable x horizon x draw) when
#'   a single shock is plotted.
#' @param col_med Colour of the posterior median line. Default `"black"`.
#' @param col_band Colour of the credible band. Default `"steelblue"`.
#' @param col_ref Colour of the overlaid responses. Default `"black"`.
#' @param alpha_band Transparency of the credible band. Default `0.3`.
#' @param xlab,ylab Axis labels. Defaults `""`.
#' @param main_cex Scaling factor for panel titles. Default `0.9`.
#' @param mfrow Panel layout passed to [graphics::par()]; default one row per
#'   variable and one column per shock.
#'
#' @return Called for its side effect (a grid of base-R plots). Invisibly
#'   returns `NULL`.
#'
#' @examples
#' sim  <- simulate_var_example(T = 120, seed = 1)
#' mcmc <- post_sim(sim$Y, p = 1L, R = 200L)
#' S    <- rbind(supply = c(1, -1, 0), demand = c(1, 1, 0), monetary = c(-1, -1, 1))
#' acc  <- sign_restrict(mcmc, S, max_tries = 500L, verbose = FALSE)
#' irfs <- collect_irfs(acc, n = 3L, p = 1L, H = 8L)
#'
#' plot_irfs(irfs, probs = c(0.16, 0.84))               # all panels, 68% band
#' plot_irfs(irfs, shocks = "monetary", mfrow = c(1, 3)) # one shock, one row
#'
#' # overlay the true responses of the simulated process
#' irf_true <- compute_irf(sim$impact_true, cbind(sim$intercept_true, sim$B_true),
#'                         n = 3L, p = 1L, H = 8L)
#' plot_irfs(irfs, ref = array(irf_true, c(dim(irf_true), 1)))
#' @export
plot_irfs <- function(irfs, vars = NULL, shocks = NULL,
                      var_names = NULL, shock_names = NULL,
                      probs = c(0.16, 0.84), ref = NULL,
                      col_med = "black", col_band = "steelblue", col_ref = "black",
                      alpha_band = 0.3, xlab = "", ylab = "", main_cex = 0.9,
                      mfrow = NULL) {

    if (!is.array(irfs) || length(dim(irfs)) != 4)
        stop("`irfs` must be the 4-dimensional array returned by collect_irfs() ",
             "(variable x horizon x shock x draw).", call. = FALSE)
    if (length(probs) != 2 || !is.numeric(probs) || any(probs < 0 | probs > 1) || probs[1] >= probs[2])
        stop("`probs` must be two increasing values in [0, 1], e.g. c(0.16, 0.84).", call. = FALSE)

    dn <- dimnames(irfs)
    if (is.null(var_names))   var_names   <- if (!is.null(dn[[1]])) dn[[1]] else paste0("Var ", seq_len(dim(irfs)[1]))
    if (is.null(shock_names)) shock_names <- if (!is.null(dn[[3]])) dn[[3]] else paste0("Shock ", seq_len(dim(irfs)[3]))
    if (length(var_names) != dim(irfs)[1])
        stop("`var_names` must have length ", dim(irfs)[1], ".", call. = FALSE)
    if (length(shock_names) != dim(irfs)[3])
        stop("`shock_names` must have length ", dim(irfs)[3], ".", call. = FALSE)

    vars   <- resolve_index(vars,   var_names,   "vars")
    shocks <- resolve_index(shocks, shock_names, "shocks")

    irfs <- unclass(irfs)[vars, , shocks, , drop = FALSE]
    var_names   <- var_names[vars]
    shock_names <- shock_names[shocks]

    n     <- dim(irfs)[1]
    H1    <- dim(irfs)[2]       # H + 1
    n_sh  <- dim(irfs)[3]
    horiz <- 0:(H1 - 1)

    if (!is.null(ref)) {
        if (!is.array(ref) || !(length(dim(ref)) %in% c(3, 4)))
            stop("`ref` must be a 4-d array (variable x horizon x shock x draw) or, for a single shock, ",
                 "a 3-d array (variable x horizon x draw).", call. = FALSE)
        if (length(dim(ref)) == 3) {
            if (n_sh != 1)
                stop("A 3-d `ref` can only be overlaid when a single shock is plotted.", call. = FALSE)
            ref <- array(ref, c(dim(ref)[1:2], 1, dim(ref)[3]))
            ref <- ref[vars, , , , drop = FALSE]
        } else {
            ref <- unclass(ref)[vars, , shocks, , drop = FALSE]
        }
        if (dim(ref)[2] != H1)
            stop("`ref` has ", dim(ref)[2], " horizons but `irfs` has ", H1, ".", call. = FALSE)
    }

    irf_med <- apply(irfs, c(1, 2, 3), stats::median)
    irf_lo  <- apply(irfs, c(1, 2, 3), stats::quantile, probs = probs[1])
    irf_hi  <- apply(irfs, c(1, 2, 3), stats::quantile, probs = probs[2])
    if (!is.null(ref)) {
        ref_med <- apply(ref, c(1, 2, 3), stats::median)
        ref_lo  <- apply(ref, c(1, 2, 3), stats::quantile, probs = probs[1])
        ref_hi  <- apply(ref, c(1, 2, 3), stats::quantile, probs = probs[2])
    }

    col_fill <- grDevices::adjustcolor(col_band, alpha.f = alpha_band)
    if (is.null(mfrow)) mfrow <- c(n, n_sh)

    op <- graphics::par(mfrow = mfrow, mar = c(2.5, 2.5, 1.8, 0.5),
                        mgp = c(1.4, 0.4, 0), tcl = -0.25, cex.axis = 0.75)
    on.exit(graphics::par(op))

    for (j in seq_len(n)) {
        for (k in seq_len(n_sh)) {
            ylim <- range(irf_lo[j, , k], irf_hi[j, , k],
                          if (!is.null(ref)) c(ref_lo[j, , k], ref_hi[j, , k]))
            plot(horiz, irf_med[j, , k], type = "n", ylim = ylim,
                 xlab = xlab, ylab = ylab)
            graphics::polygon(c(horiz, rev(horiz)),
                              c(irf_lo[j, , k], rev(irf_hi[j, , k])),
                              col = col_fill, border = NA)
            graphics::abline(h = 0, lty = 2, col = "grey50")
            if (!is.null(ref)) {
                graphics::lines(horiz, ref_lo[j, , k],  col = col_ref, lty = 3)
                graphics::lines(horiz, ref_hi[j, , k],  col = col_ref, lty = 3)
                graphics::lines(horiz, ref_med[j, , k], col = col_ref, lty = 2, lwd = 1.5)
            }
            graphics::lines(horiz, irf_med[j, , k], col = col_med, lwd = 1.5)
            title_k <- if (n_sh == 1) var_names[j] else paste0(var_names[j], " <- ", shock_names[k])
            graphics::title(title_k, cex.main = main_cex)
        }
    }

    invisible(NULL)
}

# integer or character selection against a vector of names (internal)
resolve_index <- function(x, names, arg) {
    if (is.null(x)) return(seq_along(names))
    if (is.character(x)) {
        idx <- match(x, names)
        if (anyNA(idx))
            stop("`", arg, "` contains names not found: ", paste(x[is.na(idx)], collapse = ", "),
                 ". Available: ", paste(names, collapse = ", "), ".", call. = FALSE)
        return(idx)
    }
    if (!is.numeric(x) || anyNA(x) || any(x != as.integer(x)) || any(x < 1 | x > length(names)))
        stop("`", arg, "` must be names or integer indices between 1 and ", length(names), ".", call. = FALSE)
    as.integer(x)
}
