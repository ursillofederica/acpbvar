#' Plot posterior impulse responses with credible bands
#'
#' Produces a grid of base-R plots displaying the posterior median impulse
#' response of each variable to each identified shock, together with a
#' shaded credible band defined by the quantiles in `probs`.
#'
#' @param irfs A four-dimensional array of impulse responses as returned by
#'   [collect_irfs()].
#' @param var_names Optional character vector of variable names of length
#'   `n`; defaults to `"Var i"`.
#' @param shock_names Optional character vector of shock names of length
#'   `n_shocks`; defaults to `"Shock k"`.
#' @param probs Two-element numeric vector of quantiles for the credible
#'   band. Default `c(0.16, 0.84)`.
#' @param col_med Colour of the posterior median line. Default `"black"`.
#' @param col_band Colour of the credible band. Default `"steelblue"`.
#' @param alpha_band Transparency of the credible band. Default `0.3`.
#' @param ylab Y-axis label. Default empty.
#' @param main_cex Scaling factor for panel titles. Default `0.9`.
#'
#' @return Called for its side effect (a grid of base-R plots). Invisibly
#'   returns `NULL`.
#'
#' @export
plot_irfs <- function(irfs, var_names = NULL, shock_names = NULL,
                      probs = c(0.16, 0.84), col_med = "black",
                      col_band = "steelblue", alpha_band = 0.3,
                      ylab = "", main_cex = 0.9) {

    if (!is.array(irfs) || length(dim(irfs)) != 4)
        stop("`irfs` must be the 4-dimensional array returned by collect_irfs() ",
             "(variable x horizon x shock x draw).", call. = FALSE)
    if (length(probs) != 2 || !is.numeric(probs) || any(probs < 0 | probs > 1) || probs[1] >= probs[2])
        stop("`probs` must be two increasing values in [0, 1], e.g. c(0.16, 0.84).", call. = FALSE)
    if (!is.null(var_names) && length(var_names) != dim(irfs)[1])
        stop("`var_names` must have length ", dim(irfs)[1], ".", call. = FALSE)
    if (!is.null(shock_names) && length(shock_names) != dim(irfs)[3])
        stop("`shock_names` must have length ", dim(irfs)[3], ".", call. = FALSE)

    n     <- dim(irfs)[1]
    H1    <- dim(irfs)[2]       # H + 1
    n_sh  <- dim(irfs)[3]
    horiz <- 0:(H1 - 1)

    if (is.null(var_names))   var_names   <- paste0("Var ", seq_len(n))
    if (is.null(shock_names)) shock_names <- paste0("Shock ", seq_len(n_sh))

    irf_med <- apply(irfs, c(1, 2, 3), stats::median)
    irf_lo  <- apply(irfs, c(1, 2, 3), stats::quantile, probs = probs[1])
    irf_hi  <- apply(irfs, c(1, 2, 3), stats::quantile, probs = probs[2])

    col_fill <- grDevices::adjustcolor(col_band, alpha.f = alpha_band)

    op <- graphics::par(mfrow = c(n, n_sh), mar = c(2.5, 2.5, 1.8, 0.5),
                        mgp = c(1.4, 0.4, 0), tcl = -0.25, cex.axis = 0.75)
    on.exit(graphics::par(op))

    for (j in seq_len(n)) {
        for (k in seq_len(n_sh)) {
            ylim <- range(irf_lo[j, , k], irf_hi[j, , k])
            plot(horiz, irf_med[j, , k], type = "n", ylim = ylim,
                 xlab = "", ylab = ylab)
            graphics::polygon(c(horiz, rev(horiz)),
                              c(irf_lo[j, , k], rev(irf_hi[j, , k])),
                              col = col_fill, border = NA)
            graphics::lines(horiz, irf_med[j, , k], col = col_med, lwd = 1.5)
            graphics::abline(h = 0, lty = 2, col = "grey50")
            graphics::title(paste0(var_names[j], " <- ", shock_names[k]),
                            cex.main = main_cex)
        }
    }

    invisible(NULL)
}
