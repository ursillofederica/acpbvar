#' Print, summary and plot methods for acpbvar objects
#'
#' [post_sim()] returns an object of class `"acpbvar_post"`,
#' [sign_restrict()] an object of class `"acpbvar_signs"` and
#' [collect_irfs()] an object of class `"acpbvar_irf"`. All three remain the
#' plain list or array documented in the respective function (subsetting an
#' `"acpbvar_irf"` with `[` gives a plain array); the classes only add the
#' methods below.
#'
#' @param x,object An object of the corresponding class.
#' @param horizons Integer vector of horizons to summarise (default all).
#' @param probs Quantiles to report (summary, default `c(0.16, 0.5, 0.84)`)
#'   or the two band quantiles (plot, default `c(0.16, 0.84)`).
#' @param vars,shocks Variables and shocks to summarise or plot: integer
#'   indices or names. Default `NULL` means all.
#' @param digits Number of digits when printing a summary.
#' @param ... For `plot`, further arguments passed to [plot_irfs()] (for
#'   example `ref`, `col_band`, `xlab`, `mfrow`); ignored otherwise.
#'
#' @return `print` methods return their argument invisibly. `summary` returns
#'   an object of class `"summary.acpbvar_irf"`: a list, one element per
#'   selected shock, each a list of `probs x horizons` matrices, one per
#'   selected variable. `plot` returns `NULL` invisibly.
#'
#' @examples
#' sim  <- simulate_var_example(T = 120, seed = 1)
#' mcmc <- post_sim(sim$Y, p = 1L, R = 200L, var_names = sim$var_names)
#' mcmc
#' S    <- rbind(supply = c(1, -1, 0), demand = c(1, 1, 0), monetary = c(-1, -1, 1))
#' acc  <- sign_restrict(mcmc, S, max_tries = 500L, verbose = FALSE)
#' acc
#' irfs <- collect_irfs(acc, n = 3L, p = 1L, H = 8L)
#' irfs
#' summary(irfs, vars = "output", shocks = "monetary", horizons = c(0, 4, 8))
#' plot(irfs, shocks = "monetary", mfrow = c(1, 3))
#' @name acpbvar-methods
NULL

#' @rdname acpbvar-methods
#' @export
print.acpbvar_post <- function(x, ...) {
    sp <- x$spec
    cat("acpbvar posterior: ", sp$n, " variables, ", sp$p, " lags, ", sp$TT,
        " observations\n", sep = "")
    cat("  ", sp$R, " exact draws; elicitation \"", sp$elicit, "\"\n", sep = "")
    cat("  empirical Bayes: kappa1 = ", format(sp$kappa[["kappa1"]], digits = 4),
        ", kappa2 = ", format(sp$kappa[["kappa2"]], digits = 4),
        " (convergence code ", sp$convergence, ")\n", sep = "")
    cat("  log marginal likelihood: ", format(sp$log_ml, digits = 6), "\n", sep = "")
    if (!is.null(sp$var_names))
        cat("  variables: ", paste(sp$var_names, collapse = ", "), "\n", sep = "")
    invisible(x)
}

#' @rdname acpbvar-methods
#' @export
print.acpbvar_signs <- function(x, ...) {
    R_tot  <- attr(x, "R"); n_acc <- length(x)
    method <- attr(x, "method"); nt <- attr(x, "n_tries")
    sn <- attr(x, "shock_names")
    cat("acpbvar sign-identified draws: ", n_acc, " of ", R_tot, " accepted (",
        sprintf("%.1f%%", 100 * n_acc / R_tot), ")\n", sep = "")
    cat("  acceptance scheme \"", method, "\"", sep = "")
    if (!is.null(sn)) cat("; identified shocks: ", paste(sn, collapse = ", "), sep = "")
    cat("\n")
    if (method == "per_draw" && length(nt) > 0)
        cat("  rotations per accepted draw: median ", stats::median(nt), ", mean ",
            format(mean(nt), digits = 3), ", max ", max(nt), "\n  (max_tries = ",
            attr(x, "max_tries"), ")\n", sep = "")
    nf <- length(attr(x, "failed_draws"))
    if (nf > 0) cat("  ", nf, " draws dropped: no admissible rotation within max_tries\n", sep = "")
    invisible(x)
}

#' @rdname acpbvar-methods
#' @export
print.acpbvar_irf <- function(x, ...) {
    d <- dim(x); dn <- dimnames(x)
    cat("acpbvar impulse responses: ", d[1], " variables x ", d[2] - 1,
        " horizons (0..", d[2] - 1, ")\n  x ", d[3], " shocks x ", d[4],
        " accepted draws\n", sep = "")
    if (!is.null(dn[[1]])) cat("  variables: ", paste(dn[[1]], collapse = ", "), "\n", sep = "")
    if (!is.null(dn[[3]])) cat("  shocks:    ", paste(dn[[3]], collapse = ", "), "\n", sep = "")
    cat("  use summary() for quantiles by horizon and plot() for the figures\n")
    invisible(x)
}

#' @rdname acpbvar-methods
#' @export
summary.acpbvar_irf <- function(object, horizons = NULL, probs = c(0.16, 0.5, 0.84),
                                vars = NULL, shocks = NULL, ...) {
    if (!is.numeric(probs) || any(probs < 0 | probs > 1))
        stop("`probs` must be quantile levels in [0, 1].", call. = FALSE)
    d <- dim(object); dn <- dimnames(object)
    var_names   <- if (!is.null(dn[[1]])) dn[[1]] else paste0("Var ", seq_len(d[1]))
    shock_names <- if (!is.null(dn[[3]])) dn[[3]] else paste0("Shock ", seq_len(d[3]))
    vars   <- resolve_index(vars,   var_names,   "vars")
    shocks <- resolve_index(shocks, shock_names, "shocks")
    H <- d[2] - 1
    if (is.null(horizons)) horizons <- 0:H
    if (!is.numeric(horizons) || anyNA(horizons) || any(horizons != as.integer(horizons)) ||
        any(horizons < 0 | horizons > H))
        stop("`horizons` must be integers between 0 and ", H, ".", call. = FALSE)
    a <- unclass(object)
    out <- lapply(shocks, function(k) {
        res <- lapply(vars, function(v) {
            q <- apply(a[v, horizons + 1, k, , drop = FALSE], 2, stats::quantile, probs = probs)
            q <- matrix(q, nrow = length(probs),
                        dimnames = list(paste0(100 * probs, "%"), paste0("h", horizons)))
            q
        })
        names(res) <- var_names[vars]
        res
    })
    names(out) <- shock_names[shocks]
    structure(out, class = "summary.acpbvar_irf")
}

#' @rdname acpbvar-methods
#' @export
print.summary.acpbvar_irf <- function(x, digits = 3, ...) {
    for (k in names(x)) {
        for (v in names(x[[k]])) {
            cat("Response of ", v, " to ", k, "\n", sep = "")
            print(round(x[[k]][[v]], digits))
            cat("\n")
        }
    }
    invisible(x)
}

#' @rdname acpbvar-methods
#' @export
plot.acpbvar_irf <- function(x, vars = NULL, shocks = NULL, probs = c(0.16, 0.84), ...) {
    plot_irfs(x, vars = vars, shocks = shocks, probs = probs, ...)
}
