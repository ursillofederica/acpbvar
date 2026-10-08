# Internal argument checks shared by the exported functions. Every check stops
# with a message that names the argument, says what was expected and what was
# received, and omits the internal call (call. = FALSE) so that the user reads
# a single sentence.

check_Y <- function(Y, min_cols = 2L, fun = "post_sim") {
    if (!is.matrix(Y) || !is.numeric(Y))
        stop("`Y` must be a numeric matrix (T x n), one column per variable.", call. = FALSE)
    if (ncol(Y) < min_cols)
        stop("`Y` must have at least ", min_cols, " columns; got ", ncol(Y), ".", call. = FALSE)
    if (anyNA(Y) || any(!is.finite(Y)))
        stop("`Y` contains missing or non-finite values; remove or impute them before calling ",
             fun, "().", call. = FALSE)
    invisible(TRUE)
}

check_count <- function(x, name, min = 1L, max = Inf) {
    if (length(x) != 1 || !is.numeric(x) || is.na(x) || x != as.integer(x) || x < min || x > max) {
        bound <- if (is.finite(max)) paste0(min, " <= ", name, " <= ", max) else paste0(name, " >= ", min)
        stop("`", name, "` must be a single integer with ", bound, "; got ",
             paste(format(x), collapse = ", "), ".", call. = FALSE)
    }
    invisible(TRUE)
}

check_positive <- function(x, name) {
    if (length(x) != 1 || !is.numeric(x) || is.na(x) || x <= 0)
        stop("`", name, "` must be a single positive number; got ", paste(format(x), collapse = ", "), ".",
             call. = FALSE)
    invisible(TRUE)
}

check_lists <- function(y_list, X_list, p) {
    if (!is.list(y_list) || !is.list(X_list) || length(y_list) != length(X_list) || length(y_list) < 1)
        stop("`y_list` and `X_list` must be lists of the same length (one element per equation), ",
             "as returned by build_yX().", call. = FALSE)
    for (i in seq_along(y_list)) {
        if (!is.numeric(y_list[[i]]) || !is.matrix(X_list[[i]]) || nrow(X_list[[i]]) != length(y_list[[i]]))
            stop("Equation ", i, ": `X_list[[", i, "]]` must be a numeric matrix with as many rows as ",
                 "`y_list[[", i, "]]` has elements.", call. = FALSE)
    }
    invisible(TRUE)
}

check_mcmc <- function(mcmc) {
    ok <- is.list(mcmc) && all(c("samples", "spec") %in% names(mcmc)) &&
        all(c("theta", "sigma2") %in% names(mcmc$samples)) &&
        all(c("n", "p", "R") %in% names(mcmc$spec))
    if (!ok)
        stop("`mcmc` must be the object returned by post_sim().", call. = FALSE)
    invisible(TRUE)
}

check_S <- function(S, n_var) {
    if (!is.matrix(S) || !is.numeric(S) || anyNA(S) || !all(S %in% c(-1, 0, 1)))
        stop("`S` must be a numeric matrix with entries in {-1, 0, 1}: one row per identified shock, ",
             "one column per variable.", call. = FALSE)
    if (ncol(S) != n_var)
        stop("`S` has ", ncol(S), " columns but the model has ", n_var, " variables.", call. = FALSE)
    if (any(rowSums(S != 0) == 0))
        stop("Every row of `S` must contain at least one non-zero restriction; row(s) ",
             paste(which(rowSums(S != 0) == 0), collapse = ", "), " are all zero.", call. = FALSE)
    invisible(TRUE)
}

check_accepted <- function(accepted) {
    if (!is.list(accepted) || length(accepted) == 0)
        stop("`accepted` is empty: no posterior draw satisfied the sign restrictions. ",
             "Check `S` or increase `max_tries` in sign_restrict().", call. = FALSE)
    need <- c("impact", "B_rf")
    if (!all(need %in% names(accepted[[1]])))
        stop("`accepted` must be the list returned by sign_restrict().", call. = FALSE)
    invisible(TRUE)
}
