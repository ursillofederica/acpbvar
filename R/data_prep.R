#' Build response and regressor matrices for the A-form VAR
#'
#' Constructs the equation-by-equation response vectors and regressor matrices
#' required by the asymmetric A-form representation of a structural VAR. For
#' equation `i`, the regressors include the lagged values of all variables, an
#' intercept, and (for `i > 1`) the contemporaneous values of variables `1`
#' through `i - 1` entered with a negative sign, so that the impact matrix is
#' lower triangular by construction.
#'
#' @param Y A `T x n` numeric matrix of observations on the `n` endogenous
#'   variables, with rows ordered chronologically.
#' @param p The VAR lag order (positive integer).
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`y_list`}{A list of length `n`; `y_list[[i]]` is the post-lag
#'       response vector for equation `i`, of length `TT = T - p`.}
#'     \item{`X_list`}{A list of length `n`; `X_list[[i]]` is the regressor
#'       matrix for equation `i`, with `TT` rows and `(i - 1) + 1 + n * p`
#'       columns (contemporaneous block, intercept, lagged variables).}
#'     \item{`n`}{Number of endogenous variables.}
#'     \item{`p`}{VAR lag order.}
#'     \item{`TT`}{Effective sample size after lag truncation.}
#'   }
#'
#' @examples
#' sim <- simulate_var_example(T = 120, seed = 1)
#' dat <- build_yX(sim$Y, p = 1L)
#'
#' length(dat$y_list)          # one response vector per equation
#' sapply(dat$X_list, ncol)    # 4, 5, 6: intercept + n*p lags, plus the
#'                             # contemporaneous variables ordered before
#'                             # the equation (entered with a minus sign)
#' @export
build_yX <- function(Y, p) {

    check_Y(Y, min_cols = 1L, fun = "build_yX")
    check_count(p, "p", min = 1L, max = nrow(Y) - 1L)

    T_full <- nrow(Y)
    n      <- ncol(Y)
    TT     <- T_full - p

    X_tilde <- matrix(1, nrow = TT, ncol = 1 + n * p)

    idx <- 2
    for (l in seq_len(p)) {
        for (j in seq_len(n)) {
            X_tilde[, idx] <- Y[(p - l + 1):(T_full - l), j]
            idx <- idx + 1
        }
    }

    y_list <- vector("list", n)
    X_list <- vector("list", n)

    for (i in seq_len(n)) {

        y_list[[i]] <- Y[(p + 1):T_full, i]

        if (i == 1) {
            X_list[[i]] <- X_tilde
        } else {
            W_i <- -Y[(p + 1):T_full, 1:(i - 1), drop = FALSE]
            X_list[[i]] <- cbind(W_i, X_tilde)
        }
    }

    list(
        y_list  = y_list,
        X_list  = X_list,
        n       = n,
        p       = p,
        TT      = TT
    )
}


#' Estimate residual variances from univariate AR fits
#'
#' For each variable, fits an autoregression (with fallback to lower orders if
#' the requested order fails) and extracts the residual variance. The resulting
#' variances enter the asymmetric Minnesota prior as scaling factors.
#'
#' @param y_list A list of `n` response vectors, as produced by [build_yX()].
#' @param ar.lags Integer; AR order attempted first. If the fit fails the
#'   function falls back to orders `2` and then `1`; if all fail, the marginal
#'   variance of the series is returned. Default is `4`.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`s2`}{Numeric vector of length `n`; residual variances.}
#'     \item{`S`}{Diagonal matrix `diag(s2)`.}
#'   }
#'
#' @examples
#' sim <- simulate_var_example(T = 120, seed = 1)
#' dat <- build_yX(sim$Y, p = 1L)
#' s2  <- build_s2(dat$y_list)$s2
#' s2                          # AR(4) residual variances, one per variable;
#'                             # they scale the Minnesota prior
#' @export
build_s2 <- function(y_list, ar.lags = 4) {
    if (!is.list(y_list) || length(y_list) < 1 || !all(vapply(y_list, is.numeric, logical(1))))
        stop("`y_list` must be a non-empty list of numeric vectors, as returned by build_yX().", call. = FALSE)
    check_count(ar.lags, "ar.lags", min = 1L)
    n <- length(y_list); s2 <- numeric(n)
    for (i in seq_len(n)){
        y_i <- y_list[[i]]
        fit <- NULL
        for (lag in c(ar.lags, 2, 1)) {
            fit <- tryCatch(
                stats::arima(y_i, order = c(lag, 0, 0), method = "ML"),
                error = function(e) NULL
            )
            if (!is.null(fit)) break
        }
        if (!is.null(fit)) {
            s2[i] <- stats::var(fit$residuals)
        } else {
            s2[i] <- stats::var(y_i)
        }
    }
    S <- diag(s2)

    list(
        s2 = s2,
        S  = S
    )
}
