#' Simulate a small structural VAR with known sign restrictions
#'
#' Generates a three-variable structural VAR(1) calibrated to deliver impulse
#' responses with the canonical macroeconomic sign pattern: a positive supply
#' shock raises output and lowers inflation; a positive demand shock raises
#' both output and inflation; a positive monetary-policy shock lowers output
#' and inflation while raising the interest rate. The function returns the
#' simulated sample together with the true structural impact matrix and the
#' true reduced-form coefficients, so that recovery can be validated against
#' the data-generating process.
#'
#' The simulation discards a burn-in of 100 observations to wash out the
#' influence of the zero initial condition.
#'
#' @param T Integer; number of post-burn-in observations to return.
#'   Default `200`.
#' @param seed Optional integer seed passed to [set.seed()]. Default `NULL`
#'   (generator left untouched).
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`Y`}{`T x 3` matrix of simulated observations on output,
#'       inflation, and the interest rate.}
#'     \item{`impact_true`}{The `3 x 3` structural impact matrix used in the
#'       DGP; columns index shocks (supply, demand, monetary).}
#'     \item{`B_true`}{The `3 x 3` reduced-form coefficient matrix on the
#'       first lag.}
#'     \item{`intercept_true`}{The `3`-vector of intercepts in the DGP.}
#'     \item{`var_names`}{Variable names.}
#'     \item{`shock_names`}{Shock names corresponding to the columns of
#'       `impact_true`.}
#'   }
#'
#' @examples
#' sim <- simulate_var_example(T = 200, seed = 42)
#' head(sim$Y)
#' sim$impact_true
#'
#' @export
simulate_var_example <- function(T = 200, seed = NULL) {
    check_count(T, "T", min = 10L)

    if (!is.null(seed)) set.seed(seed)
    n <- 3

    # True structural impact matrix (rows = variables, cols = shocks)
    #          AS     AD     MP
    impact_true <- matrix(c(
         0.50,  0.40, -0.30,   # output
        -0.30,  0.50, -0.20,   # inflation
         0.10,  0.20,  0.40    # interest rate
    ), nrow = 3, byrow = TRUE)

    # True reduced-form persistence on the first lag
    B_true <- matrix(c(
         0.70,  0.10,  0.05,
         0.10,  0.60,  0.05,
         0.05,  0.10,  0.80
    ), nrow = 3, byrow = TRUE)

    intercept_true <- c(0.05, 0.02, 0.03)

    burn    <- 100L
    T_total <- T + burn

    Y_full <- matrix(0, T_total, n)
    eps    <- matrix(stats::rnorm(T_total * n), T_total, n)

    for (t in 2:T_total) {
        Y_full[t, ] <- intercept_true +
                       B_true %*% Y_full[t - 1, ] +
                       impact_true %*% eps[t, ]
    }

    Y <- Y_full[(burn + 1):T_total, , drop = FALSE]
    colnames(Y) <- c("output", "inflation", "interest_rate")

    list(
        Y              = Y,
        impact_true    = impact_true,
        B_true         = B_true,
        intercept_true = intercept_true,
        var_names      = colnames(Y),
        shock_names    = c("supply", "demand", "monetary")
    )
}
