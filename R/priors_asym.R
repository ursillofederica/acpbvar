#' Build the asymmetric Minnesota prior in A-form
#'
#' Constructs the per-equation prior means and variances for the structural
#' coefficients of an asymmetric Minnesota prior, as required by the A-form
#' representation. Own-lag, cross-lag, and intercept components have separate
#' tightness hyperparameters; lag decay is quadratic in the lag order.
#'
#' @param y_list,p,s2 As produced by [build_yX()] and [build_s2()].
#' @param kappa1 Own-lag tightness (positive scalar). Smaller values pull
#'   own-lag coefficients toward zero (or toward one, under `unit_root_mean`).
#' @param kappa2 Cross-lag tightness (positive scalar), typically much smaller
#'   than `kappa1`.
#' @param kappa3 Intercept variance (positive scalar). Default `100` is the
#'   diffuse choice.
#' @param unit_root_mean Logical; if `TRUE`, sets the prior mean of the
#'   coefficient on the first own lag to `1` (random-walk prior); otherwise
#'   zero. Default `FALSE`.
#'
#' @return A list with per-equation prior components:
#'   `m_alpha`, `V_alpha` (contemporaneous coefficients),
#'   `m_beta`, `V_beta` (intercept and lagged coefficients),
#'   `m`, `V` (concatenations used by the conjugate sampler),
#'   `nu`, `S` (Inverse-Gamma hyperparameters for the structural variances).
#'
#' @keywords internal
build_prior_asym <- function(y_list, p, s2, kappa1, kappa2, kappa3 = 100, unit_root_mean = FALSE) {

    n <- length(y_list)
    p <- p

    m_alpha_list <- vector("list", n)
    V_alpha_list <- vector("list", n)
    m_beta_list  <- vector("list", n)
    V_beta_list  <- vector("list", n)
    m_list       <- vector("list", n)
    V_list       <- vector("list", n)

    nu <- numeric(n)
    S  <- numeric(n)

    for (i in seq_len(n)) {

        m_alpha_list[[i]] <- numeric(i - 1)
        nu[i]             <- 1 + 0.5 * i
        S[i]              <- 0.5 * s2[i]

        if (i == 1) {
            V_alpha_list[[i]] <- numeric(0)
        } else {
            V_alpha_list[[i]] <- 1 / s2[1:(i - 1)]
        }
        m_beta <- numeric(1 + n * p)

        if (unit_root_mean) {
            m_beta[1 + i] <- 1
        }

        m_beta_list[[i]] <- m_beta

        v_beta <- numeric(1 + n * p)
        v_beta[1] <- kappa3

        idx <- 2
        for (l in seq_len(p)) {
            for (j in seq_len(n)) {
                if (j == i) {
                    v_beta[idx] <- kappa1 / (l^2 * s2[i])
                } else {
                    v_beta[idx] <- kappa2 / (l^2 * s2[j])
                }
                idx <- idx + 1
            }
        }

        V_beta_list[[i]] <- v_beta

        m_list[[i]] <- c(m_alpha_list[[i]], m_beta_list[[i]])
        V_list[[i]] <- c(V_alpha_list[[i]], V_beta_list[[i]])
    }

    list(
        m_alpha = m_alpha_list,
        V_alpha = V_alpha_list,
        m_beta  = m_beta_list,
        V_beta  = V_beta_list,
        m       = m_list,
        V       = V_list,
        nu      = nu,
        S       = S
    )
}
