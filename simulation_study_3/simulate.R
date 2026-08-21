#' Simulate data for simulation study
#'
#' @param seed random number seed
#' @param N number of observations
#' @param tau number of observations
#'
#' @return data frame of simulated data. 
#' 
simulate_data <- function(seed = 1, N = 1e3, tau, beta0, beta1) {
  automsm::simulate_longitudinal_dose_response(N, tau, beta0 = beta0, beta1 = beta1, seed = seed, binary = TRUE)
}

compute_true_msm <- function(
  tau = 3,
  beta0 = -1,
  beta1 = 0.5,
  N_inner = 2e4,
  seed = 10016
) {
  set.seed(seed)

  # Enumerate all regimes in {0, 1}^tau
  regimes <- expand.grid(rep(list(0:1), tau))
  colnames(regimes) <- paste0("a", 1:tau)
  regimes$v <- rowSums(regimes)
  regimes$psi <- NA

  # For each regime a_bar, compute psi_P^{(\bar{a})} = E[Y_{\bar{a}}]
  # by MC integration over the counterfactual covariate process.
  for(r in seq_len(nrow(regimes))) {
    a_bar <- as.numeric(regimes[r, 1:tau])
    v_a <- regimes$v[r]

    L_cf <- matrix(NA, nrow = N_inner, ncol = tau)
    L_cf[, 1] <- rnorm(N_inner)
    if(tau >= 2) {
      for(t in 2:tau) {
        L_cf[, t] <- rnorm(
          N_inner,
          mean = 0.4 * L_cf[, t - 1] + 0.5 * a_bar[t - 1]
        )
      }
    }

    # Compute E[Y | L_tau, A_bar = a_bar] analytically
    p_Y <- plogis(beta0 + beta1 * v_a + 0.2 * L_cf[, tau])
    regimes$psi[r] <- mean(p_Y)
  }

  fit <- glm(psi ~ 1 + v, data = regimes, family = binomial)
  fit_lsq <- glm(psi ~ 1 + v, data = regimes)

  list(
    beta_true = coef(fit),
    beta_lsq = coef(fit_lsq),
    regime_table = regimes,
    n_inner = N_inner
  )
}
