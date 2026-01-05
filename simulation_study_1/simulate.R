#' Simulate data for simulation study
#'
#' @param seed random number seed
#' @param N number of observations
#'
#' @return data frame of simulated data. 
#' 
simulate_data <- function(seed = 1, N = 1e3, sigma = 0.1, linear = TRUE) {
  set.seed(seed)
  X1 <- rnorm(N, 0, 1)
  X2 <- rnorm(N, 0, 1)
  X3 <- rnorm(N, 0, 1)
  X4 <- rnorm(N, 0, 1)

  pscore <- plogis(0.5 * X1 - 0.5 * X2 + 0.2 * X3 - 0.1 * X4)

  A <- rbinom(N, 1, pscore)

  z <- 0.5 * X2 - 0.5 * X3
  if(linear == TRUE) {
    mu0 <- z
    mu1 <- z + 3 + 1.5 * X4

    Y0 <- rnorm(N, mu0, sigma)
    Y1 <- rnorm(N, mu1, sigma)
  }
  else {
    mu0 <- plogis(z)
    mu1 <- plogis(z + 0.5 + X4)

    Y0 <- rbinom(N, 1, mu0)
    Y1 <- rbinom(N, 1, mu1)
  }
  Y <- ifelse(A == 1, Y1, Y0)
  
  data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, A = A, Y = Y, mu0 = mu0, mu1 = mu1, Y0 = Y0, Y1 = Y1, pscore = pscore)
}
