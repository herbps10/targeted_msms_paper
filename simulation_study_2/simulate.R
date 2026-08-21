#' Simulate data for simulation study
#'
#' @param seed random number seed
#' @param N number of observations
#'
#' @return data frame of simulated data. 
#' 
simulate_data <- function(seed = 1, N = 1e3, treatments = 25, sigma = 0.1, linear = TRUE) {
  set.seed(seed)
  x1 <- runif(N, -1, 1)
  x2 <- runif(N, -1, 1)
  x3 <- runif(N, -1, 1)
  
  a <- numeric(N)
  for(i in 1:N) {
    prob <- rep(1, treatments)
    if(x1[i] < 0) {
      prob[1:(floor(treatments / 2))] <- 3
    }
    a[i] <- sample(1:treatments, 1, replace = TRUE, prob = prob)
  }
  
  mu_a <- matrix(ncol = treatments, nrow = N)
  colnames(mu_a) <- paste0("mu_", 1:treatments)
  y_a  <- matrix(ncol = treatments, nrow = N)
  colnames(y_a) <- paste0("y_", 1:treatments)
  for(treatment in 1:treatments) {
    mu_a[, treatment] <- x1 + 2 / treatment + 1 / (treatments - treatment + 1)
    y_a[, treatment]  <- rnorm(N, mu_a[, treatment], sigma)
  }
  
  mu <- map_dbl(1:N, \(index) mu_a[index, a[index]])
  y  <- map_dbl(1:N, \(index) y_a[index, a[index]])

  tibble(x1, x2, x3, a, y)
}
