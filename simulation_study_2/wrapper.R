library(automsm)
library(broom)
library(torch)

oracle_sim2 <- function(data, K = 25, sigma = 0.1) {
  n <- nrow(data)
  pi_a <- matrix(1/K, n, K)
  low <- data$x1 < 0
  pi_a[low, 1:12]  <- 3/49
  pi_a[low, 13:K]  <- 1/49
  k <- seq_len(K)
  mu_a <- outer(data$x1, rep(1, K)) +
          matrix(2/k + 1/(K - k + 1), n, K, byrow = TRUE)

  idx <- cbind(seq_len(n), match(data$a, sort(unique(data$a))))

  list(
    pi = pi_a[idx], 
    pi_a = pi_a, 
    mu = mu_a[idx], 
    mu_a = mu_a,
    condvar = rep(0, n),
    condvar_a = matrix(0, n, K)
  )
}

wrapper <- function(index, N, treatments, sigma, linear, data, path) {
  # Check cache
  if(!file.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  cached_res <- NULL
  if(file.exists(path)) {
    cached_res <- read_rds(path)
    if(sum(cached_res$N == N & cached_res$linear == linear & cached_res$treatments == treatments & cached_res$sigma == sigma) > 0) {
      print("Cache hit")
      return()
    }
    else {
      print("No cache")
    }
  }

  print(glue::glue("Starting: {index} {N} {linear} {path}"))

  n_knots <- 5
  true_beta <- lm(mu ~ -1 + splines::bs(a, knots = seq(1, treatments, length.out = n_knots)[-c(n_knots)], Boundary.knots = c(1, treatments)), data = tibble(a = 1:treatments, mu = 2 / a + 1 / (treatments - a + 1)))
  true_params <- tibble(
    term = names(coef(true_beta)),
    true_beta = coef(true_beta)
  )

  learners_trt <- c("SL.ranger", "SL.glm", "SL.glm.interaction", "SL.earth")
  learners_outcome <- c("SL.ranger", "SL.glm", "SL.glm.interaction", "SL.earth")

  nuisance <- automsm::nuisance_control(
    learners_trt = learners_trt,
    learners_outcome = learners_outcome,
    outer_folds = 5,
    epsilon = 1e-3
  )

  fit <- automsm::dose_response(
    data, 
    c("x1", "x2", "x3"), "a", "y", 
    formula = ~-1 + splines::bs(a, knots = seq(1, treatments, length.out = n_knots)[-c(n_knots)], Boundary.knots = c(1, treatments)),
    nuisance = nuisance,
    outcome_type = "continuous",
    tmle = tmle_control(criterion = "epsilon")
  )

  fit_oracle <- automsm::dose_response(
    data, 
    c("x1", "x2", "x3"), "a", "y", 
    formula = ~-1 + splines::bs(a, knots = seq(1, treatments, length.out = n_knots)[-c(n_knots)], Boundary.knots = c(1, treatments)),
    nuisance_estimates = oracle_sim2(data, treatments, sigma),
    outcome_type = "continuous",
    tmle = tmle_control(criterion = "epsilon")
  )

  res <- bind_rows(
    tidy(fit) |> mutate(nuisance = "estimated") |>
      left_join(tibble(estimator = "tmle", solved = fit$tmle$solved)),
    tidy(fit_oracle) |> mutate(nuisance = "oracle") |> 
      left_join(tibble(estimator = "tmle", solved = fit_oracle$tmle$solved))
  ) |>
    mutate(
      simulation_index = index, 
      N = N, 
      treatments = treatments,
      sigma = sigma,
      linear = linear, 
      path = path
    ) |>
    left_join(true_params)

  if(is.null(cached_res)) {
    write_rds(res, path)
  }
  else {
    write_rds(bind_rows(cached_res, res), path)
  }

  return()
}
