library(automsm)
library(broom)
library(torch)

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
    epsilon = 1e-3
  )

  fit <- automsm::dose_response(
    data, 
    c("x1", "x2", "x3"), "a", "y", 
    formula = ~-1 + splines::bs(a, knots = seq(1, treatments, length.out = n_knots)[-c(n_knots)], Boundary.knots = c(1, treatments)),
    nuisance = nuisance,
    outcome_type = "continuous",
    tmle = TRUE
  )

  res <- tidy(fit) |>
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
