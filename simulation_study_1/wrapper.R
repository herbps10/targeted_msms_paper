library(TargetedMSM)
library(broom)
library(torch)

wrapper <- function(index, N, linear, scenario, data, path) {
  # Check cache
  if(!file.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  cached_res <- NULL
  if(file.exists(path)) {
    cached_res <- read_rds(path)
    if(sum(cached_res$index == index & cached_res$N == N & cached_res$linear == linear & cached_res$scenario == scenario) > 0) return()
  }

  print(glue::glue("Starting: {index} {N} {linear} {scenario}"))

  largedat <- simulate_data(seed = 1, N = 1e4, linear = linear)
  true_params <- tibble(
    term = c("(Intercept)", "X4"),
    true_beta = coef(lm(mu1 - mu0 ~ 1 + X4, data = data))
  )

  learners_trt <- learners_outcome <- "SL.random"

  if(scenario %in% c(1, 3)) {
    learners_outcome <- c("SL.glm.interaction")
  }
  if(scenario %in% c(1, 2)) {
    learners_trt <- c("SL.glm.interaction")
  }

  fit <- TargetedMSM::treatment_effect_modification(
    data, 
    c("X1", "X2", "X3", "X4"), "A", "Y", 
    formula = ~ 1 + X4,
    learners_outcome = learners_outcome, 
    learners_trt = learners_trt,
    loss = loss_squared_error,
    working_model = working_model_linear,
    outcome_type = "binomial",
    tmle = TRUE,
    tmle_linear = FALSE,
    bayes = TRUE,
    bayes_draws = 5e3,
    bayes_prior = \(beta) sum(dnorm(as.numeric(beta), mean = 0, sd = 5, log = TRUE)),
    epsilon = 1e-5
  )

  res <- tidy(fit) |>
    mutate(scenario = scenario, simulation_index = index, N = N, linear = linear, path = path) |>
    left_join(true_params)

  if(is.null(cached_res)) {
    write_rds(res, path)
  }
  else {
    write_rds(bind_rows(cached_res, res), path)
  }

  return()
}

SL.random <- function (Y, X, newX, family, obsWeights, id, ...) 
{
  meanY <- weighted.mean(Y, w = obsWeights)
  pred <- runif(nrow(newX), 0.1, 0.9)
  fit <- list(object = meanY)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.random")
  return(out)
}

predict.SL.random <- function (object, newdata, family, X = NULL, Y = NULL, ...) 
{
    runif(nrow(newdata), 0.1, 0.9)
}
