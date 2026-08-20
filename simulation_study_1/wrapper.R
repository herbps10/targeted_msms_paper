library(automsm)
library(broom)
library(torch)

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

wrapper <- function(index, N, linear, scenario, data, path) {
  # Check cache
  if(!file.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  cached_res <- NULL
  if(file.exists(path)) {
    cached_res <- read_rds(path)
    if(sum(cached_res$N == N & cached_res$linear == linear & cached_res$scenario == scenario) > 0) {
      print("Cache hit")
      return()
    }
    else {
      print("No cache")
    }
  }

  print(glue::glue("Starting: {index} {N} {linear} {scenario}"))

  largedat <- simulate_data(seed = 1, N = 1e5, linear = linear)
  true_params <- tibble(
    term = c("(Intercept)", "X4"),
    true_beta = coef(lm(mu1 - mu0 ~ 1 + X4, data = largedat))
  )

  learners_trt <- learners_outcome <- "SL.random"

  if(scenario %in% c(1, 3)) {
    learners_outcome <- c("SL.glm.interaction")
  }
  if(scenario %in% c(1, 2)) {
    learners_trt <- c("SL.glm.interaction")
  }

  set.seed(10016)
  start_time <- Sys.time()
  fit <- automsm::cate(
    data, 
    c("X1", "X2", "X3", "X4"), "A", "Y", 
    formula = ~ 1 + X4,
    loss = loss_squared_error,
    working_model = working_model_linear,
    outcome_type = "binomial",
    tmle = tmle_control(
      criterion = "epsilon",
      fluctuation = "logistic"
    ),
    bayes = bayes_control(
      chains = 4,
      warmup = 500,
      draws = 500,
      prior = \(beta) sum(dnorm(as.numeric(beta), mean = 0, sd = 5, log = TRUE)),
      eps_max = 100,
      fluctuation = "logistic"
    ),
    nuisance = nuisance_control(
      learners_outcome = learners_outcome, 
      learners_trt = learners_trt,
      epsilon = 1e-5,
      outer_folds = 1
    )
  )
  total_time <- Sys.time() - start_time

  res <- tidy(fit) |>
    mutate(scenario = scenario, simulation_index = index, N = N, linear = linear, path = path, total_time = total_time) |>
    left_join(true_params) |>
    left_join(tibble(estimator = "bayes", rejected = mean(fit$bayes_tmle$rejected)))

  if(is.null(cached_res)) {
    write_rds(res, path)
  }
  else {
    write_rds(bind_rows(cached_res, res), path)
  }

  return()
}


