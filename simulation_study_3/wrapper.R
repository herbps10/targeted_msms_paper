library(automsm)
library(broom)
library(torch)

wrapper <- function(index, N, tau, beta0, beta1, data, path) {
  # Check cache
  print(glue::glue("Starting: {index} {N} {tau}"))

  if(!file.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  cached_res <- NULL
  if(file.exists(path)) {
    cached_res <- read_rds(path)
    if(sum(cached_res$simulation_index == index & cached_res$N == N & cached_res$tau == tau & cached_res$beta0 == beta0 & cached_res$beta1 == beta1) > 0) {
      print("Cache hit")
      return()
    }
    else {
      print("No cache")
    }
  }


  true_msm <- compute_true_msm(
    tau = tau,
    beta0 = beta0,
    beta1 = beta1,
    N_inner = 1e4,
    seed = 10016
  )

  #learners_trt <- c("SL.glm", "SL.glm.interaction")
  #learners_outcome <- c("SL.glm", "SL.glm.interaction")
  learners_trt <- c("SL.glm")
  learners_outcome <- c("SL.glm")

  # ----- automsm -----
  Ls <- lapply(1:tau, function(t) paste0("L", t))
  As <- paste0("A", 1:tau)
  regimes <- expand.grid(rep(list(c(0, 1)), tau))
  colnames(regimes) <- As

  start_time <- Sys.time()
  fit <- automsm::longitudinal_dose_response(
    data, 
    Ls,
    As,
    "Y",
    regimes = regimes,
    summary_measures = function(regimes) data.frame(v = rowSums(regimes)),
    formula = ~1 + v,
    loss = loss_cross_entropy_logit,
    outcome = "binomial",
    tmle = tmle_control(eif_tol = 0.05),
    nuisance = nuisance_control(
      learners_outcome = learners_outcome, 
      learners_trt = learners_trt
    )
  )
  automsm_time <- Sys.time() - start_time

  automsm_res <- tidy(fit) |>
    mutate(
      simulation_index = index, 
      term = ifelse(term == "(Intercept)", term, "v"),
      N = N, 
      tau = tau,
      beta0 = beta0,
      beta1 = beta1,
      path = path,
      true_values = rep(true_msm$beta_true, 2),
      time = automsm_time
    ) |>
    left_join(tibble(
      estimator = "tmle",
      term = c("(Intercept)", "v"),
      eif_mean = colMeans(fit$tmle$eif),
      iter = fit$tmle$iter,
      converged = fit$tmle$converged
    )) |>
    left_join(tibble(
      estimator = "onestep",
      term = c("(Intercept)", "v"),
      eif_mean = colMeans(fit$onestep$eif),
    )) |>
    bind_rows(tibble(
      simulation_index = index,
      term = c("(Intercept)", "v"),
      N = N,
      tau = tau,
      beta0 = beta0,
      beta1 = beta1,
      path = path,
      true_values = true_msm$beta_true,
      time = automsm_time,
      estimator = "plugin",
      estimate = fit$plugin$est
    ))

  # ----- ltmle -----
  regime_mat <- as.matrix(expand.grid(rep(list(c(0, 1)), tau)))
  colnames(regime_mat) <- As
  num_regimes <- nrow(regime_mat)

  regimes_arr <- array(0, dim = c(nrow(data), tau, num_regimes))
  for(j in 1:num_regimes) {
    regimes_arr[, , j] <- matrix(regime_mat[j, ], nrow = nrow(data), ncol = tau, byrow = TRUE)
  }

  v <- rowSums(regime_mat)
  summary.measures <- array(v, dim = c(num_regimes, 1, 1))
  dimnames(summary.measures) <- list(NULL, "v", NULL)

  start_time <- Sys.time()
  ltmle_fit <- ltmle::ltmleMSM(
    data = data[, c(rbind(unlist(Ls), As), "Y")],
    Anodes = As,
    Lnodes = unlist(Ls),
    Ynodes = "Y",
    survivalOutcome = FALSE,
    regimes = regimes_arr,
    working.msm = "Y ~ 1 + v",
    summary.measures = summary.measures,
    final.Ynodes = "Y",
    msm.weights = NULL,
    SL.library = list(Q = learners_outcome, g = learners_trt),
    estimate.time = FALSE,
    SL.cvControl = list(V = 5)
  )
  ltmle_time <- Sys.time() - start_time

  ltmle_summary <- summary(ltmle_fit)
  ltmle_coefs <- ltmle_summary$cmat

  ltmle_res <- tibble::tibble(
    simulation_index = index,
    term = rownames(ltmle_coefs),
    estimate = ltmle_coefs[, "Estimate"],
    std.error = ltmle_coefs[, "Std. Error"],
    conf.low = ltmle_coefs[, "CI 2.5%"],
    conf.high = ltmle_coefs[, "CI 97.5%"],
    estimator = "ltmle",
    N = N,
    tau = tau,
    beta0 = beta0,
    beta1 = beta1, 
    path = path,
    true_values = true_msm$beta_true,
    time = ltmle_time
  )

  res <- bind_rows(automsm_res, ltmle_res)

  if(is.null(cached_res)) {
    write_rds(res, path)
  }
  else {
    write_rds(bind_rows(cached_res, res), path)
  }

  return()
}
