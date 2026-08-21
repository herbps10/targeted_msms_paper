library(dplyr)
library(readr)
library(tidyr)
#library(kableExtra)

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_3")
source(glue::glue("{basepath}/env.R"))
#source(glue::glue("{basepath}/../R/manuscript.R"))

results_path <- Sys.getenv("SIMULATION_RESULTS_PATH")
if(results_path == "") stop("Please set SIMULATION_RESULTS_PATH environment variable.")

simulations <- read_rds(glue::glue("{results_path}/simulation_results.rds")) 

results_summarized <- simulations |>
  mutate(estimate = ifelse(is.nan(estimate), NA, estimate)) |>
  mutate(error = true_values - estimate, covered = conf.low <= true_values & conf.high >= true_values) |>
  group_by(estimator, N, tau, beta0, beta1, term, msm_loss) |>
  summarize(n = n(), observed_sd = sd(estimate), beta_mean = mean(estimate), se = mean(std.error), me = mean(error, na.rm = TRUE), mse = mean(error^2, na.rm = TRUE), mae = mean(abs(error), na.rm = TRUE), coverage = mean(covered, na.rm = TRUE), na = mean(is.na(estimate)), mean_time = mean(time)) |>
  ungroup()

write_rds(results_summarized, glue::glue("{results_path}/results_summarized.rds"))

remove_dups <- \(x) {
  x[x == lag(x)] <- ""
  x
}

tab <- results_summarized |>
  select(msm_loss, tau, N, estimator, term, coverage, mae) |>
  filter(estimator %in% c("ltmle", "tmle", "onestep")) |>
  mutate(
    estimator = case_when(
      estimator == "ltmle" ~ "LTMLE",
      estimator == "tmle" ~ "TMLE",
      estimator == "onestep" ~ "Onestep"
    )
  ) |>
  pivot_wider(names_from = "estimator", values_from = c("coverage", "mae")) |>
  arrange(term, tau, N)  |>
  mutate_at(vars(starts_with("coverage")), scales::percent_format(accuracy = 0.1)) |>
  mutate_at(vars(starts_with("mae")), scales::number_format(accuracy = 0.01)) |>
  mutate_at(vars(tau, N), remove_dups) |>
  knitr::kable(format = "latex", booktabs = TRUE, linesep = "")
