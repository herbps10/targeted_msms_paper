library(tidyverse)
library(kableExtra)

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_2")
source(glue::glue("{basepath}/env.R"))
#source(glue::glue("{basepath}/../R/manuscript.R"))

results_path <- Sys.getenv("SIMULATION_RESULTS_PATH")
if(results_path == "") stop("Please set SIMULATION_RESULTS_PATH environment variable.")

simulations <- read_rds(glue::glue("{results_path}/simulation_results.rds")) 

results_summarized <- simulations |>
  mutate(estimate = ifelse(is.nan(estimate), NA, estimate)) |>
  mutate(error = true_beta - estimate, covered = conf.low <= true_beta & conf.high >= true_beta) |>
  group_by(estimator, N, linear, term) |>
  summarize(n = n(), se = mean(std.error), me = mean(error, na.rm = TRUE), mse = mean(error^2, na.rm = TRUE), mae = mean(abs(error), na.rm = TRUE), coverage = mean(covered, na.rm = TRUE), na = mean(is.na(estimate))) |>
  ungroup()

remove_dups <- \(x) {
  x[x == lag(x)] <- ""
  x
}

tab <- results_summarized |>
  select(N, estimator, term, coverage, mae) |>
  mutate(
    term = stringr::str_sub(term, -1),
    estimator = case_when(
      estimator == "tmle" ~ "TMLE",
      estimator == "onestep" ~ "Onestep"
    )
  ) |>
  pivot_wider(names_from = "estimator", values_from = c("coverage", "mae")) |>
  arrange(N)  |>
  mutate_at(vars(starts_with("coverage")), scales::percent_format(accuracy = 0.1)) |>
  mutate_at(vars(starts_with("mae")), scales::number_format(accuracy = 0.01)) |>
  mutate_at(vars(N), remove_dups) |>
  knitr::kable(format = "latex", booktabs = TRUE, linesep = "")
