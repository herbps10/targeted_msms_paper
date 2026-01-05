library(tidyverse)
library(kableExtra)

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_1")
source(glue::glue("{basepath}/env.R"))
#source(glue::glue("{basepath}/../R/manuscript.R"))

results_path <- Sys.getenv("SIMULATION_RESULTS_PATH")
if(results_path == "") stop("Please set SIMULATION_RESULTS_PATH environment variable.")

simulations <- read_rds(glue::glue("{results_path}/simulation_results.rds")) 

results_summarized <- simulations |>
  mutate(error = true_beta - estimate, covered = conf.low <= true_beta & conf.high >= true_beta) |>
  group_by(scenario, estimator, N, linear, term) |>
  summarize(n = n(), se = mean(std.error), me = mean(error, na.rm = TRUE), mse = mean(error^2, na.rm = TRUE), mae = mean(abs(error), na.rm = TRUE), coverage = mean(covered, na.rm = TRUE), na = mean(is.na(estimate)))
