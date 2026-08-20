#
# Simulation Study
#

library(tidyr)
library(dplyr)
library(purrr)
library(readr)

torch::install_torch()

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_1")

# Load simulation files
source(glue::glue("{basepath}/simulate.R"))
source(glue::glue("{basepath}/wrapper.R"))
source(glue::glue("{basepath}/env.R"))

cache_path <- Sys.getenv("SIMULATION_CACHE_PATH")
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if(cache_path == "") stop("Please set SIMULATION_CACHE_PATH environment variable.")
if(task_id == "") stop("Task id not set. Please set SLURM_ARRAY_TASK_ID, or run simulations through a Slurm job array, which will set this environment variable for you.")

task_id <- as.numeric(task_id)

N_simulations <- 5
simulations <- expand_grid(
  index = ((task_id - 1) * N_simulations + 1):(task_id * N_simulations),
  scenario = 1:4,
  N = c(100, 250, 5e2, 750, 1e3),
  linear = FALSE
)

simulations <- simulations |>
  mutate(
    seed = index
  ) |>
  mutate(
    data = pmap(list(seed, N, linear), simulate_data, sigma = 0.1),
    path = glue::glue("{cache_path}/{index}.rds"),
    fits = pmap(list(index, N, linear, scenario, data, path), wrapper)
  )
