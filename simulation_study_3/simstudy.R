#
# Simulation Study 3
#

library(tidyr)
library(dplyr)
library(purrr)
library(readr)

torch::install_torch()

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_3")

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
  N = c(250, 500, 1e3, 2.5e3, 5e3),
  tau = c(1:5),
  beta0 = -1,
  beta1 = 0.5,
  index = ((task_id - 1) * N_simulations + 1):(task_id * N_simulations),
)

simulations <- simulations |>
  mutate(
    seed = index
  ) |>
  mutate(
    data = pmap(list(seed, N, tau, beta0, beta1), simulate_data),
    path = glue::glue("{cache_path}/{task_id}.rds"),
    fits = pmap(list(index, N, tau, beta0, beta1, data, path), wrapper)
  )

