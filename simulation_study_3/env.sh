#/bin/bash

# Path to main simulation study R script
export SIMULATION_STUDY_SCRIPT_PATH="/gpfs/home/susmah01/targeted_msms_paper/simulation_study_3/simstudy.R"

# Path to cache folder for temporary results
export SIMULATION_CACHE_PATH="/gpfs/scratch/susmah01/targeted_msms_paper/cache_3"

# Path to log folder
export SIMULATION_LOG_PATH="/gpfs/scratch/susmah01/targeted_msms_paper/logs_3"

# Path to results folder for final results
export SIMULATION_RESULTS_PATH="/gpfs/home/susmah01/targeted_msms_paper/simulation_study_3/results"

# Slurm job configuration

export JOB_NAME="msm3" 

export PARTITION=cpu_short # SLURM partition

export MAIL_USER=herbert.susmann@nyulangone.org

export CPUS_PER_TASK=1 # CPU to assign to each task in the array

export MEM_PER_CPU=2GB # Memory assigned per CPU

export TIME=12:00:00 # Job time limit

export ARRAY=1-100 # Array indexes to run
