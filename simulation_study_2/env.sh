#/bin/bash

# Path to main simulation study R script
export SIMULATION_STUDY_SCRIPT_PATH="/gpfs/home/$USER/targeted_msms_paper/simulation_study_2/simstudy.R"

# Path to cache folder for temporary results
export SIMULATION_CACHE_PATH="/gpfs/scratch/$USER/targeted_msms_paper/cache_two"

# Path to log folder
export SIMULATION_LOG_PATH="/gpfs/scratch/$USER/targeted_msms_paper/logs_two"

# Path to results folder for final results
export SIMULATION_RESULTS_PATH="/gpfs/home/$USER/targeted_msms_paper/simulation_study_2/results"

# Slurm job configuration

export JOB_NAME="msm2" 

export PARTITION=cpu_medium # SLURM partition

export MAIL_USER=herbert.susmann@nyulangone.org

export CPUS_PER_TASK=1 # CPU to assign to each task in the array

export MEM_PER_CPU=2GB # Memory assigned per CPU

export TIME=24:00:00 # Job time limit

export ARRAY=1-250 # Array indexes to run
