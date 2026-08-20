#!/bin/bash
#SBATCH --job-name=SimulationStudy
##SBATCH --mail-user=user@test.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --output=/log/path
#SBATCH --partition cpu_short
#SBATCH --array=1-100


module load apptainer/1.4.5

source env.sh

for var in $(compgen -v | grep '^SIMULATION_'); do
  export "APPTAINERENV_${var}=${!var}"
done

export APPTAINERENV_R_LIBS_USER=/gpfs/data/diazi07lab/R/
export APPTAINERENV_SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID"
export APPTAINERENV_TMPDIR=/gpfs/scratch/susmah01/tmp/

apptainer exec \
  --bind /gpfs/home/susmah01 \
  --bind /gpfs/data/diazi07lab \
  --bind /gpfs/scratch/susmah01 \
  /gpfs/data/diazi07lab/containers/ml_latest.sif \
  Rscript $SIMULATION_STUDY_SCRIPT_PATH
