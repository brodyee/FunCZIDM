#!/bin/bash
#SBATCH --job-name=infantCaseStudy
#SBATCH --output=output_%a.out       # Output file
#SBATCH --error=error_%a.err         # Error file
#SBATCH --ntasks=1                   # Number of tasks per job (1 task per job)
#SBATCH --cpus-per-task=1            # Number of CPUs per task
#SBATCH --array=0-3%4 
#SBATCH --time=02:30:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --mem=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brody.erlandson@colostate.edu

module load miniforge
conda activate func

# List of seeds
SEEDS=(18342 33261 87808 67235)
# Get the seed for the current task
SEED=${SEEDS[$SLURM_ARRAY_TASK_ID]}

Rscript caseStudy.R -i 85000 -b 75000 --seed $SEED --toReturn RA --outputFolder "caseStudy" --thin 10 --df 4
