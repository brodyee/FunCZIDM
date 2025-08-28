#!/bin/bash -l
#SBATCH --job-name=infantCaseStudy
#SBATCH --output=output_%a.out       # Output file
#SBATCH --error=error_%a.err         # Error file
#SBATCH --ntasks=1                   # Number of tasks per job (1 task per job)
#SBATCH --cpus-per-task=1            # Number of CPUs per task
#SBATCH --time=00:30:00
#SBATCH --mem=10g

module load miniforge
conda activate func

SCRIPT=$(Rscript -e "cat(system.file('scripts', 'chainsTraceplots.R', package = \
                                     'FunCZIDM'))")
Rscript "$SCRIPT"

SCRIPT=$(Rscript -e "cat(system.file('scripts', 'combineOutputs.R', package = \
                                     'FunCZIDM'))")
Rscript "$SCRIPT"