#!/bin/bash -l
#SBATCH --output=output_%a.out      # Output file
#SBATCH --error=error_%a.err        # Error file
#SBATCH --ntasks=1                  # Number of tasks per job (1 task per job)
#SBATCH --cpus-per-task=1           # Number of CPUs per task
#SBATCH --array=0-99%100            # Array job specification (100 tasks total, 100 running concurrently)

module load miniforge
conda activate func

# List of seeds (You can adjust these as needed)
SEEDS=(85744 39350 62007 27231 17228 23170 28217 21995 92637 33582 82945 97198 \
       21490 65985 13704 51186 92608 19003 57905 71554 13601 22752 83194 63649 \
       52799 33528 18501 99162 93396 12349 39975 94940 66026 99507 87312 10869 \
       76079 87705 42204 82733 24992 60546 32898 87808 68099 82185 49493 33261 \
       59667 54386 69836 36964 62281 78712 61596 45089 66142 92100 95047 18342 \
       97134 76054 59190 16956 52012 97448 72537 74642 64659 74339 89717 90645 \
       55079 84367 59699 46207 20271 17741 29556 59305 83988 97163 44867 64963 \
       36372 34887 66534 27667 11008 60419 60361 72188 53605 96824 33170 23939 \
       69498 71559 15470 58697)

# Get the seed for the current task
SEED=${SEEDS[$SLURM_ARRAY_TASK_ID]}

# Run your R script with the seed as an argument
SCRIPT=$(Rscript -e "cat(system.file('scripts', 'simulationStudy.R', package = \
                                     'FunCZIDM'))")
Rscript "$SCRIPT" --seed $SEED -n 50 -c $C -p 10 -i 10000 -b 5000 -l $L \
                  -h $H --propSigma .1 --thin 5 --kappaShape 100 \
                  --kappaRate 900 --noZI "$NOZI" --noRE "$NORE" \
                  --noZIRE "$NOZIRE"