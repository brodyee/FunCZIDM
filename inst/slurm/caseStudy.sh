#!/bin/bash
set -euo pipefail

# User info - These should be set before running the script
email=""
qos=""
partition=""
mailType="ALL"

# Runs the script to sample from the model for the case study
JOBID=$(sbatch --parsable --qos="$qos" --partition="$partition" \
               --mail-type="$mailType" --mail-user="$email" caseStudy.sbatch)   

# Now after the sampling is done, we run the combining and traceplotting scripts
sbatch --dependency=afterok:${JOBID} --qos="$qos" --partition="$partition" \
       --mail-type="$mailType" --mail-user="$email" combineAndTraceplot.sbatch