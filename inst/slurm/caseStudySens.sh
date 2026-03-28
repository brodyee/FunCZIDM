#!/bin/bash
set -euo pipefail

SCRIPT=caseStudySens.sbatch
SCRIPT2=combineAndTraceplot.sbatch
submitSlurmScript() {
  # User info - These should be set before running the script
  local email=""
  local qos=""
  local partition=""
  local mailType="ALL"

  # Variables for each simulation, will be passed into the function
  local name="$1"; shift
  local krate="$1"; shift
  local kshape="$1"; shift
  local pshape="$1"; shift
  local prate="$1"; shift
  local talpha="$1"; shift
  local tbeta="$1"; shift
  local a="$1"; shift
  local intVar="$1"; shift

  # Submit the job to SLURM
  JOBID=$(sbatch --job-name="sim${name}" \
                 --qos="$qos" \
                 --partition="$partition" \
                 --mail-type="$mailType" \
                 --mail-user="$email" \
                 --export=ALL,NAME="$name",KRATE="$krate",KSHAPE="$kshape",PSHAPE="$pshape",PRATE="$prate",TALPHA="$talpha",TBETA="$tbeta",A="$a",INTVAR="$intVar" \
                 "$SCRIPT" | awk '{print $4}')

  # Now after the sampling is done, we run the combining and traceplotting scripts
  sbatch --dependency=afterok:${JOBID} \
         --qos="$qos" \
         --partition="$partition" \
         --mail-type="$mailType" \
         --mail-user="$email" \
         --export=ALL,NAME="$name" \
         "$SCRIPT2"

}

# Submitting each simulation scenario
                # Name             KRATE KSHAPE PSHAPE PRATE TALPHA TBETA  A     intVar
submitSlurmScript CS_Orig           900.  100.   3.0    9.0   0.01  10.0   1.0   1.0   
submitSlurmScript CS_unifTh         900.  100.   3.0    9.0   1.0   1.0    1.0   1.0
submitSlurmScript CS_midTh          900.  100.   3.0    9.0   1.0   5.0    1.0   1.0
submitSlurmScript CS_weakInfoKappa  9.0   1.0    3.0    9.0   0.01  10.0   1.0   1.0
submitSlurmScript CS_highShrKappa   900.  450.   3.0    9.0   0.01  10.0   1.0   1.0
submitSlurmScript CS_lowShrKappa    900.  10.    3.0    9.0   0.01  10.0   1.0   1.0
submitSlurmScript CS_highShrPhi     900.  100.   3.0    3.0   0.01  10.0   1.0   1.0 
submitSlurmScript CS_lowShrPhi      900.  100.   1.0    9.0   0.01  10.0   1.0   1.0
submitSlurmScript CS_highGlobParam  900.  100.   3.0    9.0   0.01  10.0   10.0  1.0
submitSlurmScript CS_lowGlobParam   900.  100.   3.0    9.0   0.01  10.0   0.1   1.0
submitSlurmScript CS_lowIntVar      900.  100.   3.0    9.0   0.01  10.0   1.0   0.5  
submitSlurmScript CS_highIntVar     900.  100.   3.0    9.0   0.01  10.0   1.0   2.0  
submitSlurmScript CS_vLowIntVar     900.  100.   3.0    9.0   0.01  10.0   1.0   0.1  
submitSlurmScript CS_vHighIntVar    900.  100.   3.0    9.0   0.01  10.0   1.0   10.0 
