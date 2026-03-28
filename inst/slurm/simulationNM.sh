#!/bin/bash
set -euo pipefail

SCRIPT=simNM.sbatch

submitSlurmScript() {
  # User info - These should be set before running the script
  local email=""
  local qos=""
  local partition=""
  local mailType="ALL"

  # Variables for each simulation, will be passed into the function
  local name="$1"; shift
  local cval="$1"; shift
  local nozi="$1"; shift
  local nore="$1"; shift
  local nozire="$1"; shift

  # Making the directory for output files
  mkdir "${name}${cval}"
  cd "${name}${cval}"
  cp ../"$SCRIPT" .

  # Set the time and memory allocation and low and high counts based on the 
  # number of categories
  case "$cval" in
    50) time="00:30:00"; mem="2g"; low=150; high=450 ;;
    250) time="01:00:00"; mem="4g"; low=750; high=1500 ;;
    500) time="04:30:00"; mem="7g"; low=1500; high=2500 ;;
    1000) time="10:00:00"; mem="12g"; low=2500; high=3500 ;;
    *) echo "cval not used in original simulations"; exit 1 ;;
  esac

  # Submit the job to SLURM
  sbatch --job-name="sim${name}${cval}" \
         --time="$time" \
         --mem="$mem" \
         --qos="$qos" \
         --partition="$partition" \
         --mail-type="$mailType" \
         --mail-user="$email" \
         --export=ALL,C="$cval",L="$low",H="$high",NOZI="$nozi",NORE="$nore",NOZIRE="$nozire"\
         "$SCRIPT"
  
  cd ..
}

# Submitting each simulation scenario
submitSlurmScript NMSimFunCZIDM 50 false false false
submitSlurmScript NMSimFunCDM   50 true  false false
submitSlurmScript NMSimZIDM     50 false true  false
submitSlurmScript NMSimDM       50 false false true
submitSlurmScript NMSimFunCZIDM 250 false false false
submitSlurmScript NMSimFunCZIDM 500 false false false
submitSlurmScript NMSimFunCZIDM 1000 false false false