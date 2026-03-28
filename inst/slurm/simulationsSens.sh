#!/bin/bash
set -euo pipefail

SCRIPT=simSens.sbatch

submitSlurmScript() {
  # User info - These should be set before running the script
  local email=""
  local qos=""
  local partition=""
  local mailType="ALL"

  # Variables for each simulation, will be passed into the function
  local name="$1"; shift
  local cval="$1"; shift
  local krate="$1"; shift
  local kshape="$1"; shift
  local pshape="$1"; shift
  local prate="$1"; shift
  local talpha="$1"; shift
  local tbeta="$1"; shift
  local a="$1"; shift
  local intVar="$1"; shift

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
         --export=ALL,C="$cval",L="$low",H="$high",KRATE="$krate",KSHAPE="$kshape",PSHAPE="$pshape",PRATE="$prate",TALPHA="$talpha",TBETA="$tbeta",A="$a",INTVAR="$intVar" \
         "$SCRIPT"
  
  cd ..
}

# Submitting each simulation scenario
                # Name           C  KRATE KSHAPE PSHAPE PRATE TALPHA TBETA  A     intVar
submitSlurmScript Orig          250  900.  100.   3.0    9.0   0.01  10.0   1.0   1.0   
submitSlurmScript unifTh        250  900.  100.   3.0    9.0   1.0   1.0    1.0   1.0
submitSlurmScript midTh         250  900.  100.   3.0    9.0   1.0   5.0    1.0   1.0
submitSlurmScript weakInfoKappa 250  9.0   1.0    3.0    9.0   0.01  10.0   1.0   1.0
submitSlurmScript highShrKappa  250  900.  450.   3.0    9.0   0.01  10.0   1.0   1.0
submitSlurmScript lowShrKappa   250  900.  10.    3.0    9.0   0.01  10.0   1.0   1.0
submitSlurmScript highShrPhi    250  900.  100.   3.0    3.0   0.01  10.0   1.0   1.0 
submitSlurmScript lowShrPhi     250  900.  100.   1.0    9.0   0.01  10.0   1.0   1.0
submitSlurmScript highGlobParam 250  900.  100.   3.0    9.0   0.01  10.0   10.0  1.0
submitSlurmScript lowGlobParam  250  900.  100.   3.0    9.0   0.01  10.0   0.1   1.0
submitSlurmScript lowIntVar     250  900.  100.   3.0    9.0   0.01  10.0   1.0   0.5  
submitSlurmScript highIntVar    250  900.  100.   3.0    9.0   0.01  10.0   1.0   2.0  
submitSlurmScript vLowIntVar    250  900.  100.   3.0    9.0   0.01  10.0   1.0   0.1  
submitSlurmScript vHighIntVar   250  900.  100.   3.0    9.0   0.01  10.0   1.0   10.0 
