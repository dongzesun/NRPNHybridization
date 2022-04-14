#!/bin/bash -
#SBATCH -J Hybrid                   # Job Name
#SBATCH -o out           # Output file name
#SBATCH -e err               # Error file name
#SBATCH -n 24                  # Number of cores
#SBATCH --ntasks-per-node 24        # number of MPI ranks per node
#SBATCH -t 24:0:00   # Run time
#SBATCH -A sxs                # Account name
#SBATCH --no-requeue

# DO NOT MANUALLY EDIT THIS FILE! See `MakeSubmit.py update -h`.
# This is for submitting a batch job on 'wheeler'.
umask 0022
. bin/this_machine.env || echo 'Not using env file.'
set -x
export PATH=$(pwd -P)/bin:$PATH

python Hybridization.py
