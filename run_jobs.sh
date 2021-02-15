#!/bin/bash
#SBATCH --array=1-100
#SBATCH --time=05:59:00           # time (HH:MM:SS)
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --error=err/%j.err
#SBATCH --output=err/%j.out
#SBATCH --mail-user=asad.haris5862@gmail.com # Send email updates to you or someone else
#SBATCH --mail-type=BEGIN         # send an email in all cases (job started, job ended, job aborted)


## ARGUMENTS ARE:
## 1. seed
## 2. n: Sample Size $1
## 3. sigma: SD of residual
## 4. sigma_e: SD of error in data

module load StdEnv/2020 r/4.0.2

export R_LIBS=~/local/R_libs/
R CMD BATCH --no-save --no-restore "--args $SLURM_ARRAY_TASK_ID $1 1 0.5" simulations_test1.R bla.Rout

#### SBATCH --account=def-rwplatt   # replace this with your own account
#### SBATCH --ntasks=4              # number of processes
#### SBATCH --mem-per-cpu=2048M      # memory; default unit is megabytes
#### SBATCH --time=0-01:15           # time (DD-HH:MM)

