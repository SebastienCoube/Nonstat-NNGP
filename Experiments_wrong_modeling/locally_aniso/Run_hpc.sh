#!/bin/bash
#SBATCH --partition=regular
#SBATCH --job-name=R_JOB
#SBATCH --mem=30gb
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


module load R/4.2.1-foss-2022a

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

Rscript $1 >& $1.out

