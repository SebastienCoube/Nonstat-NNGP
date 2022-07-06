#!/bin/bash
#
#### R JOB ######################
#
#  Usage: sbatch sR.sh myscript.R 
#
#################################
#
#
#
########################################### BEGIN SLURM OPTIONS  ##################################
#
#
# JOB AND OUTPUT FILE NAMES. Uncomment if needed.
#   If not specified, the default output file name will be slurm-NNNN.out, 
#    where NNNN is the Slurm job id.
#SBATCH --job-name=NONSTAT_NNGP
#SBATCH --output=acceptance_output_delayed.out
#
# MAIL NOTIFICATIONS (at job start, step, end or failure). Uncomment if needed.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastien.coube@univ-pau.fr
#
# EXECUTION IN THE CURRENT DIRECTORY
#SBATCH --workdir=.
#
# PARTITION
# The following command provides the state of available partitions:
#   sinfo
#SBATCH --partition=standard
#
# ACCOUNT
# The following command provides the list of accounts you can use:
#   sacctmgr list user withassoc name=your_username format=user,account,defaultaccount
#SBATCH --account=uppa
#
# JOB MAXIMAL WALLTIME. Format: D-H:M:S
#SBATCH --time=5-0:00:00
#
# CORES NUMBER 
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#
# MEMORY PER CORE (MB)
#SBATCH --mem-per-cpu=2000
#
#
########################################### END SLURM OPTIONS ##################################
#
#
#
######################## BEGIN UNIX COMMANDS  #########################

# MODULES
module purge
module load R/4.1.1

# COMMANDS
export TMPDIR=$SCRATCHDIR
time R CMD BATCH $1

######################### END UNIX COMMANDS #########################

