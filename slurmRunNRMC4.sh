#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH --account=zi             #account name
#SBATCH -J AxonMotionNode           #job name
#SBATCH --time=17:00:00                  #time the job will take to run
#SBATCH -N 1                 #if use one node, do not specify memory requirement
#SBATCH --exclusive

module load matlab
echo "launching a matlab run"
date

#Command to execute Matlab code
matlab -nosplash -nodisplay -nodesktop -r "runNRMC4, quit" 

# End of script