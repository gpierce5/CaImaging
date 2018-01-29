#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH --account=zi             #account name
#SBATCH -J AxonMotion           #job name
#SBATCH --time=17:00:00                  #time the job will take to run
#SBATCH --exclusive              #something you need for using a node
#SBATCH --nodes=1                   #use one node, do not specify memory requirement

module load matlab
echo "launching a matlab run"
date

#Command to execute Matlab code
matlab -nosplash -nodisplay -nodesktop -r "runNRMC, quit" 

# End of script