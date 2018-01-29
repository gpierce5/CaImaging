#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH --account=zi             #account name
#SBATCH -J AxonMotion1           #job name
#SBATCH --time=17:00:00                  #time the job will take to run
#SBATCH -c 16                  #if use one node, do not specify memory requirement
#SBATCH --mem-per-cpu=8gb        # The memory the job will use per cpu core.

module load matlab
echo "launching a matlab run"
date

#Command to execute Matlab code
matlab -nosplash -nodisplay -nodesktop -r "runNRMC3, quit" 

# End of script