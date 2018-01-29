#!/bin/sh
#
# Simple Matlab submit script for Slurm.
#
#
#SBATCH --account=zi             #account name
#SBATCH -J AxonMotion           #job name
#SBATCH --time=17:00:00                  #time the job will take to run
#SBATCH --mem=128gb        #memory per node
#SBATCH -N 1                   #use one node


echo "Launching a Matlab run"
date

#Command to execute Matlab code
matlab -nosplash -nodisplay -nodesktop -r "runNRMC(24)" # > matoutfile

# End of script