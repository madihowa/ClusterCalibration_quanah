#!/bin/bash
#SBATCH --job-name=step3_trial1
#SBATCH --output=quanah_log/%x.o%j
#SBATCH --error=quanah_log/%x.e%j
#SBATCH --partition quanah
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=5

# actual code to execute
./run.sh "step3_trial1"

