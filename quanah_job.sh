#!/bin/bash
#SBATCH --job-name=testRun
#SBATCH --output=quanah_log/%x.o%j
#SBATCH --error=quanah_log/%x.e%j
#SBATCH --partition quanah
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --mem-per-cpu=3994MB
#SBATCH --mail-user=<madison.howard@ttu.edu>
#SBATCH --mail-type=ALL

# actual code to execute
./run.sh "testRun"

