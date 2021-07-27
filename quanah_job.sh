#!/bin/bash
#SBATCH --job-name=initial_inputs2372021
#SBATCH --output=quanah_log/%x.o%j
#SBATCH --error=quanah_log/%x.e%j
#SBATCH --partition nocona
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --mem-per-cpu=3994MB
#SBATCH --mail-user=<madison.howard@ttu.edu>
#SBATCH --mail-type=ALL
# actual code to execute
./run.sh "initial_inputs"
