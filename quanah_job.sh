#!/bin/bash
#SBATCH --job-name=new_data_cluster_SECOND_TIME_PTD
#SBATCH --output=quanah_log/%x.o%j
#SBATCH --error=quanah_log/%x.e%j
#SBATCH --partition quanah
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=5
#SBATCH --mail-user=<madison.howard@ttu.edu>
#SBATCH --mail-type=ALL

# actual code to execute
./run.sh "new_data_cluster_SECOND_TIME_PTD"

