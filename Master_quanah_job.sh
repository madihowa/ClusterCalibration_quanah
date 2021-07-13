#!/bin/bash
#SBATCH --job-name=felix_run_my_code_Master_only
#SBATCH --output=quanah_log/%x.o%j
#SBATCH --error=quanah_log/%x.e%j
#SBATCH --partition quanah
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=5
#SBATCH --mail-user=<madison.howard@ttu.edu>
#SBATCH --mail-type=ALL

# actual code to execute
python Master.py  "felix_run_my_code_Master_only"

