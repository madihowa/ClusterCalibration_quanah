#!/bin/bash
#SBATCH --job-name=fixed_predict_NN
#SBATCH --output=quanah_log/%x.o%j
#SBATCH --error=quanah_log/%x.e%j
#SBATCH --partition quanah
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=5
#SBATCH --mail-user=<madison.howard@ttu.edu>
#SBATCH --mail-type=ALL

# actual code to execute
./setUpQuanah.sh
python predict.py Results_2021-07-05_run_7/ & 
python predict.py Results_2021-07-05_run_8/ &


