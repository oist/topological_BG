#!/bin/bash

#SBATCH --time=2-00:00:00 
#SBATCH --partition=compute 
#SBATCH --mem-per-cpu=4G 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10 
#SBATCH --job-name=sBCBG_2021_05_03_16_58_22_xp000000
#SBATCH --input=none
#SBATCH -e bg_job_e.out 
#SBATCH --output=bg_job_o.out 
module load openmpi.gcc/4.0.3 
module load python/3.7.3 
source /apps/unit/DoyaU/nest-2.20-install/bin/nest_vars.sh 
echo "path: " $PATH 
echo "pythonpath: " $PYTHONPATH 
echo "cores: " $OMP_NUM_THREADS 
srun --mpi=pmi2 python3 stim_all_model.py 
