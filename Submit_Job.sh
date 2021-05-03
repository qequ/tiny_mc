#!/bin/bash

#SBATCH --job-name=TinyMCAlvaroLautaro
#SBATCH --nodes=1
#SBATCH --gres=gpu:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28

### Environment setup
. /etc/profile

export MKL_NUM_THREADS=28

srun ./automate_compilation.sh
