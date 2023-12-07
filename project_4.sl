#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --partition=gpu
#SBATCH --nodelist=d14-15
#SBATCH --output=project_4.out

./async_read_pthreads 4
