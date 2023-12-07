#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --partition=gpu
#SBATCH --nodelist=d14-17
#SBATCH --output=project_1.out

./async_read_pthreads 1
