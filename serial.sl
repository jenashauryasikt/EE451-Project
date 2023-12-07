#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --nodelist=d14-18 
#SBATCH --output=project_serial.out
#SBATCH --error=project_serial.err

./serial
