#!/bin/bash
#SBATCH --job-name="hw5"
#SBATCH --partition="windfall"
#SBATCH --output="hw51.o%j"
#SBATCH --error="hw51.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -t 00:30:00

export OMP_NUM_THREADS=8
./hw5



