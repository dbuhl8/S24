#!/bin/bash
#SBATCH --job-name="hw5"
#SBATCH --partition="Instruction"
#SBATCH --output="hw51.o%j"
#SBATCH --error="hw51.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 00:30:00

module load ohpc hb hp-gnu
module load gnu/5.4.0
module load openmpi

export OMP_NUM_THREADS=16
./hw5



