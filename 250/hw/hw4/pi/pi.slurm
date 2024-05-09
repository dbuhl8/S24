#!/bin/bash
#SBATCH --job-name="Pi"
#SBATCH --partition="Instruction"
#SBATCH --output="pi1.o%j"
#SBATCH --error="pi1.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -t 00:05:00

module load ohpc hb hp-gnu
module load gnu/5.4.0
module load openmpi

mpirun -np 24 pi_mpi

