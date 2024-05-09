#!/bin/bash
#SBATCH --job-name="ring"
#SBATCH --partition="Instruction"
#SBATCH --output="ring1.o%j"
#SBATCH --error="ring1.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -t 00:05:00

module load ohpc hb hp-gnu
module load gnu/5.4.0
module load openmpi

mpirun -np 12 ring_mpi
