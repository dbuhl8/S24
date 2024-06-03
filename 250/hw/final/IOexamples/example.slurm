#!/bin/bash
#SBATCH --job-name="example"
#SBATCH --partition="Instruction"
#SBATCH --output="example1.o%j"
#SBATCH --error="example1.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH -t 00:05:00

module load ohpc hb hp-gnu
module load gnu/5.4.0
module load openmpi

mpif90 example.f90 -o exam
mpirun -np 6 exam
