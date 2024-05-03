#!/bin/bash
#SBATCH --job-name="AM250Hw3"
#SBATCH --partition="Instruction"
#SBATCH --output="HW3mpirun1.o%j"
#SBATCH --error="HW3mpirun1.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:05:00

mpirun -np 4 hello_mpi_hb
