#!/bin/bash
#SBATCH --job-name="AM250Hw3"
#SBATCH --partition="128x24"
#SBATCH --output="HW3omprun1.o%j"
#SBATCH --error="HW3omprun1.e%j"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:05:00

export OMP_NUM_THREADS=8
./hello_omp_hb
