#!/bin/bash

#SBATCH -p Instruction
#SBATCH -J Hello6
#SBATCH -e job6.err
#SBATCH -o job6.out
#SBATCH -N 2
#SBATCH --ntasks-per-node 4
#SBATCH -t 00:05:00

export OMP_NUM_THREADS=8
./hello_omp_hb

