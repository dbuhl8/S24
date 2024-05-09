#!/bin/bash
#SBATCH --job-name="SendRec"
#SBATCH --partition="Instruction"
#SBATCH --output="sendrec1.o%j"
#SBATCH --error="sendrec1.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH -t 00:05:00

module load ohpc hb hp-gnu
module load gnu/5.4.0
module load openmpi

mpirun -np 2 sendrec_mpi
