#!/bin/bash
#SBATCH --job-name="Latency"
#SBATCH --partition="windfall"
#SBATCH --output="ping1.o%j"
#SBATCH --error="ping1.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH -t 00:05:00

module load ohpc hb hp-gnu
module load gnu/5.4.0
module load openmpi

mpirun -np 2 latency.exe
