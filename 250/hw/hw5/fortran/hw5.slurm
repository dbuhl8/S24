#!/bin/bash
#SBATCH --job-name="Hello"
#SBATCH --partition="Instruction"
#SBATCH --output="hello1.o%j"
#SBATCH --error="hello1.e%j"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:05:00

module load ohpc hb hp-gnu
module load gnu/5.4.0
module load openmpi



