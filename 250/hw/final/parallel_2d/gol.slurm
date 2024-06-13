#!/bin/bash
#SBATCH --job-name="GameOfLife"
#SBATCH --partition="Instruction"
#SBATCH --output="gol1.out"
#SBATCH --error="gol1.err"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:05:00

mpirun -np 8 driver.exe <IC
