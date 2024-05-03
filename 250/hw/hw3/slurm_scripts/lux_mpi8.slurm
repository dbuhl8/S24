#!/bin/bash
#SBATCH  --job-name="AM250Hw3"
#SBATCH  --partition="windfall"
#SBATCH  --output="lux_mpi_2node_8cpus.out"
#SBATCH  --error="lux_mpi_2node_8cpus.err"
#SBATCH  --nodes=2
#SBATCH  --ntasks-per-node=4
#SBATCH  -t 00:05:00

module load gcc
module load devtoolset-9
module load openmpi/4.1.5

mpirun -np 8 hello_mpi_lux

