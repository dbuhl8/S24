#!/bin/bash
#SBATCH  --job-name="AM250Hw3"
#SBATCH  --partition="windfall"
#SBATCH  --output="lux_omp_1node_4cpus.out"
#SBATCH  --error="lux_omp_1node_4cpus.err"
#SBATCH  --nodes=1
#SBATCH  --ntasks-per-node=4
#SBATCH  -t 00:05:00

module load gcc
module load devtoolset-9
module load openmpi/4.1.5

export OMP_NUM_THREADS=4
./hello_omp_lux

