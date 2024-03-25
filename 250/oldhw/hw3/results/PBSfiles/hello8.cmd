

#PBS -q newest
#PBS -N hello8
#PBS -l nodes=2:ppn=4
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=8
./hello_omp_grape

