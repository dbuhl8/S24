
#PBS -q newest
#PBS -N hello7
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=4
./hello_omp_grape
