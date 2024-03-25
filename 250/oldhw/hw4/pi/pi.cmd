


#PBS -q newest
#PBS -N pi
#PBS -o piO
#PBS -j oe
#PBS -l nodes=5:ppn=10
#PBS -l walltime=00:05:00


cd $PBS_O_WORKDIR
mpif90 -o pi_mpi pi.f90
mpirun -np 50 pi_mpi

