


#PBS -q newest
#PBS -N ring
#PBS -o ringO
#PBS -j oe
#PBS -l nodes=4:ppn=3
#PBS -l walltime=00:02:00

cd $PBS_O_WORKDIR
mpif90 -o ring_mpi ring.f90
mpirun -np 12 ring_mpi
