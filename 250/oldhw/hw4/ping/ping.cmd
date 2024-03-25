


#PBS -q newest
#PBS -N ping
#PBS -o pingO
#PBS -j oe
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:01:00

cd $PBS_O_WORKDIR
mpif90 -o ping_mpi ping.f90
mpirun -np 2 ping_mpi
