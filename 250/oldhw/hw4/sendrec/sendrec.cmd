


#PBS -q newest 
#PBS -N sendrec
#PBS -o sendrecO
#PBS -j oe
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
mpif90 -o sendrec_mpi sendrec.f90
mpirun -np 2 sendrec_mpi
