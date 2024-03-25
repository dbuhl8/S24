


#PBS -q newest 
#PBS -N hello
#PBS -j oe
#PBS -l nodes=2:ppn=8
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
mpirun -np 16 hello_mpi
