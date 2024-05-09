!File: helloworld/f90
!Author: Dante Buhl
!Dependencies: MPI

program helloworld

  include 'mpif.h'

  integer ierr, myid, numprocs

  ! Start MPI 
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

  ! Print Hello from each processor
  print "(A, I3, A, I3)","Hello from processor #",myid," out of ",numprocs

  ! Stop MPI
  call MPI_FINALIZE(ierr)

end program helloworld
