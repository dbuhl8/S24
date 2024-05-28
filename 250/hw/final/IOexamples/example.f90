! This is a file that Nic sent to me so that I could get an example of working
! MPI IO. It apparently comes from Indiana University. 


program  create_file
  
  !**************************************************************************
  !  This is a Fortran 90 program to write data directly to a file by each
  !  member of an MPI group.  It is suitable for large jobs which will not
  !  fit into core memory (such as "out of core" solvers)  
  !
  !  Copyright by the Trustees of Indiana University 2005
  ***************************************************************************

  USE MPI
  
  implicit none

  integer, parameter :: kind_val = 4
  integer, parameter  :: filesize = 40
  integer :: realsize = 4
  integer ::  rank, ierr, fh, nprocs, num_reals
  integer ::  i, region
  real (kind = kind_val) :: datum
  integer, dimension (MPI_STATUS_SIZE) :: status
  integer (kind = MPI_OFFSET_KIND) :: offset, empty

!  Set filename to output datafile

  character (len = *), parameter :: filename = "/u/ac/rays/new_data.dat"
  real (kind = kind_val), dimension ( : ), allocatable  :: bucket

!  Basic MPI set-up

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!  Sanity print

  print*, "myid is ", rank

!  Carve out a piece of the output file and create a data bucket

  empty = 0
  region = filesize / (nprocs )
  offset = ( region * rank )
  allocate (bucket(region))

!  There is no guarantee that an old file will be clobbered,  so wipe out any previous output file

  if (rank .eq. 0) then
          call MPI_File_delete(filename, MPI_INFO_NULL, ierr)
  endif

!  Set the file handle to an initial value (this should not be required)

  fh = 0

!  Open the output file

  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL, fh, ierr)

!  Wait on everyone to catch up.

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!  Do some work and fill up the data bucket

  call random_seed()

  do i = 1, region

    call random_number(datum)

    bucket(i) = datum * 1000000. * (rank + 1)

    print *, " bucket  ",i ,"= ", bucket(i)
  enddo

!  Basic "belt and suspenders insurance that everyone's file pointer is at the beginning of the output file.

  call MPI_FILE_SET_VIEW(fh, empty, MPI_REAL4, MPI_REAL4, 'native', MPI_INFO_NULL, ierr)

!  Send the data bucket to the output file in the proper place

  call MPI_FILE_WRITE_AT(fh, offset, bucket, region, MPI_REAL4, status, ierr)

!  Wait on everyone to finish and close up shop

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  call MPI_FILE_CLOSE(fh, ierr)

  call MPI_FINALIZE(ierr)

end  program  create_file
