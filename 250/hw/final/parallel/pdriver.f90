!File: pdriver.f90
!Author: Dante Buhl
!Dependencies: gameoflife
!Description: A driver routine to run the game of life

! NOTE: Needs an input file giving some initial condition

program pdriver
  ! Parallelized Version
  use pgameoflife
  !use MPI
  
  implicit none
  

  ! General Util Vars
  integer, parameter :: ki=selected_int_kind(15)
  integer :: i, j, k,fn=10
  
  ! Game of Life Vars
  integer, allocatable :: A(:,:), G(:,:)
  integer :: m, n, nt
  
  ! MPI Variable
  integer :: ie, id, np, info
  integer(kind=MPI_OFFSET_KIND) :: offset
  !type(MPI_INFO) :: info
  !type(MPI_FILE) :: fn
  integer, allocatable :: counts(:), stat(:), disp_vec(:)
 
  ! MPI Initialization
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)  

  allocate(counts(np),disp_vec(np),stat(np))

  ! simulate the game of life
  if (id.eq.0) then
    call readmat(A, m, n) 
  end if
  
  call MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ie)
  call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ie)

  ! Decompose the tasks and distribute to available processors 
  call mpi_decomp_1d(id, np, n, counts)
  counts = counts*m

  disp_vec(1) = 0
  do i = 1, np-1
    disp_vec(i+1) = disp_vec(i)+counts(i)
  end do 
  allocate(G(m, num_tasks+2))
  G = 0
  !call preadmat(G(:,2:num_tasks+1), m, id)

  ! Scatter Columns to other processors
  call MPI_SCATTERV(A, counts, disp_vec, MPI_INTEGER, G(:,2:(num_tasks+1)),&
    num_tasks*m, MPI_INTEGER, 0, MPI_COMM_WORLD, ie)

  ! Attempting MPI Parallel IO
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'help.dat',MPI_MODE_CREATE+MPI_MODE_WRONLY,&
    MPI_INFO_NULL, fn, ie)
  offset = 0
  call MPI_FILE_SET_VIEW(fn,offset,MPI_INT,MPI_INT,"native",MPI_INFO_NULL,ie)
  if (id.eq.0) then
    call MPI_FILE_WRITE(fn, id, 1, MPI_INTEGER,MPI_STATUS_IGNORE,ie)
  end if
  !call MPI_FILE_WRITE_ALL(fn,id,1,MPI_BYTE,MPI_STATUS_IGNORE,ie)
  ! END MPI Paralle IO

  !call MPI_BARRIER(MPI_COMM_WORLD, ie)
  ! ------------------------- Sequential Component
  !call pprintmat(G(:,2:num_tasks+1), m, id)
  !nt = 200
  !do i = 1, nt
    !call pupdate_bound_1d(G, n, id, np, i*100)
    !call update(G, m, num_tasks+2)
    !call pprintmat(G(:,2:num_tasks+1), m, id)
  !end do 
  ! ------------------------- END Sequential Component
  call MPI_FILE_CLOSE(fn, ie)
  ! Stop MPI
  call MPI_FINALIZE(ie)
  deallocate(counts, disp_vec, G)

  ! Branched Logic so this works while still debugging 
  !if (id.eq.0) then 
    !deallocate(A)
    !open(fn+1, file='params.dat')
      !write(fn+1, "(3I6)") m, n,nt+1
    !close(fn+1)
  !end if
end program pdriver


