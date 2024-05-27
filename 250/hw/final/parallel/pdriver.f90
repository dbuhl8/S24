!File: pdriver.f90
!Author: Dante Buhl
!Dependencies: gameoflife
!Description: A driver routine to run the game of life

! NOTE: Needs an input file giving some initial condition

program pdriver
  ! Parallelized Version
  use pgameoflife
  use MPI
  
  implicit none

  ! General Util Vars
  integer, parameter :: ki=selected_int_kind(15)
  integer :: i, j, k,fn=10
  
  ! Game of Life Vars
  integer, allocatable :: A(:,:), G(:,:)
  integer :: m, n, nt
  
  ! MPI Variable
  integer :: ie, id, np
  integer, allocatable :: counts(:), stat(:), disp_vec(:)
 
  ! MPI Initialization
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)  

  allocate(counts(np),disp_vec(np),stat(np))
  
  ! simulate the game of life
  if (id.eq.0) then
    call readmat(A, m, n) 
    !call getprobsize(m, n)
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
    !do i = 0, np-1
      !if(id.eq.i) then
        !call printmat(G(:,2:num_tasks+1), m, num_tasks)
        !print *, " "
      !end if
      !call MPI_BARRIER(MPI_COMM_WORLD, ie)
    !end do 

    ! Branched Logic so this works while still debugging 
    !open(fn, file='gof.dat') 
    !checkid = 2
    !call MPI_FILE_OPEN(MPI_COMM_WORLD,"gof.dat",MPI_MODE_CREATE+MPI_MODE_WRONLY,& 
      !MPI_INFO_NULL,fn,ie)
    !if(id.eq.checkid) then
      !print *, trim(str(id))//": Got here "
    !end if
    !call MPI_TYPE_CREATE_SUBARRAY(2, (/m, n/), (/m, num_tasks/),&
      !(/0,disp_vec(id+1)/m/), MPI_ORDER_FORTRAN, MPI_INTEGER, viewtype, ie)
    !if(id.eq.checkid) then
      !print *, trim(str(id))//": Got here "
    !end if
    !call MPI_TYPE_COMMIT(viewtype, ie)
    !if(id.eq.checkid) then
      !print *, trim(str(id))//": Got here "
    !end if
    !call MPI_TYPE_CREATE_SUBARRAY(2,(/m,num_tasks+2/),(/m,num_tasks/),(/0,1/),&
      !MPI_ORDER_FORTRAN,MPI_INTEGER,writetype, ie)
    !if(id.eq.checkid) then
      !print *, trim(str(id))//": Got here "
    !end if
    !call MPI_TYPE_COMMIT(writetype, ie)
    !if(id.eq.checkid) then
      !print *, trim(str(id))//": Got here "
    !end if
    !call MPI_BARRIER(MPI_COMM_WORLD, ie)
    !if(id.eq.checkid) then
      !print *, trim(str(id))//": Got here "
    !end if
    !call MPI_FILE_SET_VIEW(fn,int(0,ki),MPI_INTEGER,viewtype,"native",&
      !MPI_INFO_NULL,ie)
    !if(id.eq.checkid) then
      !print *, trim(str(id))//": Got here "
    !end if
    !call MPI_FILE_WRITE_ALL(fn,G,num_tasks*m,writetype,stat,ie)
  ! END MPI Paralle IO

  call MPI_BARRIER(MPI_COMM_WORLD, ie)
  ! ------------------------- Sequential Component
  call pprintmat(G(:,2:num_tasks+1), m, id)
  nt = 200
  do i = 1, nt
    call pupdate_bound_1d(G, n, id, np, i*100)
    call update(G, m, num_tasks+2)
    call pprintmat(G(:,2:num_tasks+1), m, id)
  end do 
  ! ------------------------- END Sequential Component

  ! Stop MPI
  call MPI_FINALIZE(ie)
  deallocate(counts, disp_vec, G)

  ! Branched Logic so this works while still debugging 
  if (id.eq.0) then 
    deallocate(A)
    open(fn+1, file='params.dat')
      write(fn+1, "(3I6)") m, n,nt+1
    close(fn+1)
  end if
end program pdriver


