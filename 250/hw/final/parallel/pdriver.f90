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
  integer :: i, j, k,fn=10
  
  ! Game of Life Vars
  integer, allocatable :: A(:,:), G(:,:)
  integer :: m, n, nt
  
  ! MPI Variable
  integer :: ie, id, np
  integer, allocatable :: counts(:), disp_vec(:)
 
  ! MPI Initialization
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)  

  allocate(counts(np),disp_vec(np))

  
  ! Initial Condition 
    !m = 10
    !n = 10 
    !allocate(A(m, n))


    !open(fn, file='IC.dat')

    !do j = 1, n
      !do i = 1, m
        !A(i, j) = anint(rand(0))
      !end do
    !end do

    !write(fn, "('# ', I6, I6)") m, n
      
  ! simulate the game of life
  if (id.eq.0) then
    call readmat(A)
  end if
  call MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ie)
  call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ie)
 
  ! Decompose the tasks and distribute to available processors 
  call mpi_decomp_1d(id, np, n, counts)
  if (id.eq.0) then
    disp_vec(1) = 1
    do i = 1, np-1
      disp_vec(i+1) = disp_vec(i)+counts(i)-1
    end do 
  end if
  allocate(G(m, num_tasks+2))

  ! Scatter Columns to other processors
  call MPI_SCATTERV(A, counts, disp_vec, MPI_INTEGER, G(:,2:(num_tasks+1)),&
    num_tasks, MPI_INTEGER, 0, MPI_COMM_WORLD, ie)

  ! Branched Logic so this works while still debugging 
  open(fn, file='gof.dat')

  ! ------------------------- Sequential Component
  call pwritemat(G(:,2:num_tasks+1), num_tasks, id, fn)
  nt = 1000
  do i = 1, nt
    call pupdate_bound_1d(G, n, id, np, i)
    call update(G, m, num_tasks+2)
    call pwritemat(G(:,2:num_tasks+1), num_tasks, id, fn)
  end do 
  ! ------------------------- END Sequential Component

  close(fn)

  ! Stop MPI
  call MPI_FINALIZE(ie)
  deallocate(counts, disp_vec, G)

  ! Branched Logic so this works while still debugging 
  if (id.eq.0) then 
    deallocate(A)
    open(fn+1, file='params.dat')
      write(fn+1, "(3I6)") n, m, nt+1
    close(fn+1)
  end if
end program pdriver


