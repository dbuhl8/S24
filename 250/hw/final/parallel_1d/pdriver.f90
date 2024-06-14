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
 
  ! Row Decomp (if yes, set this to true)
  logical :: rowdecomp = .true.
  ! General Util Vars
  integer, parameter :: ki=selected_int_kind(15)
  integer :: i, j, k,fn=10
  
  ! Game of Life Vars
  integer, allocatable :: A(:,:), G(:,:)
  integer :: m, n, nt
  
  ! MPI Variable
  integer :: ie, id, np, info
  integer(kind=MPI_OFFSET_KIND) :: offset
  integer :: req
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
  if (rowdecomp) then
    call mpi_decomp_1d('R',id, np, m,n, counts)
    allocate(G(num_tasks+2, n))
    G = 0
    disp_vec(1) = 1
    do i = 1, np-1
      disp_vec(i+1) = disp_vec(i)+counts(i)
    end do 
    ! send into to other processors
    if (id.eq.0) then
      do i = 1, np-1
        call MPI_SEND(A(disp_vec(i+1):disp_vec(i+1)+counts(i+1)-1,:),&
        counts(i+1)*n,MPI_INTEGER,i,rtag(i,0),MPI_COMM_WORLD,ie)
      end do 
      G(2:num_tasks+1,:) = A(1:counts(1),:)
    else
      call MPI_RECV(G(2:num_tasks+1,:),counts(id+1)*n,MPI_INTEGER,&
        0,rtag(id,0),MPI_COMM_WORLD,stat,ie)
    end if
  else 
    call mpi_decomp_1d('C',id, np, m, n, counts)
    counts = counts*m
    disp_vec(1) = 0
    do i = 1, np-1
      disp_vec(i+1) = disp_vec(i)+counts(i)
    end do 
    allocate(G(m, num_tasks+2))
    G = 0

    ! Scatter Columns to other processors
    call MPI_SCATTERV(A, counts, disp_vec, MPI_INTEGER, G(:,2:num_tasks+1),&
      num_tasks*m, MPI_INTEGER, 0, MPI_COMM_WORLD, ie)
  end if

  ! ------------------------- Sequential Component
  
  if(id.eq.0) then
    open(fn, file='gol.dat')
    call writemat(A, m, n, fn)
  end if
  
  nt = 1000
  do i = 1, nt
    if (rowdecomp) then
      call pupdate_bound_1d(G, m, id, np, i)
      call MPI_BARRIER(MPI_COMM_WORLD, ie)
      call update(G, num_tasks+2, n)
      ! IO
      if (id.eq.0) then
        do j = 1, np-1
          call MPI_RECV(A(disp_vec(j+1):disp_vec(j+1)+counts(j+1)-1,:),&
          counts(j)*n,MPI_INTEGER,j,rtag(j,i),MPI_COMM_WORLD,stat,ie)
        end do 
        A(1:counts(1),:) = G(2:num_tasks+1,:)
      else
        call MPI_SEND(G(2:num_tasks+1,:),num_tasks*n,MPI_INTEGER,&
          0,rtag(id,i),MPI_COMM_WORLD,ie)
      end if
      if (id.eq.0) then
        call writemat(A, m, n, fn)
      end if
    else
      ! column decomp
      call pupdate_bound_1d(G, m, id, np, i*100)
      call MPI_BARRIER(MPI_COMM_WORLD, ie)
      call update(G, m, num_tasks+2)
      call MPI_GATHERV(G(:,2:num_tasks+1), num_tasks*m, MPI_INTEGER, A, counts,&
        disp_vec, MPI_INTEGER, 0, MPI_COMM_WORLD, ie)
      if (id.eq.0) then
        call writemat(A, m, n, fn)
      end if
    end if
    ! End IO
  end do 
  ! ------------------------- END Sequential Component
  ! Stop MPI

  ! Branched Logic so this works while still debugging 
  if (id.eq.0) then 
    close(fn)
    deallocate(A)
    open(fn+1, file='params.dat')
      write(fn+1, "(3I6)") m, n,nt+1
    close(fn+1)
  end if

  deallocate(counts, disp_vec, G)
  call MPI_FINALIZE(ie)

end program pdriver


