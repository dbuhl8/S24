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
  integer :: ie, id, np, info, body_type, core_type
  integer(kind=MPI_OFFSET_KIND) :: offset
  integer, allocatable :: counts(:,:), stat(:), disp_vec(:,:)
 
  ! MPI Initialization
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)  

  allocate(counts(np,2),disp_vec(np,2),stat(np))

  ! simulate the game of life
  if (id.eq.0) then
    call readmat(A, m, n) 
    !call getprobsize(m, n)
  end if
  
  call MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ie)
  call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ie)

  ! Decompose the tasks and distribute to available processors 
  call mpi_decomp_2d(id, np, m, n, counts)
  if (id.eq.0) then
    do i = 1, np
      print *, i, counts(i,1), counts(i,2)
    end do 
  end if
  allocate(G(counts(id+1,1)+2,counts(id+1,2)+2))
  G = 0.

  do i = 1, nm
    do j = 1, nn
      k = nm*(j-1)+i
      if(i .eq. 1) then
        disp_vec(k,1) = 0
      else
        disp_vec(k,1) = disp_vec(k-1,1)+counts(k-1,1)
      end if
      if(j .eq. 1) then
        disp_vec(k,2) = 0
      else 
        disp_vec(k,2) = disp_vec(k-nm,2)+counts(k-nm,2)
      end if
    end do 
  end do 

  if (id.eq.0) then
    do i = 1, np
      print *, i, disp_vec(i,1), disp_vec(i,2)
    end do 
  end if

  ! Create subarray MPI types (not sure if edge and corner type are needed yet)
  ! Edge Type
  ! Corner Type

  ! Body Type
  !                           ndims glbsz   subsz       starts      
  call MPI_TYPE_CREATE_SUBARRAY(2,(/m,n/),counts(id+1,:),disp_vec(id+1,:),&
    MPI_ORDER_FORTRAN, MPI_INTEGER, body_type, ie)
  !     order            datatype   new type   ie 
  call MPI_TYPE_COMMIT(body_type,ie)
  ! Global Data
  call MPI_TYPE_CREATE_SUBARRAY(2,counts(id+1,:)+2,counts(id+1,:),(/2,2/),&
    MPI_ORDER_FORTRAN, MPI_INTEGER, core_type,ie)
  call MPI_TYPE_COMMIT(core_type,ie)
  print *, "got past type declaration ", id

  !                     comm       filename     mode         info      
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'ICp',MPI_MO   DE_RDONLY,MPI_INFO_NULL,&
    fn,ie)
  print *, "got past open", id
  offset = 0
  call MPI_FILE_SET_VIEW(fn,offset,MPI_INTEGER,body_type,"native",MPI_INFO_NULL,ie)
  print *, "got set view ", id
  !                      fh buf counts, type, data, ie
  call MPI_FILE_READ_ALL(fn,G,nm*nn,core_type,stat,ie)
  !call MPI_FILE_READ_ALL(fn,int(disp_vec(id+1,:),ki), G,nm*nn,body_type,stat,ie)
  !
  print *, "got past read ", id

  if (id .eq. 0) then
    call printmat(G, counts(id+1,1)+2, counts(id+1,2)+2)
  end if
  print *, "got past print", id

  call MPI_FILE_CLOSE(fn, ie)
  print *, "got past close", id
  fn = fn + 1


  !call MPI_FILE_OPEN(MPI_COMM_WORLD,'gol.dat',MPI_CREATE+MPI_WRITE,&
    !MPI_INFO_NULL, fn)


  ! Scatter Columns to other processors
  !call MPI_SCATTERV(A, counts(:,1)*counts(:,2), disp_vec, MPI_INTEGER, G,&
    !nm*nn, core_type, 0, MPI_COMM_WORLD, ie)

  ! ------------------------- Sequential Component
  !if(id.eq.0) then
    !open(fn, file='gol.dat')
    !call writemat(A, m, n, fn)
  !end if
  
  !nt = 1000
  !do i = 1, nt
    !call pupdate_bound_1d(G, m, id, np, i*100)
    !call MPI_BARRIER(MPI_COMM_WORLD, ie)
    !write (fn+id+2, *) " "
    !call update(G, m, num_tasks+2)
    ! IO Step
    !call MPI_GATHERV(G(:,2:num_tasks+1), num_tasks*m, MPI_INTEGER, A, counts,&
      !disp_vec, MPI_INTEGER, 0, MPI_COMM_WORLD, ie)
    !if (id.eq.0) then
      !call writemat(A, m, n, fn)
    !end if
    ! End IO
  !end do 
  ! ------------------------- END Sequential Component

  !if (id.eq.0) then 
    !close(fn)
    !deallocate(A)
    !open(fn+1, file='params.dat')
      !write(fn+1, "(3I6)") m, n,nt+1
    !close(fn+1)
  !end if

  call MPI_FINALIZE(ie)
  deallocate(counts, disp_vec, stat, G)
  !deallocate(counts, disp_vec, stat)

end program pdriver


