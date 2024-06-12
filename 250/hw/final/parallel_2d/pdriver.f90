!File: pdriver.f90
!Author: Dante Buhl
!Dependencies: gameoflife
!Description: A driver routine to run the game of life

! NOTE: Needs an input file giving some initial condition
! i.e. mpirun -np 4 driver.exe < IC

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
  integer, allocatable :: cnt(:,:), stat(:), dv(:,:)
 
  ! MPI Initialization
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)  

  ! Allocating Memory
  allocate(cnt(np,2),dv(np,2),stat(np))

  ! Reading in Initial Condition
    if (id.eq.0) then
      call readmat(A, m, n) 
      !call getprobsize(m, n)
    end if
  !   ---------------------------------------
 
  ! Broadcasting global array size accross processors 
  call MPI_BCAST(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ie)
  call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ie)

  ! Decompose the tasks and distribute to available processors 
  call mpi_decomp_2d(id, np, m, n, cnt)

  ! Allocating local memory
  allocate(G(cnt(id+1,1)+2,cnt(id+1,2)+2))
  G = 0

  ! Compute the Disp_vec array
    do i = 1, nm
      do j = 1, nn
        k = nm*(j-1)+i
        if(i .eq. 1) then
          dv(k,1) = 0
        else
          dv(k,1) = dv(k-1,1)+cnt(k-1,1)
        end if
        if(j .eq. 1) then
          dv(k,2) = 0
        else 
          dv(k,2) = dv(k-nm,2)+cnt(k-nm,2)
        end if
      end do 
    end do 
  !   ---------------------------------------

  ! ---------Sending info to processors-----------
    dv = dv+1
    if(id.eq.0) then
      do i = 2, np
        ! Debug print statements to check indices
          !print *, "Sending array ", i-1
          !print *, "Sending Indices ", dv(i,1), dv(i,1)+cnt(i,1)-1,&
            !dv(i,2), dv(i,2)+cnt(i,2)-1
        ! Debug ----------------------------------------

        call MPI_SEND(A(dv(i,1):dv(i,1)+cnt(i,1)-1,&
          dv(i,2):dv(i,2)+cnt(i,2)-1),&
          product(cnt(i,:)), MPI_INTEGER, i-1,i,MPI_COMM_WORLD,ie)
      end do 
      ! Giving itself its own portion
      G(2:cnt(1,1)+1,2:cnt(1,2)+1) = A(1:cnt(1,1),1:cnt(1,2))
    else 
      ! Receive matrix on other processors
      call MPI_RECV(G(2:cnt(id+1,1)+1,2:cnt(id+1,2)+1),product(cnt(id+1,:)),&
        MPI_INTEGER, 0,id+1,MPI_COMM_WORLD,stat,ie)
    end if
  !   ---------------------------------------

  ! Debug: print statements for debuggins
    !if (id.eq.3) then
      !do i = 1, np
        !print *, i, dv(i,1), dv(i,2)
        !print *, i, cnt(i,1), cnt(i,2)
      !end do 
    !end if
    !print *, id, dv(id+1,1), dv(id+1,2)
    !print *, id, cnt(id+1,1), cnt(id+1,2)

    ! MPI IO / MPI_SCATTER Not currently working  
    ! Create subarray MPI types (not sure if edge and corner type are needed yet)
    ! Body Type
    !                           ndims glbsz   subsz       starts      
    !call MPI_TYPE_CREATE_SUBARRAY(2,(/m,n/),cnt(id+1,:),dv(id+1,:),&
      !MPI_ORDER_FORTRAN, MPI_INTEGER, body_type, ie)
    !     order            datatype   new type   ie 
    !call MPI_TYPE_COMMIT(body_type,ie)
    ! Global Data
    !call MPI_TYPE_CREATE_SUBARRAY(2,cnt(id+1,:)+2,cnt(id+1,:),(/2,2/),&
      !MPI_ORDER_FORTRAN, MPI_INTEGER, core_type,ie)
    !call MPI_TYPE_COMMIT(core_type,ie)

    !                     comm       filename     mode         info      
    !call MPI_FILE_OPEN(MPI_COMM_WORLD,'ICp',MPI_MODE_RDONLY,MPI_INFO_NULL,&
      !fn,ie)
    !offset = 0
    !call MPI_FILE_SET_VIEW(fn,offset,MPI_INTEGER,body_type,"native",MPI_INFO_NULL,ie)
    !                      fh buf cnt, type, data, ie
    !call MPI_FILE_READ_ALL(fn,G,nm*nn,core_type,stat,ie)
    !call MPI_FILE_READ_AT_ALL(fn,int(dv(id+1,:),ki), G,nm*nn,core_type,stat,ie)
    !

    !call MPI_FILE_CLOSE(fn, ie)
    !print *, "got past close", id
    !fn = fn + 1


    !call MPI_FILE_OPEN(MPI_COMM_WORLD,'gol.dat',MPI_CREATE+MPI_WRITE,&
      !MPI_INFO_NULL, fn)


    ! Scatter Columns to other processors
    !call MPI_SCATTERV(A, cnt(:,1)*cnt(:,2), dv, MPI_INTEGER, G,&
      !nm*nn, core_type, 0, MPI_COMM_WORLD, ie)
    !call MPI_SCATTERV(A, cnt(:,1)*cnt(:,2), dv, MPI_INTEGER, G,&
      !nm*nn, core_type, 0, MPI_COMM_WORLD, ie)
    ! Making sure arrays are dsitributed correctly.
    !call MPI_BARRIER(MPI_COMM_WORLD,ie)
    !if (id .eq. 0) then
      !print *, " ", id
      !call printmat(G, cnt(id+1,1)+2, cnt(id+1,2)+2)
      !open(fn+id+1, file=''//trim(str(id))//'.dat')
      !call writemat(G, ntm+2, ntn+2,fn+id+1)
      !close(fn+id+1)
      !print *, " "
    !end if
    !call MPI_BARRIER(MPI_COMM_WORLD,ie)
    !if (id .eq. 1) then
      !print *, " ", id
      !call printmat(G, cnt(id+1,1)+2, cnt(id+1,2)+2)
      !open(fn+id+1, file=''//trim(str(id))//'.dat')
      !call writemat(G, ntm+2, ntn+2,fn+id+1)
      !close(fn+id+1)
      !print *, " "
    !end if
    !call MPI_BARRIER(MPI_COMM_WORLD,ie)
    !if (id .eq. 2) then
      !print *, " ", id
      !call printmat(G, cnt(id+1,1)+2, cnt(id+1,2)+2)
      !open(fn+id+1, file=''//trim(str(id))//'.dat')
      !call writemat(G, ntm+2, ntn+2,fn+id+1)
      !close(fn+id+1)
      !print *, " "
    !end if
    !call MPI_BARRIER(MPI_COMM_WORLD,ie)
    !if (id .eq. 3) then
      !print *, " ", id
      !call printmat(G, cnt(id+1,1)+2, cnt(id+1,2)+2)
      !open(fn+id+1, file=''//trim(str(id))//'.dat')
      !call writemat(G, ntm+2, ntn+2,fn+id+1)
      !close(fn+id+1)
      !print *, " "
    !end if
    !call MPI_BARRIER(MPI_COMM_WORLD,ie)
    !print *, "got past print", id

    ! Making sure mapping array is computed correctly
    !print *, "ID: ", id, ", pmap :", pmap
    !print *, "ID: ", id, ", Num Row Tasks: ", ntm, ", Num Col Tasks:", ntn
    !open(fn+id+1, file=''//trim(str(id))//'.dat')
    !call writemat(G, ntm+2, ntn+2, fn+id+1)
    !write(fn+id+1, *) " "
  ! Debug ---------------------------------------
 
  ! ------------------------- Output and Update (timestepping)
    if(id.eq.0) then
      open(fn, file='gol.dat')
      call writemat(A, m, n, fn)
    end if
   
    nt = 1000
    do i = 1, nt
      ! there is a memory bug in pupdate_bound_2d
      call pupdate_bound_2d(G, id, i)
      !print *, "id:",id,",got here"
      call MPI_BARRIER(MPI_COMM_WORLD, ie)
      !print *, "id:",id,",got here"
      call update_2d(G, ntm+2, ntn+2)
      !print *, "id:",id,",got here"
      call MPI_BARRIER(MPI_COMM_WORLD, ie)
      !print *, "id:",id,",got here"
      ! IO Step
      ! This will probably be replaced with send/recv instead of gathers for now
      !call MPI_GATHERV(G(:,2:num_tasks+1), num_tasks*m, MPI_INTEGER, A, cnt,&
        !dv, MPI_INTEGER, 0, MPI_COMM_WORLD, ie)
      if (id.eq.0) then
        do j = 1, np-1
          call MPI_RECV(A(dv(j+1,1):dv(j+1,1)+cnt(j+1,1)-1,dv(j+1,2):dv(j+1,2)&
            +cnt(j+1,2)-1),product(cnt(j+1,:)),MPI_INTEGER,j,rtag(j,i),&
            MPI_COMM_WORLD,stat,ie)
          print *, "Recevied from processor :", j, " on step :", i
        end do 
        A(1:cnt(1,1),1:cnt(1,2)) = G(2:cnt(1,1)+1,2:cnt(1,2)+1)
      else 
        call MPI_SEND(G(2:ntm+1,2:ntn+1),ntm*ntn,MPI_INTEGER,0,rtag(id,i),&
          MPI_COMM_WORLD,ie)
      end if
      if (id.eq.0) then
        call writemat(A, m, n, fn)
      end if
      ! Debug ---------------------------------------
        !call writemat(G, ntm+2, ntn+2, fn+id+1)
        !write(fn+id+1, *) "  "
        !write(fn+id+1, *) " i :", i
        !write(fn+id+1, *) " "
      ! Debug ---------------------------------------
      ! End IO
    end do 
  ! ------------------------- END TIMESTEPPING --------------

  ! Finising data output for matlab
  if (id.eq.0) then 
    close(fn)
    ! just in case I forget to uptick fn
    fn = fn + np + 1
    open(fn, file='params.dat')
      write(fn, "(3I6)") m, n,nt+1
    close(fn)
  end if

  ! Deallocating and closing
    deallocate(cnt, dv, stat, G)
    if(id.eq.0) then
      deallocate(A)
    end if
  ! -------------------------
  ! Ending MPI
  call MPI_FINALIZE(ie)
end program pdriver

