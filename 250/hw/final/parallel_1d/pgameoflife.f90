!File: gameoflife.f90
!Author: Dante Buhl
!Dependencies:
!Description: A module containing all of the source code to run the game of life

module pgameoflife
  ! Divided into serial and parallel versions
  use MPI
  implicit none

  ! General Utility Vars
  integer :: sm, em, sn, en, num_tasks, num_procs
  integer :: nm, nn, npm, npn, ntm, ntn
  integer, dimension(8) :: pmap
  integer :: row_type, row_type2
  integer :: g2s_type, s2s_type
  logical :: row_decomp
  
  ! Subroutines and Functions
  contains
    subroutine getprobsize(m,n)
      ! Reads m,n from stin file
      integer :: n, m
      read(*,"(2I6)") m, n
    end subroutine getprobsize

    subroutine readmat(A, m, n)
      ! Returns an allocated array A, of size m x n with data read in from the
      ! default input file.
      implicit none
      
      integer, allocatable :: A(:,:)
      integer :: m, n, i

      read(*, "(2I6)") m, n
      allocate(A(m,n))

      do i = 1, n
        read(*, "("//trim(str(m))//"I1)") A(:,i)
      end do 
    end subroutine

    subroutine printmat(A, m, n)
      implicit none

      integer :: A(:,:), m, n, i

      do i = 1, m
        print "("//trim(str(n))//"(I1, ' '))", A(i,:)
      end do 
    end subroutine printmat

    subroutine writemat(A, m, n, fn) 
      !writes to file number, fn
      integer :: A(:,:), m, n, fn, i

      do i = 1, n
        write(fn, "("//trim(str(m))//"(I1, ' '))") A(:,i)
      end do 
    end subroutine

    subroutine update(A, m, n)
      ! Note that this needs to be upgraded for non square matrices. 
      implicit none
      integer, intent(in) :: m,n
      integer :: A(:, :)
      integer, dimension(m, n) :: B
      integer, dimension(n) :: nplus, nminus
      integer, dimension(m) :: mplus, mminus
      integer :: i, j

      !print *, "got here"

      do i = 1, n
        nplus(i) = i+1
        nminus(i) = i-1
      end do 
      nplus(n) = 1
      nminus(1) = n
      do i = 1, m
        mplus(i) = i+1
        mminus(i) = i-1
      end do 
      mplus(m) = 1
      mminus(1) = m

      B = 0.0

      if (row_decomp) then
        do j = 1, n
          do i = 2, m-1
            B(i, j) = A(mplus(i),nplus(j))+A(mplus(i),j)+A(mplus(i),nminus(j)) &
              + A(i,nplus(j))+A(i,nminus(j))+A(mminus(i),nplus(j)) &
              + A(mminus(i),j) + A(mminus(i),nminus(j))
          end do 
        end do 

        do j = 1, n
          do i = 2, m-1
            if (B(i, j) .eq. 3.0) then
              A(i, j) = 1.0
            else if (B(i,j) .eq. 2.0) then
              !  do nothing
            else
              A(i, j) = 0.0
            end if
          end do
        end do
      else 
        do j = 2, n-1
          do i = 1, m
            B(i, j) = A(mplus(i),nplus(j))+A(mplus(i),j)+A(mplus(i),nminus(j)) &
              + A(i,nplus(j))+A(i,nminus(j))+A(mminus(i),nplus(j)) &
              + A(mminus(i),j) + A(mminus(i),nminus(j))
          end do 
        end do 

        do j = 2, n-1
          do i = 1, m
            if (B(i, j) .eq. 3.0) then
              A(i, j) = 1.0
            else if (B(i,j) .eq. 2.0) then
              !  do nothing
            else
              A(i, j) = 0.0
            end if
          end do
        end do
      end if
    end subroutine update

    subroutine update_2d(A, m, n)
      ! Note that this needs to be upgraded for non square matrices. 
      implicit none
      integer, intent(in) :: m,n
      integer :: A(:, :)
      integer, dimension(m, n) :: B
      integer, dimension(n) :: nplus, nminus
      integer, dimension(m) :: mplus, mminus
      integer :: i, j

      do i = 1, n
        nplus(i) = i+1
        nminus(i) = i-1
      end do 
      nplus(n) = 1
      nminus(1) = n
      do i = 1, m
        mplus(i) = i+1
        mminus(i) = i-1
      end do 
      mplus(m) = 1
      mminus(1) = m

      B = 0.0

      do j = 2, n-1
        do i = 2, m-1
          B(i, j) = A(mplus(i),nplus(j))+A(mplus(i),j)+A(mplus(i),nminus(j)) &
            + A(i,nplus(j))+A(i,nminus(j))+A(mminus(i),nplus(j)) &
            + A(mminus(i),j) + A(mminus(i),nminus(j))
        end do 
      end do 

      do j = 2, n-1
        do i = 2, m-1
          if (B(i, j) .eq. 3.0) then
            A(i, j) = 1.0
          else if (B(i,j) .eq. 2.0) then
            !  do nothing
          else
            A(i, j) = 0.0
          end if
        end do
      end do
    end subroutine update_2d

    ! This routine might never be used. I like the idea of a MPI_BCast instead
    ! of parallel IO
    subroutine preadmat(A, m, id)
      ! Performs some sort of parallelized reading of the matrices
      implicit none
      integer :: A(:,:), i, j, m, ie, id

      do i = 0, num_procs-1
        if (id.eq.i) then
          do j = 1, num_tasks
            print *, i, j
            read(*, "("//trim(str(m))//"I1)") A(:,j)
            print *, A(:,j)
          end do 
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ie) 
      end do 
    end subroutine preadmat

    subroutine pprintmat(A, m, id)
      ! Performs some sort of parallelized write of the matrices
      implicit none
      integer :: A(:,:), id, m, i, j, ie
      do i = 0, num_procs-1
        if (id.eq.i) then
          do j = 1, num_tasks
            print "("//trim(str(m))//"(I1, ' '))", A(:,j)
          end do 
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ie) 
      end do 
    end subroutine pprintmat

    subroutine pwritemat(A, m, id, fn, stat)
      ! Performs some sort of parallelized write of the matrices
      implicit none
      integer :: A(:,:), id, m, fn, i, j, np, ie, stat(:)
      do i = 0, num_procs-1
        if (id.eq.i) then
          do j = 1, num_tasks
            !write(fn, "("//trim(str(m))//"(I1, ' '))") A(:,j)
            call MPI_FILE_WRITE(fn, A(:,j), m, MPI_INT, stat, ie)
          end do 
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ie) 
      end do 
    end subroutine pwritemat

    subroutine pupdate_bound_1d(A, m, id, np, tg)
      ! After updating individual matrices, need to update ghost cells (involves
      ! the other processors)
      implicit none
      integer :: A(:,:), id, np, m, tg, stat, ie, rtg

#ifdef HB 
      if (row_decomp) then
        rtg = tg
        ! top
        call MPI_ISEND(A(2,1),1,row_type,pmap(1),ftag(id,pmap(1),rtg,1),&
          MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(1,1),1,row_type,pmap(1),ftag(pmap(1),id,rtg,2),&
          MPI_COMM_WORLD, stat,ie)
        ! bottom
        call MPI_ISEND(A(num_tasks+1,1),1,row_type,pmap(2),ftag(id,pmap(2),rtg,2),&
          MPI_COMM_WORLD, stat,ie)
        CALL MPI_IRECV(A(num_tasks+2,1),1,row_type,pmap(2),ftag(pmap(2),id,rtg,1),&
          MPI_COMM_WORLD, stat,ie)
      else 
        if (id.eq.0) then ! If root processor
          ! Prev
          CALL MPI_ISEND(A(:,2),m,MPI_INTEGER,np-1,tg+2*id+1,MPI_COMM_WORLD,stat,ie)
          CALL MPI_IRECV(A(:,1),m,MPI_INTEGER,np-1,tg+2*id+2,MPI_COMM_WORLD,stat,ie)
          ! Next 
          CALL MPI_ISEND(A(:,num_tasks+1),m,MPI_INTEGER,1,tg+2*id+4,MPI_COMM_WORLD,stat,ie)
          CALL MPI_IRECV(A(:,num_tasks+2),m, MPI_INTEGER,1,tg+2*id+3,MPI_COMM_WORLD,stat,ie)
        else if (id.eq.np-1) then ! If last processor
          ! Prev
          CALL MPI_ISEND(A(:,2),m,MPI_INTEGER,id-1,tg+2*id+1,MPI_COMM_WORLD,stat,ie)
          CALL MPI_IRECV(A(:,1),m,MPI_INTEGER,id-1,tg+2*id+2,MPI_COMM_WORLD,stat,ie)
          ! Next
          CALL MPI_ISEND(A(:,num_tasks+1),m,MPI_INTEGER,0,tg+2,MPI_COMM_WORLD,stat,ie)
          CALL MPI_IRECV(A(:,num_tasks+2),m, MPI_INTEGER,0,tg+1,MPI_COMM_WORLD,stat,ie)
        else  ! Any processor in the middle
          ! Prev
          CALL MPI_ISEND(A(:,2),m,MPI_INTEGER,id-1,tg+2*id+1,MPI_COMM_WORLD,stat,ie)
          CALL MPI_IRECV(A(:,1),m,MPI_INTEGER,id-1,tg+2*id+2,MPI_COMM_WORLD,stat,ie)
          ! Next
          call MPI_ISEND(A(:,num_tasks+1),m,MPI_INTEGER,id+1,tg+2*id+4,MPI_COMM_WORLD,stat,ie)
          call MPI_IRECV(A(:,num_tasks+2),m, MPI_INTEGER,id+1,tg+2*id+3,MPI_COMM_WORLD,stat,ie)
        end if
      end if
#else
      if (row_decomp) then
         ! top
        call MPI_ISEND(A(2,1),1,row_type,pmap(1),ftag(id,pmap(1),rtg,1),&
          MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(1,1),1,row_type,pmap(1),ftag(pmap(1),id,rtg,2),&
          MPI_COMM_WORLD, stat,ie)
        ! bottom
        call MPI_ISEND(A(num_tasks+1,1),1,row_type,pmap(2),ftag(id,pmap(2),rtg,2),&
          MPI_COMM_WORLD, stat,ie)
        CALL MPI_IRECV(A(num_tasks+2,1),1,row_type,pmap(1),ftag(pmap(1),id,rtg,2),&
          MPI_COMM_WORLD, stat,ie)
      else 
        if (id.eq.0) then ! If root processor
          ! Prev
          CALL MPI_ISENDRECV(A(:,2),m,MPI_INTEGER,np-1,tg+2*id+1,A(:,1),m,&
            MPI_INTEGER,np-1,tg+2*id+2,MPI_COMM_WORLD,stat,ie)
          ! Next 
          CALL MPI_ISENDRECV(A(:,num_tasks+1),m,MPI_INTEGER,1,tg+2*id+4,&
            A(:,num_tasks+2),m, MPI_INTEGER,1,tg+2*id+3,MPI_COMM_WORLD,stat,ie)
        else if (id.eq.np-1) then ! If last processor
          ! Prev
          CALL MPI_ISENDRECV(A(:,2),m,MPI_INTEGER,id-1,tg+2*id+1,A(:,1),m,&
            MPI_INTEGER,id-1,tg+2*id+2,MPI_COMM_WORLD,stat,ie)
          ! Next
          CALL MPI_ISENDRECV(A(:,num_tasks+1),m,MPI_INTEGER,0,tg+2,&
            A(:,num_tasks+2),m, MPI_INTEGER,0,tg+1,MPI_COMM_WORLD,stat,ie)
        else  ! Any processor in the middle
          ! Prev
          CALL MPI_ISENDRECV(A(:,2),m,MPI_INTEGER,id-1,tg+2*id+1,A(:,1),m,&
            MPI_INTEGER,id-1,tg+2*id+2,MPI_COMM_WORLD,stat,ie)
          ! Next
          CALL MPI_ISENDRECV(A(:,num_tasks+1),m,MPI_INTEGER,id+1,tg+2*id+4,&
            A(:,num_tasks+2),m, MPI_INTEGER,id+1,tg+2*id+3,MPI_COMM_WORLD,stat,ie)
        end if
      end if
#endif
    end subroutine pupdate_bound_1d

    subroutine pupdate_bound_2d(A, id, tg)
      ! After updating individual matrices, need to update ghost cells (involves
      ! the other processors)
      implicit none
      integer :: A(:,:), id, tg, stat, ie
      integer, dimension(ntn) :: helper, helper2
      integer, dimension(ntm) :: helper3, helper4

      ! pmap is length 8 array for each processor and contains the information
      ! as to which processor to ask for which information. 
      !        ---------------------
      ! pmap = |L|R|U|D|UL|UR|DL|DR|
      !        ---------------------
      ! Starting communication channels
      !
#ifdef HB
      ! Edges
        ! Left
        CALL MPI_ISEND(A(2:ntm+1,2),ntm,MPI_INTEGER,pmap(1),&
          ftag(id,pmap(1),tg,1),MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(2:ntm+1,1),ntm,MPI_INTEGER,pmap(1),&
          ftag(pmap(1),id,tg,2),MPI_COMM_WORLD,stat,ie)
        ! right
        CALL MPI_ISEND(A(2:ntm+1,ntn+1),ntm,MPI_INTEGER,pmap(2),&
          ftag(id,pmap(2),tg,2),MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(2:ntm+1,ntn+2),ntm, MPI_INTEGER,pmap(2),&
          ftag(pmap(2),id,tg,1),MPI_COMM_WORLD,stat,ie)
        ! Top
        CALL MPI_ISEND(A(2,2),1,row_type,pmap(3),&
          ftag(id,pmap(3),tg,3),MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(1,2),1,row_type,pmap(3),&
          ftag(pmap(3),id,tg,4),MPI_COMM_WORLD,stat,ie)
        ! Bottom
        CALL MPI_ISEND(A(ntm+1,2),1,row_type,pmap(4),&
          ftag(id,pmap(4),tg,4),MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(ntm+2,2),1,row_type,pmap(4),&
          ftag(pmap(4),id,tg,3),MPI_COMM_WORLD,stat,ie)
      ! Corners
        ! Upper Left
        CALL MPI_ISEND(A(2,2),1,MPI_INTEGER,pmap(5),&
          ftag(id,pmap(5),tg,5),MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(1,1),1,MPI_INTEGER,pmap(5),&
          ftag(pmap(5),id,tg,8),MPI_COMM_WORLD,stat,ie)
        ! Upper Right
        CALL MPI_ISEND(A(2,ntn+1),1,MPI_INTEGER,pmap(6),&
          ftag(id,pmap(6),tg,6),MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(1,ntn+2),1,MPI_INTEGER,pmap(6),&
          ftag(pmap(6),id,tg,7),MPI_COMM_WORLD,stat,ie)
        ! Down Left
        CALL MPI_ISEND(A(ntm+1,2),1,MPI_INTEGER,pmap(7),&
          ftag(id,pmap(7),tg,7),MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(ntm+2,1),1,MPI_INTEGER,pmap(7),&
          ftag(pmap(7),id,tg,6),MPI_COMM_WORLD,stat,ie)
        ! Down Right
        CALL MPI_ISEND(A(ntm+1,ntn+1),1,MPI_INTEGER,pmap(8),&
          ftag(id,pmap(8),tg,8),MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(ntm+2,ntn+2),1,MPI_INTEGER,pmap(8),&
          ftag(pmap(8),id,tg,5),MPI_COMM_WORLD,stat,ie)
      ! Done
#else
      ! Edges
        ! Left
        !                sendbuf counts type     dest       tag  recvbuf counts 
        CALL MPI_ISENDRECV(A(2:ntm+1,2),ntm,MPI_INTEGER,pmap(1),&
          ftag(id,pmap(1),tg,1),A(2:ntm+1,1),ntm,MPI_INTEGER,pmap(1),&
          ftag(pmap(1),id,tg,2),MPI_COMM_WORLD,stat,ie)
        !    type      src     tag       comm      stat   ie
        ! Right
        CALL MPI_ISENDRECV(A(2:ntm+1,ntn+1),ntm,MPI_INTEGER,pmap(2),&
          ftag(id,pmap(2),tg,2),A(2:ntm+1,ntn+2),ntm, MPI_INTEGER,pmap(2),&
          ftag(pmap(2),id,tg,1),MPI_COMM_WORLD,stat,ie)
        ! Top
        CALL MPI_ISENDRECV(A(2,2),1,row_type,pmap(3),&
          ftag(id,pmap(3),tg,3),A(1,2),1,row_type,pmap(3),&
          ftag(pmap(3),id,tg,4),MPI_COMM_WORLD,stat,ie)

        ! Bottom
        CALL MPI_ISENDRECV(A(ntm+1,2),1,row_type,pmap(4),&
          ftag(id,pmap(4),tg,4),A(ntm+2,2),1,row_type,pmap(4),&
          ftag(pmap(4),id,tg,3),MPI_COMM_WORLD,stat,ie)
      ! Corners
        ! Upper Left
        CALL MPI_ISENDRECV(A(2,2),1,MPI_INTEGER,pmap(5),&
          ftag(id,pmap(5),tg,5),A(1,1),1,MPI_INTEGER,pmap(5),&
          ftag(pmap(5),id,tg,8),MPI_COMM_WORLD,stat,ie)
        ! Upper Right
        CALL MPI_ISENDRECV(A(2,ntn+1),1,MPI_INTEGER,pmap(6),&
          ftag(id,pmap(6),tg,6),A(1,ntn+2),1,MPI_INTEGER,pmap(6),&
          ftag(pmap(6),id,tg,7),MPI_COMM_WORLD,stat,ie)
        ! Down Left
        CALL MPI_ISENDRECV(A(ntm+1,2),1,MPI_INTEGER,pmap(7),&
          ftag(id,pmap(7),tg,7),A(ntm+2,1),1,MPI_INTEGER,pmap(7),&
          ftag(pmap(7),id,tg,6),MPI_COMM_WORLD,stat,ie)
        ! Down Right
        CALL MPI_ISENDRECV(A(ntm+1,ntn+1),1,MPI_INTEGER,pmap(8),&
          ftag(id,pmap(8),tg,8),A(ntm+2,ntn+2),1,MPI_INTEGER,pmap(8),&
          ftag(pmap(8),id,tg,5),MPI_COMM_WORLD,stat,ie)
      ! Done
#endif
    end subroutine pupdate_bound_2d

    subroutine mpi_decomp_1d(char1,id, np, m,n, counts)
      ! Decomposes the data given to the program along columns (pencils)
      ! according to the number of CPU's
      implicit none
      integer :: id, np, n, extras, i, counts(:),ie,m
      character :: char1

      if (char1 .eq. 'R') then
        ! row wise decomposition
        row_decomp = .true.
        num_procs = np
        num_tasks = m/np 
        extras = mod(m, np)
        counts = num_tasks
        if (extras .ne. 0) then
          do i = 1, extras
            if (id .eq. np-i) then
              num_tasks = num_tasks+1
            end if
            counts(np-i+1) = num_tasks+1
          end do 
        end if
        if (id .eq. 0) then
          pmap(1) = np-1
          pmap(2) = id+1
        else if (id.eq.np-1) then
          pmap(1) = id-1
          pmap(2) = 0
        else
          pmap(1) = id-1
          pmap(2) = id+1
        end if
        ! create MPI row_type
        call MPI_TYPE_VECTOR(n,1,num_tasks+2,MPI_INTEGER,row_type,ie)
        call MPI_TYPE_COMMIT(row_type,ie) 
      else 
        ! column decomposition
        num_procs = np
        num_tasks = n/np
        extras = mod(n, np) 
        counts = num_tasks
        if (extras .ne. 0) then
          do i = 1, extras
            if(id.eq.np-i) then
              num_tasks = num_tasks+1 
            end if
            counts(np-i+1) = num_tasks+1
          end do
        end if
      end if
    end subroutine mpi_decomp_1d

    subroutine mpi_decomp_2d(id, np, m, n, counts, dv)
      ! Decomposes the data given to the program into subgrids according to the
      ! number of CPU's
      implicit none
      integer :: id, np, m, n, ex_m, ex_n, counts(:,:), dv(:,:)
      integer :: i, j, k, ie
      integer, allocatable :: ind_array(:)

      num_procs = np
      ! pmap is length 8 array for each processor and contains the information
      ! as to which processor to ask for which information. 
      !        ---------------------
      ! pmap = |L|R|U|D|UL|UR|DL|DR|
      !        ---------------------

      ! Cases for Domain Decomposition
      if (int(sqrt(real(np)))**2 .eq. np) then 
        ! if np is a perfect square then break domain into "squares"

        ! i.e. if there are 9 processors, the domain will be decomposed into the
        ! following grid, labeled by processor
        ! 
        ! -------
        ! |1|4|7| 
        ! -------
        ! |2|5|8| 
        ! -------
        ! |3|6|9| 
        ! -------

        ! these two are the default values for the subarray
        nm =  int(sqrt(real(np)))
        nn =  nm

        ! Mapping -----------------------------------------------
          ! Left Mappin
          if(id.lt.nm) then
            pmap(1) = id+nm*(nn-1) 
            ! UL, DL
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(5) = id+nm*(nn)-1
              pmap(7) = id+nm*(nn-1)+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(5) = id+nm*(nn-1)-1
              pmap(7) = id+nm*(nn-2)+1
            ! else middle
            else 
              pmap(5) = id+nm*(nn-1)-1
              pmap(7) = id+nm*(nn-1)+1
            end if
          else
            pmap(1) = id-nm
            ! UL, DL
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(5) = id-1
              pmap(7) = id-nm+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(5) = id-nm-1
              pmap(7) = id-2*nm+1
            ! else middle
            else 
              pmap(5) = id-nm-1
              pmap(7) = id-nm+1 
            end if
          end if

          ! Right Mapping
          if(id.ge.nm*(nn-1)) then
            pmap(2) = id-nm*(nn-1) 
            ! UR, DR
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(6) = id-nm*(nn-2)-1
              pmap(8) = id-nm*(nn-1)+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(6) = id-nm*(nn-1)-1
              pmap(8) = id-nm*(nn)+1
            ! else middle
            else 
              pmap(6) = id-nm*(nn-1)-1
              pmap(8) = id-nm*(nn-1)+1
            end if
          else
            pmap(2) = id+nm
            ! UR, DR
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(6) = id+2*nm-1
              pmap(8) = id+nm+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(6) = id+nm-1
              pmap(8) = id+1
            ! else middle
            else 
              pmap(6) = id+nm-1
              pmap(8) = id+nm+1
            end if
          end if

          ! Up Mapping
          if(mod(id,nm).eq.0) then
            pmap(3) = id+nm-1
          else
            pmap(3) = id-1
          end if

          ! Down Mapping
          if(mod(id+1,nm).eq.0) then
            pmap(4) = id-nm+1
          else
            pmap(4) = id+1
          end if
        ! Mapping -----------------------------------------------

        allocate(ind_array(nm))
        do i = 1, nm
          ind_array(i) = 1 + (i-1)*nm
        end do 

        ! Agglomeration ----------------------------------------
          counts(:,1) = m/nm
          counts(:,2) = n/nn
          ex_m = mod(m, nm) 
          ex_n = mod(n, nn) 

          if (ex_m .ne. 0) then
            do i = 1, ex_m
              counts((nm-i)+ind_array,1) = counts((nm-i)+ind_array,1) + 1
            end do
          end if
          if (ex_n .ne. 0) then
            do i = 1, ex_n
              counts(np-nm*i+1:np-nm*(i-1),2)=counts(np-nm*i+1:np-nm*(i-1),2)+1
            end do
          end if
          ntm = counts(id+1,1)
          ntn = counts(id+1,2)
          ! creating disp_vec
          do i = 1, nm
            do j = 1, nn
              k = nm*(j-1)+i
              if(i .eq. 1) then
                dv(k,1) = 0
              else
                dv(k,1) = dv(k-1,1)+counts(k-1,1)
              end if
              if(j .eq. 1) then
                dv(k,2) = 0
              else 
                dv(k,2) = dv(k-nm,2)+counts(k-nm,2)
              end if
            end do 
          end do 
        ! Agglomeration ----------------------------------------
    
        deallocate(ind_array)

        ! Create MPI types for row vector sends
        call MPI_TYPE_VECTOR(ntn,1,ntm+2,MPI_INTEGER,row_type,ie)
        call MPI_TYPE_COMMIT(row_type,ie)
        sm = dv(id+1,1)
        sn = dv(id+1,2)
        em = sm + counts(id+1,1)-1
        en = sn + counts(id+1,2)-1
        !print *, m,n,ntm,ntn,sm,sn
        !call MPI_TYPE_CREATE_SUBARRAY(2, (/m,n/), (/ntm,ntn/), (/sm,sn/), &
          !MPI_ORDER_FORTRAN, MPI_INTEGER, g2s_type, ie)
        !call MPI_TYPE_COMMIT(g2s_type,ie)
        !call MPI_TYPE_CREATE_SUBARRAY(2, (/ntm+2,ntn+2/), (/ntm,ntn/),(/2,2/),&
          !MPI_ORDER_FORTRAN, MPI_INTEGER, s2s_type,ie)
        !call MPI_TYPE_COMMIT(s2s_type,ie)
      else if (mod(np,2).eq.0) then
        ! if np is even, divide the domain in half, and then split rows

        ! -----
        ! |1|4| 
        ! -----
        ! |2|5| 
        ! -----
        ! |3|6| 
        ! -----

        nm = np/2
        nn = 2

        ! Mapping -----------------------------------------------
          ! Left Mappin
          if(id.lt.nm) then
            pmap(1) = id+nm*(nn-1) 
            ! UL, DL
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(5) = id+nm*(nn)-1
              pmap(7) = id+nm*(nn-1)+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(5) = id+nm*(nn-1)-1
              pmap(7) = id+nm*(nn-2)+1
            ! else middle
            else 
              pmap(5) = id+nm*(nn-1)-1
              pmap(7) = id+nm*(nn-1)+1
            end if
          else
            pmap(1) = id-nm
            ! UL, DL
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(5) = id-1
              pmap(7) = id-nm+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(5) = id-nm-1
              pmap(7) = id-2*nm+1
            ! else middle
            else 
              pmap(5) = id-nm-1
              pmap(7) = id-nm+1 
            end if
          end if

          ! Right Mapping
          if(id.ge.nm*(nn-1)) then
            pmap(2) = id-nm*(nn-1) 
            ! UR, DR
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(6) = id-nm*(nn-2)-1
              pmap(8) = id-nm*(nn-1)+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(6) = id-nm*(nn-1)-1
              pmap(8) = id-nm*(nn)+1
            ! else middle
            else 
              pmap(6) = id-nm*(nn-1)-1
              pmap(8) = id-nm*(nn-1)+1
            end if
          else
            pmap(2) = id+nm
            ! UR, DR
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(6) = id+2*nm-1
              pmap(8) = id+nm+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(6) = id+nm-1
              pmap(8) = id+1
            ! else middle
            else 
              pmap(6) = id+nm-1
              pmap(8) = id+nm+1
            end if
          end if

          ! Up Mapping
          if(mod(id,nm).eq.0) then
            pmap(3) = id+nm-1
          else
            pmap(3) = id-1
          end if

          ! Down Mapping
          if(mod(id+1,nm).eq.0) then
            pmap(4) = id-nm+1
          else
            pmap(4) = id+1
          end if
        ! Mapping -----------------------------------------------

        allocate(ind_array(nn))
        do i = 1, nn
          ind_array(i) = 1 + (i-1)*nm
        end do 

        ! Agglomeration ----------------------------------------
          counts(:,1) = m/nm
          counts(:,2) = n/nn
          ex_m = mod(m, nm) 
          ex_n = mod(n, nn) 

          if (ex_m .ne. 0) then
            do i = 1, ex_m
              counts((nm-i)+ind_array,1) = counts((nm-i)+ind_array,1) + 1
            end do
          end if
          if (ex_n .ne. 0) then
            do i = 1, ex_n
              counts(np-nm*i+1:np-nm*(i-1),2)=counts(np-nm*i+1:np-nm*(i-1),2)+1
            end do
          end if
          ntm = counts(id+1,1)
          ntn = counts(id+1,2)
          ! creating disp_vec
          do i = 1, nm
            do j = 1, nn
              k = nm*(j-1)+i
              if(i .eq. 1) then
                dv(k,1) = 0
              else
                dv(k,1) = dv(k-1,1)+counts(k-1,1)
              end if
              if(j .eq. 1) then
                dv(k,2) = 0
              else 
                dv(k,2) = dv(k-nm,2)+counts(k-nm,2)
              end if
            end do 
          end do 
        ! Agglomeration ----------------------------------------

        ! Mpi Derivated Datatypes for Row Comm
        call MPI_TYPE_VECTOR(ntn,1,ntm+2,MPI_INTEGER,row_type,ie)
        call MPI_TYPE_COMMIT(row_type,ie)

        deallocate(ind_array)
      else if (mod(np,3).eq.0) then
        ! if np is divisible by 3, then decompose the domain into 3 columns and
        ! then split those columns to rows. 

        ! -------
        ! |0|2|5|
        ! -------
        ! |1|3|6|
        ! -------

        nm = np/3
        nn = 3

        ! Mapping -----------------------------------------------
          ! Left Mappin
          if(id.lt.nm) then
            pmap(1) = id+nm*(nn-1) 
            ! UL, DL
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(5) = id+nm*(nn)-1
              pmap(7) = id+nm*(nn-1)+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(5) = id+nm*(nn-1)-1
              pmap(7) = id+nm*(nn-2)+1
            ! else middle
            else 
              pmap(5) = id+nm*(nn-1)-1
              pmap(7) = id+nm*(nn-1)+1
            end if
          else
            pmap(1) = id-nm
            ! UL, DL
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(5) = id-1
              pmap(7) = id-nm+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(5) = id-nm-1
              pmap(7) = id-2*nm+1
            ! else middle
            else 
              pmap(5) = id-nm-1
              pmap(7) = id-nm+1 
            end if
          end if

          ! Right Mapping
          if(id.ge.nm*(nn-1)) then
            pmap(2) = id-nm*(nn-1) 
            ! UR, DR
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(6) = id-nm*(nn-2)-1
              pmap(8) = id-nm*(nn-1)+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(6) = id-nm*(nn-1)-1
              pmap(8) = id-nm*(nn)+1
            ! else middle
            else 
              pmap(6) = id-nm*(nn-1)-1
              pmap(8) = id-nm*(nn-1)+1
            end if
          else
            pmap(2) = id+nm
            ! UR, DR
            ! if top row
            if(mod(id,nm).eq.0) then
              pmap(6) = id+2*nm-1
              pmap(8) = id+nm+1
            ! else if bottom row
            else if(mod(id+1,nm).eq.0) then
              pmap(6) = id+nm-1
              pmap(8) = id+1
            ! else middle
            else 
              pmap(6) = id+nm-1
              pmap(8) = id+nm+1
            end if
          end if

          ! Up Mapping
          if(mod(id,nm).eq.0) then
            pmap(3) = id+nm-1
          else
            pmap(3) = id-1
          end if

          ! Down Mapping
          if(mod(id+1,nm).eq.0) then
            pmap(4) = id-nm+1
          else
            pmap(4) = id+1
          end if
        ! Mapping -----------------------------------------------

        allocate(ind_array(nn))
        do i = 1, nn
          ind_array(i) = 1 + (i-1)*nm
        end do 

        ! Agglomeration ----------------------------------------
          counts(:,1) = m/nm
          counts(:,2) = n/nn
          ex_m = mod(m, nm) 
          ex_n = mod(n, nn) 

          if (ex_m .ne. 0) then
            do i = 1, ex_m
              counts((nm-i)+ind_array,1) = counts((nm-i)+ind_array,1) + 1
            end do
          end if
          if (ex_n .ne. 0) then
            do i = 1, ex_n
              counts(np-nm*i+1:np-nm*(i-1),2)=counts(np-nm*i+1:np-nm*(i-1),2)+1
            end do
          end if
          ntm = counts(id+1,1)
          ntn = counts(id+1,2)
          ! creating disp_vec
          do i = 1, nm
            do j = 1, nn
              k = nm*(j-1)+i
              if(i .eq. 1) then
                dv(k,1) = 0
              else
                dv(k,1) = dv(k-1,1)+counts(k-1,1)
              end if
              if(j .eq. 1) then
                dv(k,2) = 0
              else 
                dv(k,2) = dv(k-nm,2)+counts(k-nm,2)
              end if
            end do 
          end do 
        ! Agglomeration ----------------------------------------

        ! MPI Derivated Datatypes for Row Comm
        call MPI_TYPE_VECTOR(ntn,1,ntm+2,MPI_INTEGER,row_type,ie)
        call MPI_TYPE_COMMIT(row_type,ie)

        deallocate(ind_array)
      else
        ! if np is odd and doesn't satisfy any other criterium, then decompose
        ! into 2 columns, and then split those columnds evenly to rows, and for
        ! the last processor give it a couple full rows but with less rows than
        ! the column splot rows. 
         
        print *, " This decomposition method is not yet integrated yet. Be"//&
          " careful"

        ! NOTE: THIS DECOMPOSITION METHOD IS NOT YET INPLEMENTED BE AWARE IT IS
        ! NOT FUNCTIONAL YET

        !nm = 2
        !nn = (np-1)/2
         
        !ex_m = mod(m, nm)
        !ex_n = mod(n, nn)
        !if (ex_m .ne. 0) then

        !end if
        !if (ex_n .ne. 0) then
          ! give extras to the last processor

        !end if
      end if
    end subroutine mpi_decomp_2d

    function ftag(src,dest,step,dir) 
      ! Returns a tag for MPI Communication based on the following function
      ! tag = **** | ** | ** |  * 
      !       step | src|dst | dir
      ! tag = 1*10**(11) + src*10**8 + dest*10**5+step
      ! Assumptions: less than 100 processors, less than 100000 steps
      implicit none
      integer :: ftag, src, dest, step, dir
      ftag = step*10**5 + src*10**3 + dest*10+dir
    end function ftag

    function rtag(src,step) 
      ! Returns a tag for MPI Communication based on the following function
      ! tag = 1|***|***|******
      !         src|dst| step
      ! tag = 1*10**(11) + src*10**8 + dest*10**5+step
      ! Assumptions: less than 100 processors, less than 100000 steps
      implicit none
      integer :: rtag, src, step
      rtag = src + step*10**2
    end function rtag

    character(len=20) function str(k)
      ! "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
    end function str
end module pgameoflife


