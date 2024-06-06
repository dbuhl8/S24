!File: gameoflife.f90
!Author: Dante Buhl
!Dependencies:
!Description: A module containing all of the source code to run the game of life

module pgameoflife
  ! Divided into serial and parallel versions
  use MPI
  implicit none

  ! General Utility Vars
  integer :: sx, ex, sy, ey, num_tasks, num_procs
  integer :: nm, nn, npm, npn
  
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

      do j = 1, n
        do i = 1, m
          B(i, j) = A(mplus(i),nplus(j))+A(mplus(i),j)+A(mplus(i),nminus(j)) &
            + A(i,nplus(j))+A(i,nminus(j))+A(mminus(i),nplus(j)) &
            + A(mminus(i),j) + A(mminus(i),nminus(j))
        end do 
      end do 

      do j = 1, n
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
    end subroutine update

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
      integer :: A(:,:), id, np, m, tg, stat, ie

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
      !call MPI_BARRIER(MPI_COMM_WORLD, ie)
    end subroutine pupdate_bound_1d

    subroutine mpi_decomp_1d(id, np, n, counts)
      ! Decomposes the data given to the program along columns (pencils)
      ! according to the number of CPU's
      implicit none
      integer :: id, np, n, extras, i, counts(:)

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
    end subroutine mpi_decomp_1d

    subroutine mpi_decomp_2d(id, np, m, n, counts)
      ! Decomposes the data given to the program into subgrids according to the
      ! number of CPU's
      implicit none
      integer :: id, np, m, n, ex_m, ex_n, counts(:,:), i
      integer, allocatable :: ind_array(:)

      num_procs = np

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

        allocate(ind_array(nm))
        do i = 1, nm
          ind_array(i) = 1 + (i-1)*nm
        end do 

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
    
        !nm = counts(id+1, 1)
        !nn = counts(id+1, 2)
        deallocate(ind_array)
      
      else if (mod(np,2).eq.0) then
        ! if np is even, divide the domain in half, and then split rows
        nm = 2
        nn = np/2
         
        ex_m = mod(m, nm)
        ex_n = mod(n, nn)
        if (ex_m .ne. 0) then

        end if
        if (ex_n .ne. 0) then

        end if
         
      else if (mod(np,3).eq.0) then
        ! if np is divisible by 3, then decompose the domain into 3 columns and
        ! then split those columns to rows. 
        nm = 3
        nn = np/3
         
        ex_m = mod(m, nm)
        ex_n = mod(n, nn)
        if (ex_m .ne. 0) then

        end if
        if (ex_n .ne. 0) then

        end if

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

    character(len=20) function str(k)
      ! "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
    end function str
end module pgameoflife


