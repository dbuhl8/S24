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
  
  ! Subroutines and Functions
  contains

    subroutine getprobsize(m,n)
      ! Reads m,n from stin file
      integer :: n, m
      read(*,"(2I6)") m, n
    end subroutine getprobsize

    subroutine readmat(A)
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
    subroutine preadmat(A, n, np)
      ! Performs some sort of parallelized reading of the matrices
      implicit none
      integer :: A(:,:), np, i, j, n, ie

      do i = 1, np
        do j = 1, num_tasks
          read(*, "("//trim(str(n))//"I1)") A(:,j)
        end do 
        call MPI_BARRIER(MPI_COMM_WORLD,ie) 
      end do 
    end subroutine preadmat

    subroutine pwritemat(A, n, id, fn)
      ! Performs some sort of parallelized write of the matrices
      implicit none
      integer :: A(:,:), id, n, fn, i, j, np, ie
      do i = 1, num_procs
        if (id.eq.i) then
          do j = 1, num_tasks
            write(fn, "("//trim(str(n))//"I1)") A(:,j)
          end do 
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ie) 
      end do 
    end subroutine pwritemat

    subroutine pupdate_bound_1d(A, n, id, np, tg)
      ! After updating individual matrices, need to update ghost cells (involves
      ! the other processors)
      implicit none
      integer :: A(:,:), id, np, n, tg, stat, ie

      if (id.eq.0) then
        CALL MPI_ISEND(A(:,2),n,MPI_INTEGER,np-1,tg,MPI_COMM_WORLD,stat,ie)
        CALL MPI_ISEND(A(:,num_tasks+1),n,MPI_INTEGER,1,tg,MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(:,1), n, MPI_INTEGER,np-1,tg,MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(:,num_tasks+2),n,MPI_INTEGER,1,tg,MPI_COMM_WORLD,stat,ie)
      else if (id.eq.np-1) then
        CALL MPI_ISEND(A(:,2),n,MPI_INTEGER,np-2,tg,MPI_COMM_WORLD,stat,ie)
        CALL MPI_ISEND(A(:,num_tasks+1),n,MPI_INTEGER,0,tg,MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(:,1), n, MPI_INTEGER,np-2,tg, MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(:,num_tasks+2),n,MPI_INTEGER,0,tg,MPI_COMM_WORLD,stat,ie)
      else 
        CALL MPI_ISEND(A(:,2),n,MPI_INTEGER,id-1,tg,MPI_COMM_WORLD,stat,ie)
        CALL MPI_ISEND(A(:,num_tasks+1),n,MPI_INTEGER,id+1,tg,MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(:,1),n,MPI_INTEGER,id-1,tg,MPI_COMM_WORLD,stat,ie)
        CALL MPI_IRECV(A(:,num_tasks+2),n,MPI_INTEGER,id+1,tg,MPI_COMM_WORLD,stat,ie)
      end if
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
          counts(np-i) = num_tasks+1
        end do
      end if
    end subroutine mpi_decomp_1d

    subroutine mpi_decomp_2d(id, np, m, n)
      ! Decomposes the data given to the program into subgrids according to the
      ! number of CPU's
      implicit none
      integer :: id, np, m, n
    end subroutine mpi_decomp_2d

    character(len=20) function str(k)
      ! "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
    end function str
end module pgameoflife


