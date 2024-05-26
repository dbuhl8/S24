!File: gameoflife.f90
!Author: Dante Buhl
!Dependencies:
!Description: A module containing all of the source code to run the game of life

module gameoflife
  ! Divided into serial and parallel versions
  implicit none

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

    character(len=20) function str(k)
      ! "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
    end function str
end module gameoflife


