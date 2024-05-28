!File: driver.f90
!Author: Dante Buhl
!Dependencies: gameoflife
!Description: A driver routine to run the game of life

! NOTE: Needs an input file giving some initial condition

program driver

  use gameoflife
  
  implicit none

  ! General Util Vars
  integer :: i, j, k,fn=10
  
  ! Game of Life Vars
  integer, allocatable :: A(:,:), G(:,:)
  integer :: m, n, nt
  
 
  ! simulate the game of life
  call readmat(A, m, n)
   
  open(fn, file='gof.dat')

  ! ------------------------- Sequential Component
  call writemat(A, m, n, fn)
  nt = 1000
  do i = 1, nt
    call update(A, m, n)
    call writemat(A, m, n,fn)
  end do 
  ! ------------------------- END Sequential Component

  close(fn)

  open(fn+1, file='params.dat')
    write(fn+1, "(3I6)") m, n, nt+1
  close(fn+1)

  deallocate(A)
end program driver


