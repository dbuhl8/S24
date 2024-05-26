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
    write(fn+1, "(3I6)") n, m, nt+1
  close(fn+1)

  deallocate(A)
end program driver


