!File: hw2_driver.f90
!Author: Dante Buhl
!Requirements: NumDE.mod, LinAl.mod

program hw2

use NumDE
use LinAl

implicit none

  integer, parameter :: ma=3, na=3
  real, dimension(ma, na) :: A
  real :: tol=10.d-14
  logical :: shift=.false.

  A(1, :) = (/0.0, 10.0, -10.0/)
  A(2, :) = (/-100.0, -1.0, 0.0/)
  A(3, :) = (/0.0, 10.0, -100.0/)
  
  call eigQR(A, ma, shift, tol) 

  print *, " "
  print *, " Eigenvalues of A: "
  print *, " ---------------- "
  print *, " "
  call printmat(A, ma, na)
  print *, " "

end program hw2
