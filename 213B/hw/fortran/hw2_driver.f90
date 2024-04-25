!File: hw2_driver.f90
!Author: Dante Buhl
!Requirements: NumDE.mod, LinAl.mod

program hw2

use NumDE
use LinAl

implicit none

  integer, parameter :: ma=3, na=3
  real, dimension(ma, na) :: A
  real :: tol=10.d-8
  logical :: shift=.true.

  A(1, :) = (/0.0, 10.0, -10.0/)
  A(2, :) = (/-100.0, -1.0, 0.0/)
  A(3, :) = (/0.0, 10.0, -100.0/)
 
  !Note: This fails becaause the eigenvalues of this matrix are complex 
  ! need to upgrade my linal module to handle complex matrices/eigenvectors.
  !call eigQR(A, ma, shift, tol) 
  ! Note that the eigenvalues of A are:
  ! -0.9666 + i 30.1255
  ! -0.9666 - i 30.1255
  ! -99.0667 

  !print *, " "
  !print *, " Eigenvalues of A: "
  !print *, " ---------------- "
  !print *, " "
  !call printmat(A, ma, na)
  !print *, " "

end program hw2
