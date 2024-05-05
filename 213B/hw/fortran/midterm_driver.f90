!File: midterm_driver.f90
!Author: Dante Buhl
!Dependencies: NumDE.mod LinAl.mod

program midterm_driver

  use NumDE
  use LinAl

  implicit none

  integer, parameter :: num_points=1000, ma=4
  real, dimension(ma, ma) :: F
  real, dimension(nk, ma) :: Y
  real, dimension(nk) :: T
  real :: dt=10.d-4

  F(1, :) = (/-1.0, 3.0,-5.0, 7.0 /)     
  F(2, :) = (/ 0.0,-2.0, 4.0,-6.0 /)     
  F(3, :) = (/ 0.0, 0.0,-4.0, 6.0 /)     
  F(4, :) = (/ 0.0, 0.0, 0.0,-16.0/)     

  Y(1, :) = (/ 1.0, 1.0, 1.0, 1.0 /)

  T(0) = 0.0

  call MidtermRK3(F, Y, T, ma, num_points, dt)

end program midterm_driver
