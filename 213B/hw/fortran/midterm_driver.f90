!File: midterm_driver.f90
!Author: Dante Buhl
!Dependencies: NumDE.mod LinAl.mod

program midterm_driver

  use NumDE
  use LinAl

  implicit none

  integer, parameter :: num_points=1000, ma=4
  real, dimension(ma, ma) :: F
  real, allocatable :: Y(:, :), Y2(:, :)
  real, allocatable :: T(:), T2(:)
  real :: dt=0.3, dt_crit = 0.3387, dt_pert = 0.0001
  integer :: np1, np2,i
  real :: t_fin=10.0

  F(1, :) = (/-1.0, 3.0,-5.0, 7.0 /)     
  F(2, :) = (/ 0.0,-2.0, 4.0,-6.0 /)     
  F(3, :) = (/ 0.0, 0.0,-4.0, 6.0 /)     
  F(4, :) = (/ 0.0, 0.0, 0.0,-16.0/)     

  np1 = int(t_fin/dt_crit)
  np2 = int(t_fin/(dt_crit+dt_pert))

  allocate(Y(np1+1, ma), T(np1+1), Y2(np2+1, ma), T2(np2+1))

  Y(1, :) = (/ 1.0, 1.0, 1.0, 1.0 /)
  Y2(1, :) = (/ 1.0, 1.0, 1.0, 1.0 /)

  T(1) = 0.0
  T2(1) = 0.0

  call MidtermRK3(F, Y, T, ma, np1, dt_crit)
  call MidtermRK3(F, Y2, T2, ma, np2, dt_crit+dt_pert)
 
  open(15, file="RK3.dat")
  if (np2 .gt. np1) then 
    do i = 1, np1
      write(15, "(10(F12.8, '    '))") T2(i), Y2(i, :), T(i), Y(i,:)
    end do 
    do i = np1+1, np2
      write(15, "(5(F12.8, '    '))") T2(i), Y2(i, :)
    end do 
  else 
    do i = 1, np2
      write(15, "(10(F12.8, '    '))") T2(i), Y2(i, :), T(i), Y(i,:)
    end do 
  end if
  close(15)

end program midterm_driver
