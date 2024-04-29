!File: hw2_driver.f90
!Author: Dante Buhl
!Requirements: NumDE.mod, LinAl.mod

program hw2
  use NumDE
  use LinAl

  implicit none
  !problem 1.a
  integer, parameter :: ma=3, na=3
  real, dimension(ma, na) :: A
  real :: tol=10.d-8
  logical :: shift=.true.

  !problem 1.c
  real :: dt, tstart=0.0, tstop=10.0, dt_pert=10.d-4
  real :: dt_crit=0.0055
  real, allocatable :: Y(:,:),T(:),Y2(:,:),T2(:),Y3(:,:),T3(:)
  integer :: num_points, i
  integer,dimension(3) :: np_array
  

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

  !Run AB3 for problem 1 and return a plot.
  !very dt to be just above and below the threshold dt for absolute stability.
  dt = 10.d-4
  num_points = int((tstop-tstart)/dt)
  np_array(1) = num_points
  allocate(Y(num_points+1,3), T(num_points+1))
  Y(1,:) = (/10., 10., 10./)
  T(1) = 0.
  call AB3(A, Y, T, ma, num_points, dt)
  dt = dt_crit - dt_pert
  num_points = int((tstop-tstart)/dt)
  np_array(2) = num_points
  allocate(Y2(num_points+1,3), T2(num_points+1))
  Y2(1,:) = (/10., 10., 10./)
  T2(1) = 0.
  call AB3(A, Y2, T2, ma, num_points, dt)
  dt = dt_crit + dt_pert
  num_points = int((tstop-tstart)/dt)
  np_array(3) = num_points
  allocate(Y3(num_points+1,3), T3(num_points+1))
  Y3(1,:) = (/10., 10., 10./)
  T3(1) = 0.
  call AB3(A, Y3, T3, ma, num_points, dt)   

  open(12, file='AB3Traj.dat')
    do i = 1, np_array(3)
      write (12, "(3(F8.5, 3F16.8))") T(i),Y(i,:),T2(i),Y2(i,:),T3(i),Y3(i,:) 
    end do
    if(np_array(3) .ne. np_array(2)) then
      do i = np_array(3)+1, np_array(2)
        write (12, "(2(F8.5, 3F16.8))") T(i),Y(i,:),T2(i),Y2(i,:)
      end do 
      do i = np_array(2)+1, np_array(1)
        write (12, "(F8.5, 3F16.8)") T(i),Y(i,:)
      end do 
    else 
      do i = np_array(3)+1, np_array(1)
        write (12, "(F8.5, 3F16.8)") T(i), Y(i,:)
      end do 
    end if 
  close(12)
  deallocate(Y, T, Y2, T2, Y3, T3)
end program hw2
