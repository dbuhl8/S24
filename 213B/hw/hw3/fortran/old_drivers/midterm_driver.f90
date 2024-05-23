!File: midterm_driver.f90
!Author: Dante Buhl
!Dependencies: NumDE.mod LinAl.mod

program midterm_driver

  use NumDE
  use LinAl

  implicit none

  !Q1
  integer, parameter :: np_rk4=60000, m_rk4=4
  real, dimension(m_rk4, np_rk4+1) :: Y_rk4, Y2_rk4 
  real, dimension(np_rk4+1) ::  Y_act_rk4, T_rk4, Err_rk4, DY_act_rk4
  real :: dt_rk4, tol = 10.d-10, error_tot, prev_err
  real, dimension(3) :: v1, v2
  real, dimension(2) :: y_error, dy_error
  real, dimension(m_rk4, m_rk4) :: F_rk4
  real :: dv1, dv2, de1, de2

  !Q2
  integer, parameter :: num_points=1000, ma=4
  real, dimension(ma, ma) :: F, eye
  real, allocatable :: Y(:, :), Y2(:, :), Y3(:,:)
  real, allocatable :: T(:), T2(:), T3(:)
  real :: dt=0.001, dt_crit = 0.3387, dt_pert = 0.0001
  integer :: np1, np2,i, np3
  real :: t_fin=100.0

  !Q3 
  integer :: m3=3
  real, dimension(3, 3) :: A
  real, dimension(3, 1) :: B, X
  integer, dimension(3) :: P3
  logical :: bool=.FALSE.


  print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 1: BVP with Shooting Method "
    print *, " "
    v1 = 0.0
    v2 = 0.0
    v1(2) = 1.0
    v2(2) = 1.0
    Y_rk4 = 0.0
    Y2_rk4 = 0.0
    dt_rk4 = 1./np_rk4

    T_rk4 = 0

    F_rk4(1,:) = (/ 0.0, 1.0, 0.0, 0.0/)
    F_rk4(2,:) = (/ 0.0, 0.0, 1.0, 0.0/)
    F_rk4(3,:) = (/ 0.0, 0.0, 0.0, 1.0/)
    F_rk4(4,:) = (/ 0.0, 0.0, 0.0, 0.0/)

    Y_rk4(1, 1) = 0.0 
    Y_rk4(2, 1) = 0.0
    Y_rk4(3, 1) = 1.0
    Y_rk4(4, 1) = 1.0
    Y2_rk4(1, 1) = 0.0
    Y2_rk4(2, 1) = 0.0
    Y2_rk4(3, 1) = 0.0
    Y2_rk4(4, 1) = 0.0

    call MidtermRK4(F_rk4, Y_rk4, T_rk4, m_rk4, np_rk4, dt_rk4)

    ! compute analytical sol
    Y_act_rk4 = (T_rk4**6)/360. - (T_rk4**3)/90. + (T_rk4**2)/120.
    DY_act_rk4 = (T_rk4**5)/60. - (T_rk4**2)/30. + (T_rk4)/60.
    y_error(2) =  - Y_rk4(1, np_rk4+1) 
    dy_error(2) = - Y_rk4(2, np_rk4+1) 
    T_rk4 = 0

    call MidtermRK4(F_rk4, Y2_rk4, T_rk4, m_rk4, np_rk4, dt_rk4)

    y_error(1) = - Y2_rk4(1, np_rk4+1) 
    dy_error(1) = - Y2_rk4(2, np_rk4+1) 

    error_tot = sqrt(dy_error(2)**2+y_error(2)**2) !+ abs(y_error(2))
    !print *, error_tot, y_error(2), dy_error(2)

    do while (error_tot > tol)
       !update v
      dv1 = v1(2) - v1(1)
      dv2 = v2(2) - v2(1)
      de1 = y_error(2) - y_error(1)
      de2 = dy_error(2) - dy_error(1)
      v1(3) = v1(2) - y_error(2)*(dv1/de1) !- dy_error(2)*(dv1/de2)
      v2(3) = v2(2) - dy_error(2)*(dv2/de2) !- y_error(2)*(dv2/de1)
      print *, "V1: ", v1(3), "V2: ", v2(3)

      ! perform rk4 again
      T_rk4 = 0
      !Y_rk4(:,1) = (/0.0,0.0,v2(3),v1(3)/)
      Y_rk4(1, 1) = 0.0
      Y_rk4(2, 1) = 0.0
      Y_rk4(3, 1) = v2(3)
      Y_rk4(4, 1) = v1(3)
      call MidtermRK4(F_rk4, Y_rk4, T_rk4, m_rk4, np_rk4, dt_rk4)
      v1(1:2) = v1(2:3)
      v2(1:2) = v2(2:3)

      ! compute the error
      y_error(1) = y_error(2)
      y_error(2) =  - Y_rk4(1, np_rk4+1) 
      dy_error(1) = dy_error(2)
      dy_error(2) = - Y_rk4(2, np_rk4+1)

      ! update error
      error_tot = sqrt(dy_error(2)**2+y_error(2)**2)
      !print *, error_tot, y_error(2), dy_error(2)
    end do 
    print *, "v1: ", v1
    print *, "v2: ", v1

    !compute error from actual
    Err_rk4 = Y_act_rk4-Y_rk4(1, :)

    ! write to file 16, rk4.dat
    open(16, file="rk4.dat")
      do i = 1, np_rk4+1
        write(16,"(6(F18.12,'    '))")T_rk4(i),Y_act_rk4(i),&
                       DY_act_rk4(i),Y_rk4(1:2, i),Err_rk4(i)
      end do 
    close(16)
    print *, "V1: ", v1(3), "V2: ", v2(3)

    print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 2: Implicit RK3 Method "
    print *, " "

    F(1, :) = (/-1.0, 3.0,-5.0, 7.0 /)     
    F(2, :) = (/ 0.0,-2.0, 4.0,-6.0 /)     
    F(3, :) = (/ 0.0, 0.0,-4.0, 6.0 /)     
    F(4, :) = (/ 0.0, 0.0, 0.0,-16.0/)     

    np1 = int(t_fin/dt_crit-dt_pert)
    np2 = int(t_fin/(dt_crit+dt_pert))
    np3 = int(t_fin/dt)

    allocate(Y(ma,np1+1), T(np1+1), Y2(ma,np2+1), T2(np2+1))
    allocate(Y3(np3+1, ma), T3(np3+1))

    Y(:,1) = 1.0
    Y2(:,1) = 1.0
    !Y3(1, :) = 1.0

    T(1) = 0.0
    T2(1) = 0.0
    !T3(1) = 0.0

    call MidtermRK3(F, Y, T, ma, np1, dt_crit-dt_pert)
    call MidtermRK3(F, Y2, T2, ma, np2, dt_crit+dt_pert)
    !call AB3(F, Y3, T3, ma, np3, dt)
   
    open(15, file="RK3.dat")
    if (np2 .gt. np1) then 
      do i = 1, np1
        write(15, "(10(F24.8, '    '))") T2(i), Y2(:,i), T(i), Y(:,i)
      end do 
      do i = np1+1, np2
        write(15, "(5(F24.8, '    '))") T2(i), Y2(:,i)
      end do 
    else 
      do i = 1, np2
        write(15, "(10(F24.8, '    '))") T2(i), Y2(:,i), T(i), Y(:,i)
      end do 
    end if
    close(15)

    !open(17, file="prob2_act.dat")
      !do i = 1, np3+1
        !write(17, "(5(F16.8, '    '))") T3(i), Y3(i, :)
      !end do 
    !close(17)

    print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 3: LMM "
  print *, " "

  ! proving consistency
  A(1, :) = (/ 1.0, -1./3., -1./3./)
  A(2, :) = (/ 1.0,  1./3.,  2./3./)
  A(3, :) = (/ 1./2., -1./6., -2./3./)
  B(1,1) = 1./3.
  B(2,1) = 2.
  B(3,1) = -1./3.
  X = 0.0

  call LU(A, m3, bool, P3)

  call LUsolve(A, m3, B, X, 1, P3)

  call printmat(X, m3, 1)

  print *, "Order 2: ", X(1,1)/6. + X(2,1)/18. + X(3,1)*8./18. - 5./12.
  print *, "Order 3: ", X(1,1)/24. - X(2,1)/3. - X(3,1)*16./3. + 2./9.

end program midterm_driver
