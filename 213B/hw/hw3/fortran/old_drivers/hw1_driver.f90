

program driver

  use LinAl
  use NumDE

  implicit none

  ! problem 1
  integer, parameter :: n1=20, n2=60, n3=200, ne=200, ae = 1, be = 10000
  real, parameter :: a = 0.0, b = 1.0
  integer :: i, j, de
  real :: dx
  real, dimension(n1+4) :: X1, F1
  real, dimension(n2+4) :: X2, F2
  real, dimension(n3+1) :: Xa, Fa
  real, allocatable :: XE(:), FE(:), FaE(:), XaE(:)
  real, dimension(ne) :: E
  integer, dimension(ne) :: N
  ! problem 2
  integer, parameter :: num_points = 1000, ma = 2
  real, dimension(ma, num_points+3) :: Y
  real, dimension(ma, ma) :: F
  real, dimension(num_points+3) :: T
  real, dimension(3) :: rB, rC
  real, dimension(3, 3) :: rA
  real :: dt, tf = 10.
  real, dimension(4) :: dt_array
  real, allocatable :: YE(:, :), TE(:), ER(:), EA(:)
  integer :: np
  

  print *, " "
  print *, "--------------------------------------------------------------------------"
  print *, " "
  print *, " Question 1: Finite Difference Backwards Differentiation Method: BDF3"
  print *, " "

  ! Ploting the n=20, n=60 approx

  dx = (b-a)/n1
  do i = 0, n1+3
    X1(i+1) = dx*(i-3)
  end do
  dx = (b-a)/n2
  do i = 0, n2+3
    X2(i+1) = dx*(i-3)
  end do
  do i = 0, n3
    Xa(i+1) = i/real(n3)
  end do

  F1 = log10(2. + sin(2.*pi*X1))
  F2 = log10(2. + sin(2.*pi*X2))
  Fa = (2.*pi*cos(2*pi*Xa))/(log(10.)*(2. + sin(2.*pi*Xa)))



  dx = (b-a)/n1
  call BDF3(F1, n1, dx)
  
  dx = (b-a)/n2
  call BDF3(F2, n2, dx)

  open(11, file='df.dat')
    do i = 1, n1+1
        write(11, "(6F8.4)") Xa(i), Fa(i), X2(i+3), F2(i), X1(i+3), F1(i)
    end do
    do i = n1+2, n2+1
        write(11, "(4F8.4)") Xa(i), Fa(i), X2(i+3), F2(i)
    end do
    do i = n2+2, n3+1
        write(11, "(2F8.4)") Xa(i), Fa(i)
    end do
  close(11)

  ! Plotting the error
  de = (be-ae)/ne
  do i = 1, ne

    N(i) = ae + de*(i-1)

    allocate(XE(N(i)+4), FE(N(i)+4), FaE(N(i)+4), XaE(N(i)+4))

    dx = (b-a)/N(i)

    do j = 0, N(i)+3
        XE(j+1) = dx*(j-3)
    end do

    FE = log10(2. + sin(2.*pi*XE))
    XaE = 0.0
    do j = 0, N(i)
        XaE(i+1) = dx*i
    end do
    FaE = (2.*pi*cos(2*pi*XE))/(log(10.)*(2. + sin(2.*pi*XE)))

    call BDF3(FE, N(i), dx)

    E(i) = maxval(abs(FaE(4:N(1)+4)-FE(1:N(1)+1)))

    deallocate(XE, FE, FaE, XaE)

  end do
 
  open(12, file="e.dat") 

    do i = 1, ne
        write(12, "(I8, 2F20.15)") N(i), E(i), real(N(i))**(-3)
    end do
  close(12)


  print *, " "
  print *, "--------------------------------------------------------------------------"
  print *, " "
  print *, " Question 2: Implicit/Explicit Multistep Methods for ODEs"
  print *, " "

  T = 0.
  F(1, :) = (/ -1.0, 3.0 /)
  F(2, :) = (/ -3.0, -1.0 /)

  Y(1, 1) = -3.0
  Y(2, 1) = 1.0

  dt = 0.01

  ! call RK3 method here
  
  ! defining the RK3 Coefficient Matrices
  rA(1, :) = (/ 0.0, 0.0, 0.0/)
  rA(2, :) = (/ 1./3., 0.0, 0.0/)
  rA(3, :) = (/ 0.0, 2./3., 0.0/)
   
  rB = (/1./4., 0., 3./4./)

  rC = (/ 0.0, 1./3., 2./3./)

  ! Runge Kutta 3
  call RK3(F, Y, T, ma, num_points, rA, rB, rC, dt)

  open(12, file="rk3.dat")

    do i = 1, num_points
        write(12, "(3F8.4)") T(i), Y(:, i)
    end do

  close(12)

  open(14, file="actual.dat")
    
    do i = 1, num_points

        write(14, "(3F8.4)") T(i), (-3.*cos(3.*T(i)) + sin(3.*T(i)))*exp(-T(i)), (3.*sin(3.*T(i))+cos(3.*T(i)))*exp(-T(i))

    end do

  close(14)

  ! Adams Moulton 3
  Y = 0.
  T = 0.
  Y(1, 1) = -3.0
  Y(2, 1) = 1.0

  ! need to get first 2 points with RK3
  call RK3(F, Y, T, ma, 2, rA, rB, rC, dt)
  call AM3(F, Y, T, ma, num_points-2, dt)

  open(13, file="am3.dat")

    do i = 1, num_points
        write(13, "(3F8.4)") T(i), Y(:, i)
    end do

  close(13)

  dt_array = (/0.1, 0.05, 0.005, 0.0005/)

  do i = 1, 4

    np = tf/dt_array(i)

    allocate(YE(2, np+1), TE(np+1), ER(np+1), EA(np+1))

    YE = 0.
    TE = 0.
    YE(1, 1) = -3.0
    YE(2, 1) = 1.0

    call RK3(F, YE, TE, ma, np, rA, rB, rC, dt_array(i))
   
    YE(1, :) = YE(1, :) - (-3.*cos(3.*TE) + sin(3.*TE))*exp(-TE)
    YE(2, :) = YE(2, :) - (3.*sin(3.*TE)+cos(3.*TE))*exp(-TE)
 
    call vectwonorm(YE, ER, np)

    YE = 0.
    TE = 0.
    YE(1, 1) = -3.0
    YE(2, 1) = 1.0

    call RK3(F, YE, TE, ma, 2, rA, rB, rC, dt_array(i))
    call AM3(F, YE, TE, ma, np-2, dt_array(i))
 
    YE(1, :) = YE(1, :) - (-3.*cos(3.*TE) + sin(3.*TE))*exp(-TE)
    YE(2, :) = YE(2, :) - (3.*sin(3.*TE)+cos(3.*TE))*exp(-TE)
   
    call vectwonorm(YE, EA, np)

    open(15, file="error"//trim(str(i))//".dat")
    do j = 1, np+1
       write(15, "(3F20.15)") TE(j), ER(j), EA(j)
    end do
    close(15)

    if (i .eq. 1) then
        open(16, file="finalerror.dat")
            write(16, "(5F20.15)") dt_array(i), ER(np), real(np)**(-3), EA(np), real(np)**(-4)
        close(16)
    else
        open(16, file="finalerror.dat", status="old", position="append")
            write(16, "(5F20.15)") dt_array(i), ER(np), real(np)**(-3), EA(np), real(np)**(-4)
        close(16)
    end if
    deallocate(YE, TE, ER, EA)

  end do   

end program driver




