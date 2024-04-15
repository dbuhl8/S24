

program driver

  use LinAl
  use NumDE

  implicit none

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
    print *, N(i)

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
  
  


end program driver



