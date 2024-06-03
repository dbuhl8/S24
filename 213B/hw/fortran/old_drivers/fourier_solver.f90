!File: hw4_driver.f90
!Author: Dante Buhl
!Depencies: LinAl, NumDe

program hw4_driver

  use LinAl 
  use NumDE

  implicit none
 
  ! Gen Util Vars
  integer, parameter :: kr=kind(dble(1.0))
  integer :: i, j
  real(kind=kr) :: start, finish

  ! FFT Testing Vars 
  integer, parameter :: nx=200, nt = 20000
  real(kind=kr) :: xstart, xstop, dx, dt
  real(kind=kr), dimension(nx) :: k
  real(kind=kr), dimension(nx,nt+1) :: u, x, t
  real(kind=kr), dimension(nt+1) :: ts
  complex(kind=kr), dimension(nx,nt+1) :: u_spec
  complex(kind=kr), dimension(nx,nx) :: D2

  ! Question 1 Vars

  ! Question 2 Vars


  print *, "------------------------------------------------------------------"
  print *, " "
  print *, " Fourier Spectral Solver Outputs"
    print *, " "
    ! Initializing Domain Vars
    xstart = -25
    xstop = 25.
    dx  = (xstop-xstart)/(nx-1)
    do i = 1, nx
      x(i,1) = xstart + (i-1)*dx
    end do
    dt = 10.**(-4)
    ts(1) = 0.

    ! Initializing IC
    u(:,1) = sin(X(:,1))*exp(-((x(:,1)-10)**2)/2.0)
    !5.*(1. - x(:,1)**2)**2
    call ft1d(u(:,1), u_spec(:,1), x(:,1), k, nx, dx)
    D2 = 0.0
    call D2fourier(D2, k, nx)

    ! Timesteping
    call CRK4(D2, u_spec, ts, nx, nt, dt)

    ! calling and inverse fourier transform
    do i = 0, nt
      call ift1d(u(:,i+1), u_spec(:,i+1), x(:,1), k, nx)
      u(:,i+1) = u(:,i+1) !+ 3. + X(:,1)
    end do
    ! writing out domain matrices 
    do i = 1, nt
      x(:,i+1) = x(:,1)
    end do 
    do i = 1, nx
      t(i,:) = ts
    end do 

    ! Writing to file for plotting
    open(12, file='u_spec.dat')
      do i = 1, nx
          write(12, "("//trim(str(2*(nt+1)))//"F20.10)") u_spec(i,:)
      end do 
    close(12)
    open(15, file='u_phys.dat')
      do i = 1, nx
          write(15, "("//trim(str(nt+1))//"F20.10)") u(i,:)
      end do 
    close(15)
    open(13, file='x.dat')
      do i = 1, nx
        write(13, "("//trim(str(nt+1))//"F12.6)") x(i,:)
      end do 
    close(13)
    open(14, file='t.dat')
      do i = 1, nx
        write(14, "("//trim(str(nt+1))//"F12.6)") t(i,:)
      end do 
    close(14)

    print *, " "
  print *, " "
  print *, "------------------------------------------------------------------"
  print *, " "
  print *, " Question 1 "
    print *, " "
    print *, " "
  print *, " "
  print *, "------------------------------------------------------------------"
  print *, " "
  print *, " Question 2 "
    print *, " "
    print *, " "
  print *, " "
  print *, "------------------------------------------------------------------"
  
end program hw4_driver



