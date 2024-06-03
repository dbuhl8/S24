!File: hw4_driver.f90
!Author: Dante Buhl
!Depencies: LinAl, NumDe

program hw4_driver

  use LinAl 
  use NumDE

  implicit none
 
  ! Gen Util Vars
  integer, parameter :: kr=kind(dble(1.0))
  integer :: i, j, k, l
  real(kind=kr) :: start, finish
  integer :: fn = 10, ind

  ! Combined Vars
  real(kind=kr), parameter :: dt=0.0005
  integer, parameter :: nx=80, ny=80, na=(nx)*(ny), nt=int(2./dt), nt_p=9
  integer, parameter :: nn_p = 4
  integer, dimension(nn_p), parameter :: mse_n=(/40, 60, 100, 140/)

  ! Question 2 Vars
  real(kind=kr), allocatable :: A(:,:), TA(:)
  real(kind=kr), dimension(3,1) :: IC
  real(kind=kr), dimension(3,4) :: K4
  real(kind=kr), dimension(4) :: B, C
  real(kind=kr), dimension(ny,nx,nt_p) :: u_sol
  real(kind=kr), dimension(nx*ny,1) :: cx, cy, cu
  real(kind=kr), dimension(ny,nx) :: cmx, cmy
  real(kind=kr), dimension(nt_p) :: char_int
  integer :: temp_nt 

  ! Question 3 Vars
  real(kind=kr) :: dx, dy
  ! vectorized vars
  real(kind=kr), dimension(na, nt+1) :: U
  real(kind=kr), dimension(na, 1) :: F, G, X, Y
  real(kind=kr), dimension(nt+1) :: T
  real(kind=kr), dimension(na, na) :: D1_x, D1y, D
  real(kind=kr), dimension(nt_p) :: timeindx, int_val, t_val
  integer, dimension(1) :: indx

  ! non-vectorized vars
  real(kind=kr), dimension(nx, ny, nt+1) :: UU
  real(kind=kr), dimension(nx, ny) :: XX, YY

  real(kind=kr), dimension(9, 9) :: PDx, PDy

  ! Question 3, Part 4 Vars
  real(kind=kr), allocatable :: mse_u_fd(:,:), mse_u_ch(:,:), mse_t(:)
  real(kind=kr), allocatable :: mse_x(:,:), mse_y(:,:), mse_f(:,:), mse_g(:,:)
  real(kind=kr), allocatable :: mse_D1y(:,:), mse_D1_x(:,:), mse_D(:,:)
  real(kind=kr) :: mse_dx, mse_dy
  real(kind=kr), dimension(nn_p) :: mse_error
  real(kind=kr), dimension(nt_p) :: mpe_val
  integer :: mn, mse_nt, indt
  

  print *, " "
  print *, " Question 1 "
    print *, " "
    print *, " Woah there isn't any coding in this problem "

    call PD1xFD2_2d(PDx, 3, 3, 1.)
    call PD1yFD2_2d(PDy, 3, 3, 1.)
    print *, " PDx : " 
    call printmat(PDx, 9, 9)
    print *, " PDy : " 
    call printmat(PDy, 9, 9)
    print *, " "
  print *, " "
  print *, "------------------------------------------------------------------"
  print *, " "
  print *, " Question 2 "
    print *, " "
    
    timeindx = (/0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2./)
  
    B = (/1./6., 1./3., 1./3., 1./6. /)
    C = (/ 0.0, 1./2., 1./2., 1.0 /)
    K4 = 0.

    call vec_domain(cx, cy, nx, ny, (/0.,2*pi/),(/0.,2*pi/))

    cu = (1./(2.*pi**2))*(sin(cx+cy)**2)
    char_int(1) = sum(cu)*4.*(pi**2)/(nx*ny)
    call devectorize(cu, u_sol(:,:,1), nx, ny)
    call devectorize(cx, cmx, nx, ny)
    call devectorize(cy, cmy, nx, ny)

    do l = 2, nt_p
      temp_nt = int(timeindx(l)/dt)
      allocate(A(3, temp_nt+1), TA(temp_nt+1))
      print *, "Time Val : ", timeindx(l)
      do k = 1, ny
        do j = 1, nx
          ind = k + ny*(j-1)
          A(2,1) = cy(ind,1)
          A(1,1) = cx(ind,1)
          TA(1) = timeindx(l)
          do i = 1, temp_nt
            ! find k vec
            K4(1:2,1) = back_sys(A(1:2,i)) 
            K4(1:2,2) = back_sys(A(1:2,i) + dt*K4(1:2,1)/2.)
            K4(1:2,3) = back_sys(A(1:2,i) + dt*K4(1:2,2)/2.)
            K4(1:2,4) = back_sys(A(1:2,i) + dt*K4(1:2,3))
            !update Y  
            A(1:2,i+1) = A(1:2,i) + dt*(b(1)*K4(1:2,1) + b(2)*K4(1:2,2) +&
                                    b(3)*K4(1:2,3) + b(4)*K4(1:2,4))
            TA(i+1) = TA(i) - dt
          end do 
          ! check to see if it is doing the right thing, 
          A(1:2,1) = A(1:2, temp_nt+1)
          A(3,1) = (1./(2*pi**2))*(sin(A(1,1)+A(2,1))**2)
          TA = 0.0
          do i = 1, temp_nt
            ! find k vec
            K4(:,1) = forw_sys(A(:,i)) 
            K4(:,2) = forw_sys(A(:,i) + dt*K4(:,1)/2.)
            K4(:,3) = forw_sys(A(:,i) + dt*K4(:,2)/2.)
            K4(:,4) = forw_sys(A(:,i) + dt*K4(:,3))
            !update Y  
            A(:,i+1) = A(:,i) + dt*(b(1)*K4(:,1) + b(2)*K4(:,2) +&
                                    b(3)*K4(:,3) + b(4)*K4(:,4))
            TA(i+1) = TA(i) + dt
          end do
          ! write solution at t = t*
          cu(ind,1) = A(3,temp_nt+1)
        end do 
      end do 
      deallocate(A, TA)
      call devectorize(cu, u_sol(:,:,l), nx, ny)
      ! compute error and integral 
      char_int(l) = sum(cu)*4.*(pi**2)/(nx*ny)
      print *, " Integrates to ", char_int(l)
    end do 

    open(fn, file='char_int.dat')
      do i = 1, nt_p
        write(fn, *) timeindx(i), char_int(i)  
      end do
    close(fn)
    fn = fn+1
    open(fn, file='char_u.dat')
      do i = 1, nt_p
        call writemat(u_sol(:,:,i), ny, nx, fn)
      end do 
    close(fn)
    fn = fn+1
    open(fn, file='char_x.dat')
      call writemat(cmx, ny, nx, fn)
    close(fn)
    fn = fn+1
    open(fn, file='char_y.dat')
      call writemat(cmy, ny, nx, fn)
    close(fn)
    fn = fn+1


    
    print *, " "
  print *, " "
  print *, "------------------------------------------------------------------"
  print *, " "
  print *, " Question 3 "
    print *, " "
  
    ! Intializing Vectorized Domain
    call vec_domain(X,Y,nx,ny,(/2.*pi/80.,2*pi/),&
      (/2.*pi/80.,2.*pi/)) 

    ! Computing Differentials
    dx = 2.*pi/80.
    dy = dx

    ! Intializing forcing vectors
    F = sin(X)*sin(Y)
    G = 1. - exp(sin(X+Y)) 

    ! Initializing Initial Condition
    U(:, 1:1) = (1./(2.*(pi**2)))*(sin(X + Y)**2)
    T(1) = 0.0

    ! Initializing D1 matrices
    call PD1xFD2_2d(D1_x, nx, ny, dx)
    call PD1yFD2_2d(D1y, nx, ny, dy)

    ! Adjusting the matrices according to the problem
    do i = 1, na
      D1_x(:,i) = D1_x(:,i)*F(i,1)
      D1y(:,i) = D1y(:,i)*G(i,1)
    end do 

    ! Getting matrix for timestepping method
    D = -D1_x - D1y

    ! Calling AB2
    call AB2(D, U, T, na, nt, dt)

    ! Devectorize
    do i = 1, nt+1
      call devectorize(U(:,i:i), UU(:,:,i), nx, ny)
    end do 
    call devectorize(X, XX, nx, ny)
    call devectorize(Y, YY, nx, ny)

    
    do i = 1, nt_p
      indt = int(timeindx(i)/dt)+1
      mpe_val(i) = maxval(abs(UU(:,:,indt) - u_sol(:,:,i)))
      int_val(i) = 4*(pi**2)*sum(UU(:,:,indt))/(80**2)
    end do 

    ! Write output
    open(fn, file='fdint.dat')
      do i = 1, nt_p
        write(fn, *) timeindx(i), int_val(i)
      end do 
    close(fn)
    fn = fn+1
    open(fn, file='mpe_error.dat')
      do i = 1, nt_p
        write(fn, *) timeindx(i), mpe_val(i)
      end do 
    close(fn)
    fn = fn+1
    open(fn, file='x.dat')
      call writemat(XX, ny, nx, fn)
    close(fn)
    fn = fn+1
    open(fn, file='y.dat')
      call writemat(YY, ny, nx, fn)
    close(fn)
    fn = fn+1
    open(fn, file='u.dat')
      do i = 1, nt+1
        call writemat(UU(:,:,i), ny, nx, fn)
      end do 
    close(fn)
    fn = fn+1
    open(fn, file='t.dat')
      do i = 1, nt+1
        write(fn, '(F10.6)') T(i)
      end do
    close(fn)
    fn = fn+1
    open(fn, file='params.dat')
      write(fn, "(3I6)") nx, ny, nt+1
    close(fn)

    ! computing MSE for differnt grid sizes
    mse_nt = int(1./dt) 
    do l = 1, nn_p
      mn = mse_n(l)
      !Allocating and Initializing Domain for Solvers
      allocate(mse_u_fd(mn**2,mse_nt+1),mse_u_ch(mn**2,1),&
        mse_y(mn**2,1),mse_x(mn**2,1), mse_t(mse_nt+1))
      allocate(mse_f(mn**2,1), mse_g(mn**2,1))
      call vec_domain(mse_x,mse_y,mn,mn,(/2.*pi/mn,2.*pi/),(/2.*pi/mn,2.*pi/))

      ! Method of Characteristics Solver for t = 1.0
      allocate(A(3, mse_nt+1))
      do k = 1, mn
        do j = 1, mn
          ind = k + mn*(j-1)
          A(2,1) = mse_y(ind,1)
          A(1,1) = mse_x(ind,1)
          mse_t(1) = 1.0
          do i = 1, mse_nt
            ! find k vec
            K4(1:2,1) = back_sys(A(1:2,i)) 
            K4(1:2,2) = back_sys(A(1:2,i) + dt*K4(1:2,1)/2.)
            K4(1:2,3) = back_sys(A(1:2,i) + dt*K4(1:2,2)/2.)
            K4(1:2,4) = back_sys(A(1:2,i) + dt*K4(1:2,3))
            !update Y  
            A(1:2,i+1) = A(1:2,i) + dt*(b(1)*K4(1:2,1) + b(2)*K4(1:2,2) +&
                                    b(3)*K4(1:2,3) + b(4)*K4(1:2,4))
            mse_t(i+1) = mse_t(i) - dt
          end do 
          ! check to see if it is doing the right thing, 
          A(1:2,1) = A(1:2, mse_nt+1)
          A(3,1) = (1./(2*pi**2))*(sin(A(1,1)+A(2,1))**2)
          mse_t = 0.0
          do i = 1, mse_nt
            ! find k vec
            K4(:,1) = forw_sys(A(:,i)) 
            K4(:,2) = forw_sys(A(:,i) + dt*K4(:,1)/2.)
            K4(:,3) = forw_sys(A(:,i) + dt*K4(:,2)/2.)
            K4(:,4) = forw_sys(A(:,i) + dt*K4(:,3))
            !update Y  
            A(:,i+1) = A(:,i) + dt*(b(1)*K4(:,1) + b(2)*K4(:,2) +&
                                    b(3)*K4(:,3) + b(4)*K4(:,4))
            mse_t(i+1) = mse_t(i) + dt
          end do
          ! write solution at t = t*
          mse_u_ch(ind,1) = A(3,mse_nt+1)
        end do 
      end do 
      deallocate(A)
    
      ! Finite Difference Solver
      allocate(mse_D1_x(mn**2, mn**2),mse_D1y(mn**2, mn**2),mse_D(mn**2,mn**2))
      mse_dx = 2*pi/mn
      mse_dy = 2*pi/mn

      ! Intializing forcing vectors
      mse_F = sin(mse_X)*sin(mse_Y)
      mse_G = 1. - exp(sin(mse_X+mse_Y)) 

      ! Initializing Initial Condition
      mse_u_fd(:, 1:1) = (1./(2.*(pi**2)))*(sin(mse_X + mse_Y)**2)
      mse_t(1) = 0.0

      ! Initializing D1 matrices
      call PD1xFD2_2d(mse_D1_x, mn, mn, mse_dx)
      call PD1yFD2_2d(mse_D1y, mn, mn, mse_dy)

      ! Adjusting the matrices according to the problem
      do i = 1, mn**2
        mse_D1_x(:,i) = mse_D1_x(:,i)*mse_F(i,1)
        mse_D1y(:,i) = mse_D1y(:,i)*mse_G(i,1)
      end do 

      ! Getting matrix for timestepping method
      mse_D = -mse_D1_x - mse_D1y

      ! Calling AB2
      call AB2(mse_D, mse_u_fd, mse_t, mn**2, mse_nt, dt)

      mse_error(l) = 4*(pi**2)*sum((mse_u_fd(:,mse_nt+1)-mse_u_ch(:,1))**2)/&
        (mn**2)
      deallocate(mse_t,mse_x, mse_y, mse_f,mse_g,mse_D,mse_D1y,mse_D1_x,&
        mse_u_fd,mse_u_ch)
    end do 

    open(fn,file='mse_error.dat')  
      do i = 1, nn_p
        write(fn, *) mse_n(i), mse_error(i)
      end do 
    close(fn)

    print *, " "
  print *, " "
  print *, "------------------------------------------------------------------"

  contains

    function back_sys(y_ode) result(dydt)

      real(kind=kr) :: y_ode(:)
      real(kind=kr), dimension(2) :: dydt    

      dydt(1) = -sin(y_ode(1))*sin(y_ode(2))
      dydt(2) = -1. + exp(sin(sum(y_ode)))

    end function back_sys

    function forw_sys(u_ode) result(dydt)

      real(kind=kr) :: u_ode(:)
      real(kind=kr), dimension(3) :: dydt    

      dydt(1) = sin(u_ode(1))*sin(u_ode(2))
      dydt(2) = 1. - exp(sin(u_ode(1)+u_ode(2)))
      dydt(3) = -(cos(u_ode(1))*sin(u_ode(2)) - cos(u_ode(1)+u_ode(2))*&
        exp(sin(u_ode(1)+u_ode(2))))*u_ode(3)

    end function forw_sys
 
end program hw4_driver

