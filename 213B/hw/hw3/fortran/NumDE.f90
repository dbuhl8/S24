!File: Numerical Differential Equations Module
!Author: Dante Buhl
!Dependencies: IBVP_1DCN uses a LAPACK routine for a linear solve
!Purpose: To store numerical methods for Differential Equations

module NumDE
    
  use LinAl

  implicit none

  contains 

  subroutine BDF3(F, num_points, dx)
    ! assumes an evenly spaced grid
    
    ! At Input:
    ! F should be a vector containing at num_points + 3 entries such that each
    ! entry is F(i) = f(i), where f(x) is any function we want to differentiate.
    ! The first 3 points should be padding points, i.e. points used to calculate
    !  the first derivative points at F(4)
    ! F(1:num_points+3) = f(X(1:num_points+3)), where X(4) = x_0
    
    ! dx should be the step size between adjecent entries of.

    ! on Output: 
    ! F(1:num_points) are the values of the derivative df/dx
    ! according to the BDF3 Method

    implicit none

    real :: F(:), dx
    integer :: i, num_points

    do i = 1, num_points+1
        F(i) = (11.*F(i+3) - 18.*F(i+2) + 9.*F(i+1) - 2.*F(i))/(6.*dx)
    end do
  end subroutine BDF3 

  subroutine RK3(F, Y, T, ma, num_points, A, B, C, dt)
    implicit none
                        
    real :: Y(:, :), F(:, :), A(:, :), B(:), C(:), dt, T(:)
    integer, intent(in) :: ma, num_points
    integer :: i, j
    real, dimension(3, ma) :: K

    do i = 1, num_points
      K = 0.0
  
      ! update k vectors
      do j = 1, 3
          K(j:j ,:) = matmul(F, Y(i:i, :) + dt*matmul(K, transpose(A(j:j,:))))
      end do

      ! update y
      Y(i+1, :) = Y(i, :)
      do j = 1, 3
          Y(i+1, :) = Y(i+1, :) + dt*b(j)*K(j, :)
      end do

      T(i+1) = T(i) + dt
    end do
  end subroutine RK3

  subroutine MidtermRK3(F, Y, T, ma, num_points, dt)
    implicit none
 
    real :: Y(:, :), F(:, :), dt, T(:)
    real, dimension(3) :: B
    integer, intent(in) :: ma, num_points
    integer :: i, j
    real, dimension(ma, 3) :: K

    !vars for LU solve
    real, dimension(ma, ma) :: LU_A, eye
    real, dimension(ma, 1) :: LU_B, col
    logical :: bool
    integer, dimension(ma) :: P

    B = (/ 1./6, 2./3, 1./6 /)
 
    call ident(eye, ma)
    do i = 1, num_points
      LU_A = eye - dt*F/4.0
      call LU(LU_A, ma, bool, P)
      K = 0.0

      col=0.0
      ! update k vectors
      K(:,1) = matmul(F, Y(:,i))
      LU_B = matmul(F, Y(:,i:i) + dt*0.25*K(:, 1:1))
      call LUsolve(LU_A, ma, LU_B, K(:,2:2), 1, P)
      K(:,3) = matmul(F, Y(:,i) + dt*K(:,2))

      ! update y
      Y(:,i+1) = Y(:,i)+dt*(b(1)*k(:,1)+b(2)*k(:,2)+b(3)*k(:,3))

      T(i+1) = T(i) + dt
    end do
  end subroutine MidtermRK3

  subroutine MidtermRK4(F, Y, T, ma, num_points, dt)
    implicit none
  
    real :: F(:, :), Y(:, :), T(:), dt
    real, dimension(ma, 4) :: K
    real, dimension(4) :: B, C
    integer :: ma, num_points, i
    
    B = (/1./6., 1./3., 1./3., 1./6. /)
    C = (/ 0.0, 1./2., 1./2., 1.0 /)
    K = 0.

    do i = 1, num_points
      ! find k vec
      !if(i.eq.2) then
        !print *, " "
        !call printmat(K, ma, 4)
        !print *, " "
      !end if
      K(:,1) = matmul(F, Y(:,i))
      K(ma,1) = K(ma,1) + (T(i) + C(1)*dt)**2
      K(:,2) = matmul(F, Y(:, i) + dt*K(:,1)/2.)
      K(ma,2) = K(ma,2) + (T(i) + C(2)*dt)**2
      K(:,3) = matmul(F, Y(:, i) + dt*K(:,2)/2.)
      K(ma,3) = K(ma,3) + (T(i) + C(3)*dt)**2
      K(:,4) = matmul(F, Y(:, i) + dt*K(:,3))
      !print *, K(ma, 4)
      K(ma,4) = K(ma,4) + (T(i) + C(4)*dt)**2
      !print *, K(ma, 4)
      !if(i.eq.1) then
        !print *, " "
        !call printmat(K, ma, 4)
        !print *, " "
      !end if
      !update Y  
      Y(:,i+1) = Y(:,i) + dt*(b(1)*k(:,1) + b(2)*k(:,2) +&
                              b(3)*k(:,3) + b(4)*k(:,4))
      T(i+1) = T(i) + dt
    end do 
  end subroutine MidtermRK4

  subroutine AM3(F, Y, T, ma, num_points, dt)
    ! this subroutine can only handle linear ODE's as of right now
    ! That ism the jacobian matrix, A, is constant

    implicit none

    real :: Y(:, :), F(:, :), dt, T(:)
    integer, intent(in) :: ma, num_points
    real, dimension(ma, ma) :: Jac, eye
    real, dimension(ma, 1) :: B, X
    integer, dimension(ma) :: P
    integer :: i, j
    logical :: bool


    ! A = I - J_G(u_k+1)
    call ident(eye, ma)

    do i = 3, num_points+2
      ! implicit step, Fixed Point Iteration for u_k+1
      Jac = eye - (9.*dt/24.)*F
      call LU(Jac, ma, bool, P)
  
      B = Y(:, i:i) + (dt/24.)*(19.*matmul(F, Y(:, i:i)) - &
          5.*matmul(F, Y(:, i-1:i-1)) + matmul(F, Y(:, i-2:i-2)))
      X = 0.0
      call LUsolve(Jac, ma, B, X, 1, P)
      Y(:, i+1) = X(:, 1)

      ! AM3 Method
      Y(:, i+1) = Y(:, i) + (dt/24.0)*(9.*matmul(F, Y(:, i+1)) &
                  + 19.*matmul(F, Y(:, i)) - 5.*matmul(F, Y(:, i-1))  &
                  + matmul(F, Y(:, i-2)))
      T(i+1) = T(i) + dt
    end do
  end subroutine AM3

  subroutine AB1(F, Y, T, ma, num_points, dt)
    implicit none
    real :: F(:,:), Y(:,:),T(:),dt
    integer :: ma, num_points, i

    do i = 1, num_points
      Y(i+1,:) = Y(i+1,:)+dt*matmul(F,Y(i,:))
      T(i+1) = T(i) + dt
    end do
  end subroutine AB1

  subroutine AB2(F, Y, T, ma, num_points, dt)
    implicit none
    real :: F(:,:), Y(:,:),T(:),dt
    integer :: ma, num_points, i
    !calling Huen Method in order to get 2nd point
    call Huen(F, Y, T, ma, 1, dt)

    do i = 1, num_points-1
      Y(i+2,:) = Y(i+1,:)+(dt/2.)*(3*matmul(F,Y(i+1,:))-matmul(F,Y(i,:)))
      T(i+2) = T(i+1) + dt
    end do
  end subroutine AB2

  subroutine AB3(F, Y, T, ma, num_points, dt)
    implicit none
    real :: F(:,:), Y(:,:),T(:),dt
    integer :: ma, num_points, i
   
    !calling AB2 to get first 3 points
    call AB2(F, Y, T, ma, 2, dt)
     
    !starting the AB3 Scheme
    do i = 1, num_points-2
      Y(i+3,:) = Y(i+2,:)+(dt/12.)*(23*matmul(F,Y(i+2,:))-16*matmul(F,Y(i+1,:))&
                                      +5*matmul(F,Y(i,:)))
      T(i+3) = T(i+2) + dt
    end do 
  end subroutine AB3

  subroutine Huen(F, Y, T, ma, num_points, dt)
    implicit none
    real :: F(:,:), Y(:,:),T(:),dt
    integer :: ma, num_points, i

    do i = 1, num_points
      Y(i+1,:) = Y(i,:) + dt*matmul(F, Y(i,:))
      Y(i+1,:) = Y(i,:) + (dt/2.)*(matmul(F, Y(i,:)) + matmul(F, Y(i+1,:)))
      T(i+1) = T(i) + dt
    end do 
  end subroutine 
  
  subroutine NHuen(N, D, F, U, T, ma, num_points, dt)
    implicit none
    real :: N(:,:),D(:,:),F(:,:),U(:,:),T(:),dt
    integer :: ma, num_points, i

    do i = 1, num_points
      U(:,i+1) = U(:,i) + dt*&
        (U(:,i)*matmul(N, U(:,i)) + matmul(D, U(:,i)) + F(:,1))
      U(:,i+1) = U(:,i) + (dt/2.)*&
        ((U(:,i)*matmul(N, U(:,i)) + matmul(D, U(:,i)) + F(:,1)) +&
         (U(:,i+1)*matmul(N, U(:,i+1)) + matmul(D, U(:,i+1)) + F(:,1)))
      T(i+1) = T(i) + dt
    end do 
  end subroutine 

  subroutine FPI(F, Y)
    ! Fixed Point Interation, for a nonlinear system
    ! need to update the jacobian at each iteration (hard in fortran)
    real :: F(:, :), Y(:, :)
  end subroutine FPI

  subroutine D2FD2_2D(D, nx, ny, dx, dy)
    ! Computes D2 for an evenly spaced grid according to a second order finite
    ! difference method

    implicit none
 
    ! D should be an (nx*ny) x (nx*ny) array
    real :: D(:, :), dx, dy 
    integer :: nx, ny, i, n, m
    integer :: indx, indy
  
    D = 0.0

    indx = ny
    indy = 1
    do m = 1, ny
      do n = 1, nx
        i = m + ny*(n-1) 
        ! we require a lot of logic gates in order to make sure we are not
        ! violating the array
        D(i,i) = -2.0*(1/(dx**2) + 1/(dy**2))
        if(1 < n) then
          D(i, i-indx) = 1.0/(dx**2)
        end if
        if(n < nx) then
          D(i, i+indx) = 1.0/(dx**2)
        end if
        if(1 < m) then
          D(i, i-indy) = 1.0/(dy**2)
        end if
        if(m < ny) then
          D(i, i+indy) = 1.0/(dy**2)
        end if
      end do 
    end do 
  end subroutine

  subroutine D4FD2_1D(D, nx, dx)
    ! Computes D4 for an evenly spaced grid according to a second order finite
    ! difference method
    implicit none
 
    ! D should be an (nx) x (nx)array
    real :: D(:, :), dx
    integer :: nx, i, n
  
    D = 0.0
    do i = 1, nx
      ! we require a lot of logic gates in order to make sure we are not
      ! violating the array
      D(i,i) = 6.0
      if(1 < i) then
        D(i, i-1) = -4.0
      end if
      if(2 < i) then
        D(i, i-2) = 1.0
      end if
      if(i < nx) then
        D(i, i+1) = -4.0
      end if
      if(i < nx-1) then
        D(i, i+2) = 1.0
      end if
    end do 
    D = D/(dx**4)
  end subroutine

  subroutine D3FD2_1D(D, nx, dx)
    ! Computes D3 for an evenly spaced grid according to a second order finite
    ! difference method

    implicit none
 
    ! D should be an (nx*ny) x (nx*ny) array
    real :: D(:, :), dx
    integer :: nx, n
  
    D = 0.0

    do n = 1, nx
      if(1 < n) then
        D(n, n-1) = 1.0/(dx**3)
      end if
      if(2 < n) then
        D(n, n-2) = -1./(2.*dx**3)
      end if
      if(n < nx) then
        D(n, n+1) = -1.0/(dx**3)
      end if
      if(n < nx-1) then
        D(n, n+2) = 1./(2.*dx**3)
      end if
    end do 
  end subroutine

  subroutine D2FD2_1D(D, nx, dx)
    ! Computes D2 for an evenly spaced grid according to a second order finite
    ! difference method

    implicit none
 
    ! D should be an (nx*ny) x (nx*ny) array
    real :: D(:, :), dx
    integer :: nx, i, n
    integer :: indx
  
    D = 0.0

    indx = 1
    do i = 1, nx
      ! we require a lot of logic gates in order to make sure we are not
      ! violating the array
      D(i,i) = -2.0/(dx**2)
      if(1 < i) then
        D(i, i-indx) = 1.0/(dx**2)
      end if
      if(i < nx) then
        D(i, i+indx) = 1.0/(dx**2)
      end if
    end do 
  end subroutine

  subroutine D1FD2_1D(D, nx, dx)
    ! Computes D1 for an evenly spaced grid according to a second order finite
    ! difference method

    implicit none
 
    ! D should be an (nx*ny) x (nx*ny) array
    real :: D(:, :), dx
    integer :: nx, i
    D = 0.0
    do i = 1, nx
      ! we require a lot of logic gates in order to make sure we are not
      ! violating the array
      if(1 < i) then
        D(i, i-1) = 1./2.
      end if
      if(i < nx) then
        D(i, i+1) = -1./2.
      end if
    end do 
    D = D/dx
  end subroutine

  subroutine vec_boundary(bound,interior,corner,nc,nx, ny)
    ! Returns a populaated bound array with the indices of the boudary for a
    ! vectorized domain u
    ! bound should be an array of length (2*nx + 2*ny - 4)
    implicit none
    integer :: bound(:), nx, ny, i, j, k,n, m, interior(:)
    integer :: corner(:), l, nc(:), h
    j = 1
    k = 1
    l = 1
    h = 1
    do n = 1, nx
      do m = 1, ny
        i = m + ny*(n-1) 
        if(((n.eq.1).or.(n.eq.nx)).or.((m.eq.1).or.(m.eq.ny))) then
          bound(j) = i
          j = j + 1
        else 
          interior(k) = i
          k = k + 1
        end if
        if( ( (i.eq.1) .or. (i.eq.ny) ) .or. &
              ( (i.eq.nx*ny) .or. (i .eq. (ny*(nx-1)+1)) ) ) then
            corner(l) = i
            l = l + 1
          else 
            nc(h) = i
            h = h + 1
          end if
      end do 
    end do 
  end subroutine

  subroutine vec_domain(X, Y, nx, ny, xrange, yrange)
    !Returns a vectorized domain X, Y, similar to meshgrid but vectorized
    ! Assumes an evenly spaced grid
    ! X, Y should be arrays of length (nx*ny). 
    implicit none

    real :: X(:,:), Y(:,:), xrange(:), yrange(:)
    real :: dx, dy
    integer :: nx, ny, i, n, m

    dx = (xrange(2) - xrange(1))/(nx-1)
    dy = (yrange(2) - yrange(1))/(ny-1)

    do m = 1, ny
      do n = 1, nx
        i = m + ny*(n-1) 
        X(i,1) = xrange(1) + (n-1)*dx
        Y(i,1) = yrange(2) - (m-1)*dy
      end do 
    end do 
  end subroutine

  ! this subroutine needs to be tweaked
  subroutine meshgrid(XX, YY, X, Y, nx, ny)
    ! Returns a meshgrid of points in XX, YY 
    ! XX, YY should be (ny) x (nx) in size
    ! X should be an array of length nx
    ! Y should be an array of length ny
    implicit none

    real :: X(:), Y(:), XX(:,:), YY(:,:)
    integer :: nx, ny, i

    do i = 1, ny
      XX(i, :) = X
    end do 

    do i = 1, nx
      YY(:,i) = Y
    end do 
  end subroutine meshgrid

  subroutine devectorize(A, matA, nx, ny)
    ! Returns a devectorized version of A in matA
    ! A is an array of length nx*ny
    ! matA is a matrix of size (ny) x (nx)

    implicit none

    real :: A(:,:), matA(:,:)
    integer :: nx, ny, n, m

    do n = 1, nx
      do m = 1, ny
        matA(m, n) = A(m + ny*(n-1),1)
      end do 
    end do 
  end subroutine devectorize

  subroutine IBVP_1DCN(U, D, T, nx, nt, dt)
    !Solves the Initial Boundary Value Problem for some derivative operator D
    !and Forcing F (time-independant) using CN
    ! U should be an (nx) x (nt) 
    ! D should be (nx) x (nx)
    ! F should be (nx) x (1)
    implicit none

    real :: U(:, :), D(:, :), dt, T(:)
    integer :: nx, nt, i, j, n
    real, dimension(nx, nx) :: A, C, A_temp
    real, dimension(nx,1) :: B
    integer, dimension(nx) :: P
    logical :: bool
    integer :: info

    ! du/dt = perturbation from steady state
    ! need to update the interior at each timestep. So we have, 
    ! du/dt = D2u + f 
    
    ! CN scheme
    ! initializing global vars
    ! call some System Solver from LAPACK
    call ident(A, nx)
    C = A + (dt/2.)*D
    A = A - (dt/2.)*D
    call LU(A, nx, bool, P)
    do i = 1, nt
      !A_temp = A
      !U(:,i+1:i+1) = matmul(C, U(:, i:i))
      B = matmul(C, U(:, i:i))
      call LUsolve(A, nx, B, U(:,i+1:i+1), 1, P)
      !call dgesv(nx, 1, A_temp, nx, P, U(:,i+1:i+1), nx, info)
      T(i+1) = T(i) + dt
    end do 
  end subroutine IBVP_1DCN

  subroutine N_IC_1D_AB3(U, D, N, F, T, nx, nt, dt)
    ! Performs AB2 for a generally Nonlinear PDE. Has only order one nonlinear
    ! terms. That is the problem can be represented in the form. 
    ! dU/dt = U*(NU) + DU + F
    ! Where NU is multiplied component wise by U
    ! DU and F are linear terms.  
    implicit none

    real :: U(:,:), D(:,:), N(:,:), F(:,:), T(:), dt
    integer :: nx, nt,i

    call NHuen(N, D, F, U, T, nx, nt, dt)

    do i = 1, nt-1
      U(:,i+2) = U(:,i+1)+(dt/2.)*&
        (3*(U(:,i+1)*matmul(N, U(:,i+1)) + matmul(D,U(:,i+1)) + F(:,1))-&
        (U(:,i)*matmul(N, U(:,i)) + matmul(D,U(:,i)) + F(:,1)))
      T(i+2) = T(i+1) + dt
    end do
  end subroutine N_IC_1D_AB3

  subroutine find_empty_row(A, empty, notempty, ma, num_empty, num_notempty)
    ! Returns an array of row indices specifiying which rows of A have all zero
    ! entries. On entry empty and nonempty should be unallocated allocatable
    ! arrays
    implicit none 
    real :: A(:, :), norm
    integer :: ma, i, num_empty, num_notempty
    integer, allocatable :: empty(:), notempty(:)

    num_notempty = 0
    do i = 1, ma
      call twonorm(A(i, :), norm)
      if (norm .gt. 10.d-15) then
        num_notempty = num_notempty + 1  
      end if
    end do 
    if (num_notempty .ne. ma) then
      num_empty = ma - num_notempty
      allocate(empty(num_empty), notempty(num_notempty))

      num_empty = 0
      num_notempty = 0
      
      do i = 1, ma
        call twonorm(A(i, :), norm)
        if (norm .gt. 10.d-15) then
          num_notempty = num_notempty + 1
          notempty(num_notempty) = i
        else
          num_empty = num_empty + 1
          empty(num_empty) = i
        end if
      end do 
    end if
  end subroutine

  subroutine cosf(U, X, A, l, num_points, num_modes)
    ! Returns a fourier cosine series 
    ! On entry U should be empty and of the same size as X. 
    ! X should be a matrix filled with domain values. 
    ! Notice that A should be an array of size num_modes+1
    ! U = sum(A * cos((n*pi*X)/l))
    implicit none
    
    real :: U(:,:), X(:,:), A(:), l
    integer :: num_modes, num_points, i

    U = A(1)
    do i = 1, num_modes
      U = U + A(i+1)*cos((i*pi/l)*X)
    end do 
  end subroutine

  subroutine sinf(U, X, A, l, num_points, num_modes)
    ! Returns a fourier sin series 
    ! On entry U should be empty and of the same size as X. 
    ! X should be a matrix filled with domain values. 
    ! Notice that A should be an array of size num_modes
    ! U = sum(A * sin((n*pi*X)/l))
    implicit none
    
    real :: U(:,:), X(:,:), A(:), l
    integer :: num_modes, num_points, i

    do i = 1, num_modes
      U = U + A(i)*sin((i*pi/l)*X)
    end do 
  end subroutine

  subroutine GCLgrid_1D(X, nx)
    ! returns a grid of Gauss-Chebyshev-Lobatto points 
    ! Note that X should be a column vector of size (nx+1) x 1
    ! Also note that this grid is not vectorized in the same way that the other
    ! vectorized arrays are in vec_domain. I have to look into how this grad
    ! thing works. 

    ! 2D is going to be a nightmare. 
    
    implicit none

    real :: X(:, :)
    integer :: nx, i

    do i = 1, nx+1
      X(i,1) = -cos((i-1)*pi/nx)
    end do 
  end subroutine GCLgrid_1D

  subroutine GCLpoints(X, nx)
    ! returns a grid of Gauss-Chebyshev-Lobatto points 
    ! Note that X should be a column vector of size (nx+1) x 1
    ! Also note that this grid is not vectorized in the same way that the other
    ! vectorized arrays are in vec_domain. I have to look into how this grad
    ! thing works. 

    ! 2D is going to be a nightmare. 
    
    implicit none

    real :: X(:)
    integer :: nx, i

    do i = 1, nx+1
      X(i) = -cos((i-1)*pi/nx)
    end do 
  end subroutine GCLpoints

  subroutine GCLmeshgrid(XX, TT, nx, nt, dt, tstart)
    ! Returns a matrix XX, TT which are meshgrids of the domain
    ! XX and TT must be size (nx+1) x (nt+1) arrays

    implicit none

    real :: XX(:,:), TT(:,:), dt, tstart
    integer :: nx, nt, i, j

    do i = 1, nx+1
      XX(i,:) = -cos((i-1)*pi/nx)
    end do 
    do i = 0, nt
      TT(:,i+1) = tstart + i*dt
    end do 
  end subroutine GCLmeshgrid

  subroutine GCLdiffvec(d, nx)
    ! returns d determined by GCL spectral method
    ! d(1) = d(nx+1) = 2, d(:) = 1
    implicit none
    real :: d(:)
    integer :: nx
    d = 1.0
    d(1) = 2.0
    d(nx+1) = 2.0
  end subroutine GCLdiffvec

  subroutine D2GCL_1D(D2, nx)
    ! Returns the 'second derivative' matrix for the Gauss-Chebyshev-Lobatto
    ! coallocation method

    ! To all who wander here, this method is cursed. Look how awful this
    ! differentiation matrix is. Its so convoluted that you would have to derive
    ! it in order to understand it. 

    implicit none

    real :: D2(:,:)
    integer :: nx, i, j
    real, dimension(nx+1) :: d, x

    call GCLdiffvec(d, nx)
    call GCLpoints(x, nx)
    x = -x

    D2(1,1) = (real(nx)**4 -1.)/15.
    D2(nx+1,nx+1) = D2(1,1)

    do j = 2, nx+1
      D2(1,j)=(2.*((-1.)**(j-1))/(3.*d(j)))*&
        ((2*real(nx)**2+1)*(1.-x(j))-6)/&
        ((1.-x(j))**2)
    end do 
    do j = 1, nx
      D2(nx+1,j)=(2*((-1.)**(nx+j-1))/(3.*d(j)))*&
        ((2*real(nx)**2+1)*(1+x(j))-6)/&
        ((1+x(j))**2)
    end do 

    do j = 1, nx+1
      do i = 2, nx
        if(i.eq.j) then
          D2(i, j) = -((real(nx)**2 -1)*(1-x(i)**2) + 3)/&
            (3.*(1-x(i)**2)**2)
        else 
          D2(i, j) = ((-1.)**(i+j-2))*(x(i)**2+x(i)*x(j)-2)/&
            (d(j)*(1-x(i)**2)*(x(i) - x(j))**2)
        end if
      end do 
    end do 
  end subroutine D2GCL_1D

  subroutine GCL_2_phys(U, nx)
    ! converts U from GCL coefficients to Physical Values
    ! on entry U should be an array of length (nx+1)
    implicit none
    
    real :: U(:)
    integer :: nx, i
    real, dimension(nx+1) :: X, C
  
    C = U 
    U = 0.
    call GCLpoints(X, nx) 
    do i = 1, nx+1
      U = U + C(i)*GCLpolynomial(X, i, nx)
    end do 
  end subroutine GCL_2_phys 

  function GCLpolynomial(X, n, nx) result (F)
    ! returns the value of the m_th Gauss-Chebyshev-Lobatto Interopolating
    ! Polynomial
    implicit none
    
    integer :: n, nx, i
    real, dimension(nx+1) :: x, f, d

    call GCLdiffvec(d, nx)

    F = (-1.**(N+n+1))*sqrt(1. - X**2)*sin(nx*acos(X))/&
    (d(n)*(nx**2)*(X-X(n)))
  end function GCLpolynomial

  subroutine GCLvecdomain(X, Y, nx, ny, yrange)
    !Returns a vectorized domain X, T for the Gauss-Chebyshev-Lobatto grid
    !points
    ! X, Y should be arrays of length ((nx+1)*ny). 
    implicit none

    real :: X(:,:), Y(:,:), yrange(:)
    real :: dx, dy
    integer :: nx, ny, i, n, m

    dy = (yrange(2) - yrange(1))/(ny-1)

    do m = 1, ny
      do n = 1, nx+1
        i = m + ny*(n-1) 
        X(i,1) = -cos((n-1)*pi/nx) 
        Y(i,1) = yrange(2) - (m-1)*dy
      end do 
    end do 
  end subroutine GCLvecdomain

end module NumDE




