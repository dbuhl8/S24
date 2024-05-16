!File: Numerical Differential Equations Module
!Author: Dante Buhl
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

  subroutine FPI(F, Y)
    ! Fixed Point Interation, for a nonlinear system
    ! need to update the jacobian at each iteration (hard in fortran)
    real :: F(:, :), Y(:, :)
  end subroutine FPI

  subroutine D2FD2(D, nx, ny, dx, dy)
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

  subroutine vec_boundary(bound,interior, nx, ny)
    ! Returns a populaated bound array with the indices of the boudary for a
    ! vectorized domain u
    ! bound should be an array of length (2*nx + 2*ny - 4)
    implicit none
    integer :: bound(:), nx, ny, i, j, k,n, m, interior(:)
    j = 1
    k = 1
    do m = 1, ny
      do n = 1, nx
        i = m + ny*(n-1) 
        if(((n.eq.1).or.(n.eq.nx)).or.((m.eq.1).or.(m.eq.ny))) then
          bound(j) = i
          j = j + 1
        else 
          interior(k) = i
          k = k + 1
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
  subroutine vec_meshgrid(X, Y, nx, ny, sx, sy)
    !Returns a vectorized domain X, Y, according to a meshgrid sx, sy
    ! Assumes an evenly spaced grid
    ! X, Y should be arrays of length (nx*ny) and ordered such that 
    ! X(1) < x < X(nx), Y(1) < y < Y(ny)
    implicit none

    real :: X(:), Y(:), sx(:), sy(:)
    integer :: nx, ny, i, n, m

    do m = 1, ny
      do n = 1, nx
        i = m + ny*(n-1) 
        X(i) = sx(n)
        Y(i) = sy(ny-m+1)
      end do 
    end do 
  end subroutine

  subroutine devectorize(A, matA, nx, ny)
    ! Returns a devectorized version of A in matA
    ! A is an array of length nx*ny
    ! matA is a matrix of size (nx) x (ny)

    implicit none

    real :: A(:,:), matA(:,:)
    integer :: nx, ny, n, m

    do n = 1, nx
      do m = 1, ny
        matA(m, n) = A(m + ny*(n-1),1)
      end do 
    end do 
  end subroutine devectorize

end module NumDE




