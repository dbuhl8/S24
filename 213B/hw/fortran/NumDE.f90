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
    real, dimension(ma, 3) :: K

    do i = 1, num_points

        K = 0.0
    
        ! update k vectors
        do j = 1, 3
            K(:, j:j) = matmul(F, Y(:, i:i) + dt*matmul(K, transpose(A(j:j,:))))
        end do

        ! update y
        Y(:, i+1) = Y(:, i)
        do j = 1, 3
            Y(:, i+1) = Y(:, i+1) + dt*b(j)*K(:, j)
        end do

        T(i+1) = T(i) + dt

    end do
  end subroutine RK3

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

end module NumDE




