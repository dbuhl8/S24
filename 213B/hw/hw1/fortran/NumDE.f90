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
    ! F should be a vector containing at num_points + 3 entries such that each entry is
    ! F(i) = f(i), where f(x) is any function we want to differentiate. The first 3 points should be padding
    ! points, i.e. points used to calculate the first derivative points at F(4)
    ! F(1:num_points+3) = f(X(1:num_points+3)), where X(4) = x_0
    
    ! dx should be the step size between adjecent entries of.

    ! on Output: 
    ! F(1:num_points) are the values of the derivative df/dx according to the BDF3 Method

    implicit none

    real :: F(:), dx
    integer :: i, num_points

    do i = 1, num_points+1
        F(i) = (11.*F(i+3) - 18.*F(i+2) + 9.*F(i+1) - 2.*F(i))/(6.*dx)
    end do

  end subroutine BDF3 

  subroutine aeRK3(F, Y, A, B, ma, num_points, dt)
    ! autonomous, explicit
    ! on entry F, Y, are 2-d arrays

    implicit none
                        
    real :: Y(:, :), F(:, :), A(:, :), B(:)
    integer, intent(in) :: ma, num_points
    integer :: i, j, k
    real, dimension(ma, ma) :: K

    do i = 1, num_points

        K = 0.0
    
        ! update k vectors
        do j = 1, ma
            K(:, i) = matmul(F, Y(:, i) + dt*matmul(A(:, i:i), K))
        end do

        ! update y
        Y(:, i+1) = Y(:, i)
        do j = 1, ma
            Y(:, i+1) = Y(:, i+1) + dt*b(j)*K(:, j)
        end do

    end do

  end subroutine aeRK3

  subroutine AM3(F, Y, ma, num_points, dt, num_iter)
    ! this subroutine can only handle linear ODE's as of right now
    ! That ism the jacobian matrix, A, is constant

    implicit none

    real :: Y(:, :), A(:, :), dt
    integer, intent(in) :: ma, num_points, num_iter
    real, dimension(ma, ma) :: Jac, InvJac
    real, dimension(ma, 1) :: B, X
    integer, dimension(ma) :: P
    integer :: i, j
    logical :: bool


    ! A = I - J_G(u_k+1)
    call eye(Jac, ma)
    Jac = Jac - F
    InvJac = Jac
    call LU(InvJac, ma, bool, P)

    do j = 1, num_points

        ! implicit step, Fixed Point Iteration for u_k+1
        B = Y(:, i:i) + dt*matmul(
        X = 0.0
        do i = 1, num_iter
            ! find u_k+1^j+1 = u_k+1^j + (eye - dt*A)^-1(u_k + dt A(u_k+1^j) - u_k+1^j)
            call LUsolve(InvJac, ma, B, X, 1, P)
            B = Y(:, i:i) + dt*matmul(F, ) - Y(:, i+1:i+1) + matmul(Jac, Y(:, i:i))
            Y(:, i+1) = X
        end do
    `        
        ! updating y vectors
        Y(:, i) = !AM3 Method

    end do

  end subroutine AM3

end module NumDE




