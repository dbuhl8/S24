!File: Numerical Differential Equations Module
!Author: Dante Buhl
!Purpose: To store numerical methods for Differential Equations

module NumDE

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




end module NumDE




