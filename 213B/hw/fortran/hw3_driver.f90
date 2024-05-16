

program hw3_driver

  use LinAl
  use NumDE

  implicit none

  ! problem 1
  integer, parameter :: nx1=4, ny1=4
  real, parameter :: xstart=0., xstop=2.
  real, parameter :: ystart=0., ystop=1.
  real, dimension(nx1*ny1, nx1*ny1) :: D2
  real, dimension(nx1*ny1, 1) :: F, U, X, Y
  real, dimension(nx1, ny1) :: matF, matU, matX, matY
  integer, parameter :: num_bound = 2*nx1 + 2*ny1 - 4
  integer, parameter :: num_dom = nx1*ny1
  integer, dimension(num_bound) :: bound
  real :: dx = (xstop-xstart)/nx1
  real :: dy = (ystop-ystart)/ny1
  logical :: bool
  integer, dimension(num_dom) :: P

  ! General Utility Vars
  integer :: i, j, n, m


  print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 1: BVP with 2D Poison's Equation"
    print *, " "

    ! initializing Bounday Index arrays
    call vec_boundary(bound, nx1, ny1)

    ! initializing D2 matrix
    call D2FD2(D2, nx1, ny1, dx, dy)


    ! initializing vectorized U   
    !set boundaries to 0, inner domain to 1
    !U = 1.0
    !do i = 1, num_bound
      !U(bound(i)) = 0.0
    !end do 

    ! initalizing X, Y matrices.
    call vec_domain(X, Y, nx1, ny1, (/xstart, xstop/), (/ystart, ystop/))

    ! initializing vectorized F 
    F = -18.0 + 3*X**2 + 4*Y**2 - 8*(pi**2)*(Y**2)*sin(pi*(Y**2)) + &
      4*pi*cos(pi*(Y**2))
    ! Note: that F doesn't change per each iteration
    ! Note: F is purely a function of x, y not u
 
    !call printmat(D2, num_dom, num_dom)
    !print *, "U: "
    !print "("//trim(str(num_dom))//"F10.4)", U(:)
    !print *, "X: "
    !print "("//trim(str(num_dom))//"F10.4)", X(:)
    !print *, "Y: "
    !print "("//trim(str(num_dom))//"F10.4)", Y(:)
    !print *, "F: "
    !print "("//trim(str(num_dom))//"F10.4)", F(:)

    !Calling a LinAl Solver to compute u
    U = 0.0
    call LU(D2, num_dom, bool, P)
    call LUsolve(D2, num_dom, F, U, 1, P)

    call devectorize(U, matU, nx1, ny1)
    call devectorize(X, matX, nx1, ny1)
    call devectorize(Y, matY, nx1, ny1)

    call printmat(matU, nx1, ny1)
    call printmat(matX, nx1, ny1)
    call printmat(matY, nx1, ny1)

    ! write U, X, Y, out
    
    print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 2: IBVP with Heat Equation"
    print *, " "
 
end program hw3_driver
