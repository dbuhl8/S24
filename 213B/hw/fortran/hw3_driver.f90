! File: hw3_driver.f90
! Author: Dante Buhl
! Dependencies: NumDE, (Lapack and Lblas) or LinAl


program hw3_driver

  use LinAl 
  use NumDE

  implicit none

  ! Question 1
  integer, parameter :: nx1=81, ny1=51
  real, parameter :: xstart=0., xstop=2.
  real, parameter :: ystart=0., ystop=1.
  integer, parameter :: num_bound = 2*nx1 + 2*ny1 - 4
  integer, parameter :: num_dom = nx1*ny1
  integer, parameter :: num_int = num_dom-num_bound
  real, dimension(num_dom, num_dom) :: D2
  real, dimension(num_dom, 1) :: F, U, X, Y, G
  real, dimension(num_dom, num_int) :: Q_int, D2_int
  real, dimension(num_int, num_int) :: R_int
  real, dimension(num_int) :: Rvec ! used in QR solve
  real, dimension(num_int, 1) :: U_int, F_int
  real, dimension(ny1, nx1) :: matF, matU, matX, matY, matG
  integer, dimension(num_bound) :: bound
  integer, dimension(num_int) :: interior
  real :: dx = (xstop-xstart)/nx1
  real :: dy = (ystop-ystart)/ny1
  integer, dimension(num_dom) :: P
  
  ! Lapack variables
  real, dimension(num_int) :: S
  real :: RCOND
  real, allocatable :: work(:)
  integer :: lwork, info, RANK

  ! Question 2
  integer, parameter :: nx2=100, nt=10000
  real, dimension(nx2, nt+1) :: U2
  real, dimension(nx2, nx2) :: D
  real, dimension(nt) :: T
  real, dimension(nx2, 1) :: X2, F2
  real :: tstart, tstop
  real :: dt

  ! General Utility Vars
  real :: tol = 10.d-14
  logical :: bool
  integer :: i, j, n, m
  real :: start, finish

  print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 1: BVP with 2D Poison's Equation"
    print *, " "

    ! initializing Bounday Index arrays
    call vec_boundary(bound,interior, nx1, ny1)

    ! initializing D2 matrix
    call D2FD2_2D(D2, nx1, ny1, dx, dy)

    ! initalizing X, Y matrices.
    call vec_domain(X, Y, nx1, ny1, (/xstart, xstop/), (/ystart, ystop/))

    ! initializing vectorized F 
    F = -18.0 + 3.*X**2 + 4.*Y**2 - 8.*(pi**2)*(Y**2)*sin(pi*(Y**2)) + &
      4.*pi*cos(pi*(Y**2))
  
    ! initializing Boundary Condition
    G = 0.0
    G = 2. - X**2. - 2.*sin(pi*(Y**2.))

    call CPU_TIME(start)
    U = 0.0
    D2_int = D2(:, interior)
    allocate(work(1)) 
    lwork = -1
    info = 0
    !Calling a LinAl Solver to compute U_int

    ! DGELSS is an Lapack solver which computes the minimum norm solution to a
    ! least squares problem using SVD, QR iteration for an overdetermined system
    ! First call determines optimal size of work array 
    !call dgelss(num_dom, num_int, 1, D2_int, num_dom, F, num_dom, S, RCOND, &
                !RANK, work, lwork, info)
    !lwork = work(1)
    !deallocate(work)
    !allocate(work(lwork))
                 ! M         N    NRHS   A      LDA    B    LDB  
    !call dgelss(num_dom, num_int, 1, D2_int, num_dom, F, num_dom, S, RCOND, &
                !RANK, work, lwork, info)

    !call CPU_TIME(finish)

    !print *, finish - start

    ! Note that this is an overdetermined system
    !call householderQR(D2_int, Rvec, num_dom, num_int, bool, tol)
    !call formQ(D2_int, Q_int, num_dom, num_int) 
    !call formR(D2_int, R_int, Rvec, num_int) 

    !F_int = matmul(transpose(Q_int), F)
    !call backsub(R_int, F_int, U_int, num_int, 1) 

    !U = 0.0
    !U(interior, :) = F(1:num_int, :)

    !U = U + G

    !call devectorize(U, matU, nx1, ny1)
    !call devectorize(X, matX, nx1, ny1)
    !call devectorize(Y, matY, nx1, ny1)

    ! write U to a dat file so matlab can read it and plot the solution
    !open(10, file="u.dat") 
      !do i = 1, ny1
        !write(10, "("//trim(str(nx1))//"F30.15)") matU(i, :)
      !end do 
    !close(10)
    
    print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 2: IBVP with Heat Equation"
    print *, " "
    !Do it yeah

    tstart = 0
    tstop = 10
    dt = (tstop-tstart)/nt

    ! Initializing 1D D2 matrix using Finite Differences
    call D2FD2_1D(D, nx2, dx)

    dx = 2./(nx2-1)
    do i = 1, nx2
      X2(i,1) = -1 + (i-1)*dx
    end do 
     
    ! Initializing F vector 
    F2 = -20. + 60.*X2**2 

    ! Setting Initial Condition 
    U2 = 0.0
  
    call  IBVP_1DCN(U2, D, F2, T, nx2, nt, dt)

    F2 = (3 + X2) + 5*(1 - X2**2)**2

    ! write to out file
    open(11, file="u2.dat")
      do i = 1, nt+1
        write(11, "("//trim(str(nx2))//"F10.4)") U2(:,i) - F2(:,1)
      end do 
    close(11)

    print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
 
end program hw3_driver
