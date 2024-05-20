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
  real, dimension(num_dom, 1) :: F, U, X, Y, G, G2
  real, dimension(num_dom-4, num_int) :: D2_int
  real, dimension(num_dom-4, 1) :: F_int
  real, dimension(ny1, nx1) :: matF, matU, matX, matY, matG
  integer, dimension(4) :: corner
  integer, dimension(num_dom-4) :: nc
  integer, dimension(num_bound) :: bound
  integer, dimension(num_int) :: interior
  real :: dx = (xstop-xstart)/nx1
  real :: dy = (ystop-ystart)/ny1
  real :: norm=0, norm_2=0
  integer, dimension(num_dom) :: P
 
  ! Lapack variables
  real, dimension(num_int) :: S
  real :: RCOND
  real, allocatable :: work(:)
  integer :: lwork, info, RANK

  ! Question 2
  real, parameter :: dt=10.d-4
  real,parameter :: tstart=0, tstop=2
  integer, parameter :: nx2=100, nt=int((tstop-tstart)/dt), na = 100
  integer, parameter :: num_modes=64
  real, dimension(nx2, nx2) :: D
  real, dimension(nx2*na, 1) :: UA, XA, TA
  real, dimension(nx2, na) :: UA_mat, XA_mat, TA_mat
  real, dimension(nx2, 1) :: X2, F2
  real, dimension(num_modes) :: C
  real, allocatable :: UF(:,:), US(:,:)
  real, allocatable :: DF(:,:), DS(:,:)
  real, allocatable :: TF(:), TS(:)
  real, allocatable :: XF(:,:), FF(:,:), XS(:,:), FS(:,:)
  real, allocatable :: XXF(:,:), TTF(:,:), XXS(:,:), TTS(:,:)
  integer :: nf, ns 


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
    call vec_boundary(bound,interior,corner,nc,nx1,ny1)

    ! initializing D2 matrix
    call D2FD2_2D(D2, nx1, ny1, dx, dy)

    ! initalizing X, Y matrices.
    call vec_domain(X, Y, nx1, ny1, (/xstart, xstop/), (/ystart, ystop/))

    ! initializing Boundary Condition
    G = 2. - X**2. - 2.*sin(pi*(Y**2.))
    G2 = -2.0 + 8.*(pi**2)*(Y**2)*sin(pi*(y**2)) - 4.*pi*cos(pi*(y**2)) 

    ! initializing vectorized F 
    F = -20.0 + 3.*(X**2) + 4.*(Y**2)
    F = F - G2

    ! Preparing Matrices for Least Squares Solve
    F_int(:, 1) = F(nc,1)
    D2_int = D2(nc, interior) 
    allocate(work(1)) 
    lwork = -1
    info = 0

    ! DGELS is an Lapack solver which computes a least squares solve

    ! D2_int * U_int = F

    ! First call determines optimal size of work array 
                !! T     M         N    NRHS   A       LDA,   B,      LDB
    !call dgels('N', num_dom-4, num_int, 1, D2_int, num_dom, F_int, num_dom-4,&
      !work, lwork, info)
                ! T     M         N    NRHS   A       LDA,   B,      LDB
    call dgels('N', num_dom, num_dom, 1, D2, num_dom, F, num_dom,&
      work, lwork, info)

    ! Work is now the optimal size
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    ! Seeing how efficient this solve is
    call CPU_TIME(start)

                !! T     M         N    NRHS   A       LDA,   B,      LDB
    !call dgels('N', num_dom-4, num_int, 1, D2_int, num_dom, F_int, num_dom-4,&
      !work, lwork, info)
                 ! T     M         N    NRHS   A       LDA,   B,      LDB
    call dgels('N', num_dom, num_dom, 1, D2, num_dom, F, num_dom,&
      work, lwork, info)
   
    call CPU_TIME(finish)
    if (info .ne. 0) then
      print *, "Error in Lapack Routine"
    end if
    deallocate(work)

    print "(A, F10.4, A)", " Least Squares Solve finished in ", finish - start,&
      " seconds"

    ! Computing the 2-norm of the error
    !U = 0.0
    !F(nc, :) = F_int
    !U(interior, :) = F(interior, :)
    U = F
    U = U + G
    F = -20.0 + 3.*X**2 + 4.*Y**2 
    call D2FD2_2D(D2, nx1, ny1, dx, dy)
    F = matmul(D2, U) - F
    call twonorm(F(:,1), norm)
    print *, "Solution Error: ", norm
    F = -20.0 + 3.*X**2 + 4.*Y**2 

    ! Retuning U as a matrix
    call devectorize(U, matU, nx1, ny1)
    call devectorize(X, matX, nx1, ny1)
    call devectorize(Y, matY, nx1, ny1)
    call devectorize(F, matF, nx1, ny1)
    call devectorize(G, matG, nx1, ny1)

    !write U to a dat file so matlab can read it and plot the solution
    open(10, file="u.dat") 
      do i = 1, ny1
        write(10, "("//trim(str(nx1))//"F30.15)") matU(i, :)
      end do 
    close(10)
    open(11, file="x.dat") 
      do i = 1, ny1
        write(11, "("//trim(str(nx1))//"F30.15)") matX(i, :)
      end do 
    close(11)
    open(12, file="y.dat") 
      do i = 1, ny1
        write(12, "("//trim(str(nx1))//"F30.15)") matY(i, :)
      end do 
    close(12)
    open(13, file="f.dat") 
      do i = 1, ny1
        write(13, "("//trim(str(nx1))//"F30.15)") matF(i, :)
      end do 
    close(13)
    open(14, file="g.dat") 
      do i = 1, ny1
        write(14, "("//trim(str(nx1))//"F30.15)") matG(i, :)
      end do 
    close(14)

    
    print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 2: IBVP with 1D Heat Equation"
    print *, " "
    !Do it yeah

    ! Computing Analytical Solution
    call vec_domain(XA, TA, nx2, na, (/-1.0, 1.0/), (/0.0, 2.0/))

    C(1) = 10.*(1 - 2./3. + 1./5.)/2.
    XA = XA + 1
    do i = 1, num_modes
      C(i) = 5*((-1.)**i-1.)*((64./((pi*i)**3))-(768./((pi*i)**5)))
      UA = UA + C(i)*sin((i*pi/2.)*XA)*exp(-TA*(i*pi/2.0)**2)
    end do 
    XA = XA - 1 
    UA = UA + 3. + XA
    call devectorize(UA, UA_mat, nx2, na)
    call devectorize(XA, XA_mat, nx2, na)
    call devectorize(TA, TA_mat, nx2, na)


    !Computing Numerical Solution using Finite Differences and Crank Nicolson
    nf = 100
    !Allocating memory for the Finite Difference Method
    allocate(UF(nf, nt+1), DF(nf, nf), XF(nf, 1), FF(nf, 1), TF(nt+1))

    !Allocating memory for the GCL Spectral Method
    allocate(US(nf+1, nt+1), DS(nf+1, nf+1), XS(nf+1, 1), FS(nf+1, 1), TS(nt+1))
    
    ! Allocating Memory for the Output
    allocate(XXS(nf+1, nt+1), TTS(nf+1,nt+1))
    allocate(XXF(nf, nt+1), TTF(nf,nt+1))

    ! Initializing Domain for the spectral method
    call GCLgrid_1D(XS, nf)
    call GCLmeshgrid(XXS, TTS, nf, nt, dt, 0.)

    ! Initializing differentiation matrices for each method
    call D2FD2_1D(DF, nf, dx)
    call D2GCL_1D(DS, nf)

    ! Printing out D2 matrix to check it if is correct
    !call printmat(DS, nf+1, nf+1)
    DS(1,:) = 0.
    DS(nf+1,:) = 0.
    DF(1,:) = 0.
    DF(nf, :) = 0.

    ! Initializing domain for Finite Difference Methods
    dx = 2./(nf-1)
    do i = 1, nf
      XF(i,1) = -1 + (i-1)*dx
    end do 

    ! Setting Initial Condition 
    UF = 0.0
    UF(:,1:1) = 5*(1 - XF**2)**2
    US = 0.0
    US(:,1:1) = 5*(1 - XS**2)**2

    ! Integrating for each method 
    call  IBVP_1DCN(UF, DF, TF, nf, nt, dt)
    call  IBVP_1DCN(US, DS, TS, nf+1, nt, dt)

    ! Initializing Meshgrid for the Finite Difference Output
    call meshgrid(TTF, XXF, TF, XF(:,1), nt+1, nf)

    ! Converging back into Physical Space
    ! There is a massive bug in this routine. I'll fix it tomorrow. 
    !do i = 1, nt
      !call GCL_2_phys(US(:, i+1), nf)
    !end do

    ! Initializing Shift Vector
    FF = (3 + XF)
    FS = (3 + XS)

    ! write to out file
    open(15, file="u2.dat")
      do i = 1, na
        write(15, "("//trim(str(nx2))//"F10.4)") UA_mat(i,:)
      end do 
    close(15)
    open(16, file="x2.dat")
      do i = 1, na
        write(16, "("//trim(str(nx2))//"F10.4)") XA_mat(i,:)
      end do 
    close(16)
    open(17, file="t2.dat")
      do i = 1, na
        write(17, "("//trim(str(nx2))//"F10.4)") TA_mat(i,:)
      end do 
    close(17)
    ! write to out file
    open(15, file="us.dat")
      do i = 1, nf+1
        write(15, "("//trim(str(nt+1))//"F10.4)") US(i,:)+FS(i,1)
      end do 
    close(15)
    open(16, file="xs.dat")
      do i = 1, nf+1
        write(16, "("//trim(str(nt+1))//"F10.4)") XXS(i,:)
      end do 
    close(16)
    open(17, file="ts.dat")
      do i = 1, nf+1 
        write(17, "("//trim(str(nt+1))//"F10.4)") TTS(i,:)
      end do 
    close(17)

    open(18, file="uf.dat")
      do i = 1, nf
        write(18, "("//trim(str(nt+1))//"F10.4)") UF(i,:)+FF(i,1)
      end do 
    close(18)
    open(19, file="xf.dat")
      do i = 1, nf
        write(19, "("//trim(str(nt+1))//"F10.4)") XXF(i,:)
      end do 
    close(19)
    open(20, file="tf.dat")
      do i = 1, nf
        write(20, "("//trim(str(nt+1))//"F10.4)") TTF(i,:)
      end do 
    close(20)

    deallocate(UF, DF, FF, XF, TF)
    deallocate(US, DS, FS, XS, TS)
    deallocate(XXF, XXS, TTF, TTS) 

    print *, " "
  print *, "-------------------------------------------------------------------"
  print *, " "
  print *, " Question 3: Extra Credit"
    print *, " "


    print *, " "
  print *, "-------------------------------------------------------------------"
 
end program hw3_driver
