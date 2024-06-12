!File: final_driver.f90
!Author: Dante Buhl
!Dependencies: NumDE, LinAl

program final_driver

  use NumDE
  use LinAl

  implicit none

  ! Gen Vars
  integer, parameter :: kr=kind(dble(1.0))
  integer :: i, j, k, n, fn=10
  real(kind=kr) :: start, finish
  real(kind=kr) :: sub_start, sub_finish

  ! LAPACK Vars
  integer :: lwork, info
  real(kind=kr), allocatable :: work(:), WR(:), WI(:)
  real(kind=kr),dimension(1,1) :: rdum
  integer :: idum=1

  ! Question 1
  integer, parameter :: nb=60, nq=50000
  real(kind=kr), dimension(nb,1) :: beta,weights
  real(kind=kr) :: err,wnm, b2,e2,de
  real(kind=kr) :: l, dl, s,xs
  
  ! Question 2
  real(kind=kr), parameter :: dt = 5*1.d-5
  integer, parameter :: nt = int(2./dt), num_nr=4
  real(kind=kr), parameter :: r_start = 1.d0, r_stop = 3.d0
  integer, dimension(num_nr) :: nr_vals = (/50,100,150,200/)
  integer :: nr
  real(kind=kr) :: dr, v=0.25
  real(kind=kr), dimension(nt+1) :: TS
  real(kind=kr), allocatable :: U(:,:), R(:,:), T(:), UA(:,:), UAS(:,:),RM(:,:)
  real(kind=kr), allocatable :: D1(:,:), D2(:,:),D(:,:)
  real(kind=kr), dimension(20) :: spect, dt_crit
  integer, dimension(20) :: n_vals
  real(kind=kr), dimension(num_nr,nt+1) :: mpe_error


  ! Nonlinear solver for B values
  ! first three are given to us
  beta(1,1) = 1.548458778289446
  beta(2,1) = 3.129084015718067
  beta(3,1) = 4.703797206969750
  wnm = .75

  call cpu_time(start)
  ! for each value of B needed 
    do i = 4, nb
      ! proceed by a newtons method solve for tolerance 
      ! initial guess
      beta(i,1) = 2*beta(i-1,1) - beta(i-2,1)
      b2 = beta(i,1) - 0.25*(beta(i,1) - beta(i-1,1))
      !compute first error
      err = h_bes(beta(i,1),r_start,r_stop)
      e2 = h_bes(b2,r_start,r_stop) 
      de = (err - e2)/(beta(i,1)-b2)
      j = 1
      do while (err .gt. 1.d-15) 
        ! update Beta 
        b2 = beta(i,1) - wnm*err/de
        de = -beta(i,1)
        beta(i,1) = b2
        de = beta(i,1) + de
        ! compute error again
        e2= err
        err = h_bes(beta(i,1),r_start,r_stop)
        de = (err - e2)/de
        j = j+1
      end do 
    end do 
  ! -------------------------------------------
  call cpu_time(finish)
  print *, "Computing B values took ", finish-start, " seconds."

  call cpu_time(start)
  ! Compute Integral with a Quadrature Rule
    dl = (r_stop-r_start)/(nq-1)
    do i = 1, nb
      s = 0
      do j = 1, nq-1
        xs = r_start + dl*(j-1)
        s = s + xs*h_r(beta(i,1),xs,r_start,r_stop)*u_sol(xs,r_start,r_stop)+&
          (xs+dl)*h_r(beta(i,1),xs+dl,r_start,r_stop)*u_sol(xs+dl,r_start,r_stop)
        ! s = (y0 + 2y1 + 2y2 + ... + yN)*dx/2
      end do
      weights(i,1) = s*dl/2
    end do 
  ! -------------------------------------------
  call cpu_time(finish)
  print *, "Computing Integrals took ", finish-start, " seconds."

  call cpu_time(start)
  ! Writing out R basis functions
    open(fn, file='r.dat')
      dl = 2./200.
      do i = 1, 201
        xs = r_start + (i-1)*dl
        write (fn, '(5F20.8)') xs, h_r(beta(1,1),xs,r_start,r_stop),&
          h_r(beta(2,1),xs,r_start,r_stop),&
          h_r(beta(4,1),xs,r_start,r_stop),&
          h_r(beta(6,1),xs,r_start,r_stop)
      end do 
    close(fn)
    fn = fn+1
    open(fn, file='r_norm.dat')
      do i = 1, 30
        write(fn, '(I6, F20.16)') i, beta(i,1)
      end do 
    close(fn)
    fn = fn+1
  ! -------------------------------------------
  call cpu_time(finish)
  print *, "Basis Function Output took ", finish-start, " seconds."

  call cpu_time(start)
  ! Critical dt study ------------------------
    open(fn, file='dt.dat')
    do i = 1, 20
      n = 10 + 10*(i-1)
      n_vals(i) = n
      allocate(D(n,n),D1(n,n),D2(n,n))
      dr = (r_stop-r_start)/(n-1)
      call D1FD2_1D(D1, n, dr)
      call D2FD2_1D(D2, n, dr)
      do j = 1, nr
        D1(i,:) = D1(i,:)/R(i+1,1)
      end do 
      D = v*(D2 + D1)
      ! lapack call
      lwork = -1 
      allocate(work(1),WR(n),WI(n))
      call DGEEV('N','N',n,D,n,WR,WI,rdum,idum,rdum,idum,WORK,LWORK,info)
      lwork = work(1)
      deallocate(work)
      allocate(work(lwork))
      call DGEEV('N','N',n,D,n,WR,WI,rdum,idum,rdum,idum,WORK,LWORK,info)

      if(info .ne. 0) then
        print *, "Error in Lapack Call"
      end if

      spect(i) = maxval(sqrt(WR**2+WI**2))
      dt_crit(i) = 6.d0/(11.d0*spect(i))
      deallocate(D,D1,D2)
      deallocate(work,WR,WI)
      write(fn, '(I6,F20.11,F20.15)') n_vals(i), spect(i), dt_crit(i)
    end do 
    close(fn) 
    fn = fn + 1
  ! -------------------------------------------
  call cpu_time(finish)
  print *, "DT val study took ", finish-start, " seconds."

  call cpu_time(start)
  ! Order of convergence study 
    do i = 1, num_nr
      nr = nr_vals(i)
      dr = (r_stop-r_start)/(nr+1)
      allocate(U(nr+2,nt+1),R(nr+2,1),T(nt+1))
      allocate(D1(nr,nr),D2(nr,nr),D(nr,nr))
      allocate(UA(nr+2,1), RM(nr+2,nb))
      if (nr.eq.200) then
        allocate(UAS(nr+2,nt+1))
      end if
      do j = 1, nr+2
        R(j,1) = r_start + (j-1)*dr
      end do 
      U(:,1) = 15*((R(:,1)-r_start)**2)*((r_stop-R(:,1))**2)&
        *exp(-sin(2*R(:,1))-R(:,1))
      
      call D1FD2_1D(D1, nr, dr)
      call D2FD2_1D(D2, nr, dr)
     
      do j = 1, nr
        D1(j,:) = D1(j,:)/R(j+1,1)
      end do 
       
      D = v*(D2+D1) 
      call cpu_time(sub_start)
      call AB3(D, U(2:nr+1,:), T, nr, nt, dt)

      U(1,2:) = U(1,1)
      U(nr+2,2:) = U(nr+2,1)
      call cpu_time(sub_finish)
      print *, " Numerical Integral took ", sub_finish-sub_start, " seconds"

      ! compute max pointwise error for each time
      call cpu_time(sub_start)
      do j = 1, nb
        RM(:,j) = h_rv(beta(j,1),R,r_start,r_stop,nr+2)*weights(j,1)
      end do 
      do j = 1, nt+1
        UA(:,1) = sum(matmul(RM,exp(-(beta**2)*v*T(j))),dim=2)
        if (nr .eq.200) then
          UAS(:,j) = UA(:,1)
        end if
        mpe_error(i,j) = maxval(abs(U(2:nr+1,j) - UA(2:nr+1,1)))
      end do 
      call cpu_time(sub_finish)
      print *, " Analytical Solution took ",sub_finish-sub_start," seconds"

      ! Write Output
      if (nr .eq. 200) then
        open(fn,file='U.dat')
          call writemat(U, nr+2, nt+1, fn)
        close(fn)
        fn = fn + 1
        open(fn,file='UA.dat')
          call writemat(UAS, nr+2, nt+1, fn)
        close(fn)
        fn = fn + 1
        open(fn,file='T.dat')
          do j = 1, nt+1
            write(fn, '(F10.6)') T(j)
          end do 
        close(fn)
        fn = fn + 1
        print *, "nt = ", nt+1
      end if
      TS = T
      deallocate(D1,D2,D,U,R,T,UA,RM)
      if(nr.eq.200)then
        deallocate(UAS)
      end if
      print *, "Finished nr = ", nr
    end do 
  ! -------------------------------------------
  call cpu_time(finish)
  print *, "Grid Size study took ", finish-start, " seconds."

  call cpu_time(start)
  ! Writing Output files
    open(fn, file='mpe.dat')
      do i = 1, nt+1
        write(fn, '('//trim(str(num_nr+1))//'F32.16)') TS(i), mpe_error(:,i)
      end do 
    close(fn)
    fn = fn + 1
    open(fn, file='final_mpe.dat')
      do i = 1, num_nr
        write(fn, "(I6, F20.8)") nr_vals(i), mpe_error(i,nt+1)
      end do 
    close(fn)
  ! -------------------------------------------
  call cpu_time(finish)
  print *, "Writing OUTPUT files took ", finish-start, " seconds."

  contains
  ! Extra Functions and whatnot
    function h_bes(b, r1, r2)
      implicit none
      real(kind=kr) :: b, r1, r2, h_bes
      h_bes = bessel_j0(b*r1)*bessel_y0(b*r2)-bessel_j0(b*r2)*bessel_y0(b*r1)
    end function h_bes

    function h_hat(b, x, r1, r2)
      implicit none
      real(kind=kr) :: b, x ,r1, r2, h_hat
      h_hat = bessel_j0(b*x)*bessel_y0(b*r2)-bessel_j0(b*r2)*bessel_y0(b*x)
    end function h_hat

    function h_norm(b,r1,r2)
      implicit none
      real(kind=kr) :: b, r1, r2, h_norm
      h_norm = (2./((pi**2)*(b**2)))*((bessel_j0(b*r1))**2-(bessel_j0(b*r2))**2)/&
        ((bessel_j0(b*r1))**2)
      !h_norm = h_hat(b,x,r1,r2)/sqrt(h_norm)
    end function h_norm

    function h_r(b,x,r1,r2)
      implicit none
      real(kind=kr) :: b, r1, r2, h_r,x
      h_r = h_hat(b,x,r1,r2)/sqrt(h_norm(b,r1,r2))
    end function h_r
    function h_rv(b,x,r1,r2,np)
      implicit none
      real(kind=kr) :: b, r1, r2,x(:,:)
      integer :: np, si
      real(kind=kr), dimension(np) :: h_rv
      do si = 1, np
        h_rv(si) = h_hat(b,x(si,1),r1,r2)/sqrt(h_norm(b,r1,r2))
      end do 
    end function h_rv
    function u_sol(x,r1,r2)
      implicit none
      real(kind=kr) u_sol, x,r1,r2
      u_sol = 15*((x-r1)**2)*((r2-x)**2)*exp(-sin(2*x)-x)
    end function u_sol
    function u_solv(x,r1,r2,np)
      implicit none
      real(kind=kr) :: x(:,:),r1,r2
      integer :: np
      real(kind=kr),dimension(np,1) :: u_solv
      u_solv = 15*((x-r1)**2)*((r2-x)**2)*exp(-sin(2*x)-x)
    end function u_solv
  ! end of extra functions subroutines

end program final_driver


