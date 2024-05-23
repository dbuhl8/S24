!File: hw5.f90
!Author: Dante Buhl
!Dependencies: OpenMP, LinAl

program driver

  use omp_lib
  use LinAl

  implicit none
  integer, parameter :: kr=kind(dble(1.0))
  integer :: id, np, i, j
  integer, parameter :: nd=10000
  real(kind=kr), allocatable :: MA(:,:), MB(:,:), MC(:,:), MC2(:,:)
  real(kind=kr) :: start, finish, ratio1, int_sum1, int_sum2, ratio2

  ! int vars
  real(kind=kr) :: xstop, xstart, dx, temp_sum, min_ws, error_sum
  integer :: n
  integer :: num_threads = 8
 
  allocate(MA(nd,nd), MB(nd,nd), MC(nd,nd),MC2(nd,nd))
  call ident(MA, nd)
  do j = 1, nd
    do i = 1, nd
      MB(i, j) = i+j
    end do 
  end do 
  int_sum1 = 0.0
  int_sum2 = 0.0 
  xstart = 0.0
  xstop = 2.0
  temp_sum = 0.0
  n = 100000000
  dx = (xstop-xstart)/n 
  call cpu_time(start)
  do i = 1, n
    int_sum1 = int_sum1 + dx*((xstart + dx*(i-1))**2 + ((xstart + dx*i)**2))/2
  end do
  call cpu_time(finish)
  ratio1 = finish-start

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,id) 
    ! Question 1
    id = OMP_GET_THREAD_NUM()
    np = OMP_GET_NUM_THREADS()
    ! Print Hello from each processor
    print "(A, I3, A, I3)","Hello from processor #",id," out of ", np

    ! Question 2
    call cpu_time(start)
    !$OMP DO REDUCTION(+:temp_sum)
      do i = 1, n
        temp_sum = temp_sum + dx*((xstart+dx*(i-1))**2+((xstart + dx*i)**2))/2
      end do
    !$OMP END DO
    call cpu_time(finish)
    ratio1 = ratio1/(finish-start)
  
    ! Question 3
    !$OMP DO 
      do j = 1, nd
        do i = 1, nd
          MC2(i,j) = sum(MA(i,:)*MB(:,j))
        end do 
      end do 
    !$OMP END DO
    !$OMP WORKSHARE 
      MC = matmul(MA, MB)
    !$OMP END WORKSHARE
    !$OMP DO REDUCTION(+:error_sum)
      do j = 1, nd
        do i = 1, nd
          error_sum = error_sum + abs(MC(i,j) - MC2(i,j))
        end do
      end do
    !$OMP END DO
  !$OMP END PARALLEL
  print *, " "
  print "(2(A, F10.5))", "Parallel Sum : ", temp_sum, ", Sequential Sum : ",&
    int_sum1
  print *, " "
  print *, "Error in matmul mathods", error_sum
  print *, " "
   
  deallocate(MA, MB, MC)

end program driver
