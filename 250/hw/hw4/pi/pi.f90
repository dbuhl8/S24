!File: pi.f90
!Author: Dante Buhl
!Dependencies: MPI

program piApprox
  implicit none
  include 'mpif.h'

  ! Variable Declarations
  integer, parameter :: kr=kind(dble(1.0)), nco=1000000
  integer :: ie, id, np, i, pisum, sum_tot
  integer :: cpu0, seed_size, clock
  integer :: stat(MPI_STATUS_SIZE)
  integer, allocatable :: seed(:)
  real(kind=kr) :: pi
  real(kind=kr), dimension(nco, 2) :: co
  
  cpu0 = 0
  pisum = 0

  ! MPI Init
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)

  !calls random num within 0, 1
  call random_seed(size=seed_size)
  print *, "Seed size: ", seed_size
  allocate(seed(seed_size))
  call system_clock(count=clock)
  seed = clock + 37*(/ (i-1, i = 1, seed_size) /) + id
  call random_seed(put=seed)
  call random_number(co)
  !centers the data about origin, [-0.5, 0.5]
  co = co - 0.5
  ! finds the distance squared from origin
  co(:, 1) = co(:,1)**2 + co(:,2)**2
  ! Note that the radius of the circle considered is 0.5, so the radius squared
  ! will be 0.25

  ! Checks if they are in/on the circle
  do i = 1, nco
    if (co(i,1) .le. .25) then 
      pisum = pisum + 1
    end if
  end do

  pi = (pisum*4.)/(nco)
  print *, "Proc: ", id, ", The calculated approximation for Pi was:", pi

  sum_tot = 0
  ! Gathers the sum of pisum across all processors into tot_sum in processor 0
  ! Arguments: sendbuf,recvbuf,count,  type ,operation,root,   comm,      err
  call MPI_REDUCE(pisum, sum_tot,1,MPI_INTEGER,MPI_SUM,cpu0,MPI_COMM_WORLD,ie)

  ! On root processor, compute pi and print it out
  if (id .eq. cpu0) then
    ! Converts to probability and Mulitplies by 4!
    pi = (sum_tot*4.)/(np*nco)
    print *, "The collective calculated approximation for Pi was:", pi
  end if

  ! Stop MPI
  call MPI_FINALIZE(ie)

end program piApprox



