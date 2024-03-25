



program piApprox
implicit none
include 'mpif.h'

! Variable Declarations
integer ie, id, np, i, pisum
integer cpu0, counts
integer stat(MPI_STATUS_SIZE)
real pi
!integer, intent(in) :: n
real, allocatable ::  data(:, :), coords(:, :), vals(:)
integer, allocatable :: data2(:)

!print *, "How many darts would you like to throw? Please enter an integer:"
!read (*, *) n

cpu0 = 0
pisum = 0

! MPI Init
call MPI_INIT(ie)
call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)
call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)

if (id == cpu0) then
    allocate(data(10000, 2))
    call random_number(data)
end if
counts = FLOOR(10000./np)
allocate(coords(counts, 2))
allocate(vals(counts))
allocate(data2(np))


! Scatter Coordinate Data to processors
call MPI_SCATTER(data(:,1), counts, MPI_REAL, coords(:,1), counts, MPI_REAL, cpu0, MPI_COMM_WORLD, ie)
call MPI_SCATTER(data(:,2), counts, MPI_REAL, coords(:,2), counts, MPI_REAL, cpu0, MPI_COMM_WORLD, ie)

! Converting coordinates to radius values
vals(:) = (coords(:, 1)-.5)**2 + (coords(:, 2)-.5)**2

! Checks if they are in the circle
do i = 1, counts
    if (vals(i) .le. .25) then 
        !print *, "ID:",id,"Point:",vals(i)
        pisum = pisum + 1
    end if
end do


! Gathers all of the data back to processor 1!
call MPI_GATHER(pisum, 1, MPI_INTEGER, data2, 1, MPI_INTEGER, cpu0, MPI_COMM_WORLD, ie )

if (id == cpu0) then
    ! Converts to probability and Mulitplies by 4!
    pi = (sum(data2)*4.)/FLOAT(np*counts)
    print *, "The calculated approximation for Pi was:", pi
end if
! MPI End
call MPI_FINALIZE(ie)

stop

end program piApprox



