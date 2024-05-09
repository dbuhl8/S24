

program sendrec
  include 'mpif.h'

  integer :: id, ie, np
  integer :: msgid, src, dest, count
  integer :: buffer 
  integer :: stat(MPI_STATUS_SIZE)
  character(9) :: signal

  ! Starting MPI
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)

  ! Coordinating Send and Rec Address
  msgid = 10
  src = 0
  dest = 1
  count = 9

  ! Branched Logic to decide who sends and who receives
  if (id .eq. src) then
    ! Only the source processor will have this variable in memory before the
    ! send
    signal = 'Send this'
    call MPI_SEND(signal,count,MPI_CHARACTER,dest,msgid,MPI_COMM_WORLD,ie)
    print *, 'Processor: ', id, ' sent     :', signal
  else if (id .eq. dest) then 
    call MPI_RECV(signal,count,MPI_CHARACTER,src,msgid,MPI_COMM_WORLD,stat,ie)
    print *, 'Processor: ', id, ' received :', signal
  end if

  ! Stop MPI
  call MPI_FINALIZE(ie)

end program sendrec



