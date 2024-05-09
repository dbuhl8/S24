!File: ring.f90
!Author: Dante Buhl
!Dependencies: MPI

program ring
  implicit none

  include 'mpif.h'

  integer ie, id, np, rid, i
  integer msgid, cpu0, cpuf
  integer stat(MPI_STATUS_SIZE)
  integer, allocatable :: a(:)

  !Starting MPI
  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)

  !Getting Processor Info
  cpu0 = 0
  cpuf = np-1
  msgid = id
  allocate(a(np+1))

  !initialzing values in a(1)
  do i = 1, np
    ! Branched Logic
    if(id .eq. i-1) then
      ! Each CPU begins with a(1) = processor id
      a(1) = i-1
    end if
  end do

  ! sending to the right
  do i = 1, np
    if (id == cpu0) then
      ! Edge case, must receive from the other end
      call MPI_SEND(a(i), 1, MPI_INTEGER, id+1, msgid, MPI_COMM_WORLD, ie)
      call MPI_RECV(rid, 1, MPI_INTEGER, cpuf, cpuf, MPI_COMM_WORLD, stat, ie)
      a(i+1) = rid
      !print *, "Sent:",a(i)," to Proc:",id+1," Recv:",rid," from Proc:",cpuf
    else if (id == cpuf) then
      ! Edge case, must send to the other end
      call MPI_RECV(rid, 1, MPI_INTEGER, id-1, msgid-1, MPI_COMM_WORLD, stat,ie)
      call MPI_SEND(a(i), 1, MPI_INTEGER, 0, msgid, MPI_COMM_WORLD, ie)
      a(i+1) = rid
      !print *, "Sent:",a(i)," to Proc:",cpu0," Recv:",rid," from Proc:",id-1
    else 
      ! Regular case, sends to the right, receives from the left
      !This can be optimized since the current configuration serializes the code
      !These are blocking sends and receives. Thus processor 15 must wait to
      !recevie from processor 14 which waits for processor 13 and so on. 
      call MPI_RECV(rid, 1, MPI_INTEGER, id-1, msgid-1, MPI_COMM_WORLD, stat,ie)
      call MPI_SEND(a(i), 1, MPI_INTEGER, id+1, msgid, MPI_COMM_WORLD, ie)
      a(i+1) = rid
      !print *, "Sent:",a(i)," to Proc:",id+1," Recv:",rid," from Proc:",id-1
    end if
  end do 

  !This prints outputs from halfway across the ring to show that they have
  !received and sent the right things
  if (id == cpu0) then 
    print *, ''
    print *, 'CPU 0 Outcome:'
    do i = 1, np
      print *, 'CPU 0: ', a(i)
    end do
  end if
  if (id == 6) then 
    print *, ''
    print *, 'CPU 6 Outcome:'
    do i = 1, np
      print *, 'CPU 6: ',a(i)
    end do
  end if

  ! Stoping MPI
  call MPI_FINALIZE(ie)
end program ring




