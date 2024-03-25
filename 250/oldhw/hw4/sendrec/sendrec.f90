

program sendrec
include 'mpif.h'

integer :: id, ie, np
integer :: msgid, src, dest, count
integer :: buffer 
integer :: stat(MPI_STATUS_SIZE)
character(9) :: signal

call MPI_INIT(ie)
call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)
call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)

msgid = 10
src = 0
dest = 1
count = 9

if (id ==  src) then
    signal = 'Send this'
    call MPI_SEND(signal, count, MPI_CHARACTER, dest, msgid, &
       &  MPI_COMM_WORLD, ie)
    print *, 'Processor: ', id, ' sent ', signal
else if (id ==  dest) then 
    call MPI_RECV(signal, count, MPI_CHARACTER, src, msgid, &
       &  MPI_COMM_WORLD, stat, ie)
    print *, 'Processor: ', id, ' received ', signal
end if

call MPI_FINALIZE(ie)

stop
end program sendrec



