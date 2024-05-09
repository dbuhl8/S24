


program ping
implicit none
include 'mpif.h'

integer i
integer id, src, dest, msgid, ie, np
integer stat(MPI_STATUS_SIZE)

call MPI_INIT(ie)
call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)
call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)

i = 0
src = 0
dest = 1
msgid = 1

do while (i /= 7)
    if (mod(i,2) == 0) then
        if (id == 0) then 
            i = i + 1
            call MPI_SEND(i, 1, MPI_INTEGER, dest, msgid, MPI_COMM_WORLD, ie)
            print *, "Hey I'm processor ",id," and I've just sent ",i," TO   processor ",dest," on msg :", msgid
            msgid = msgid + 1
        else 
            call MPI_RECV(i, 1, MPI_INTEGER, src, msgid, MPI_COMM_WORLD, stat, ie)
            print *, "Hey I'm processor ",id," and I've just recv ",i," FROM processor ",src," on msg :", msgid
            msgid = msgid + 1
        end if
    else 
        if (id == 1) then 
            i = i + 1
            call MPI_SEND(i, 1, MPI_INTEGER, src, msgid, MPI_COMM_WORLD, ie)
            print *, "Hey I'm processor ",id," and I've just sent ",i," to   processor ",src," on msg :", msgid
            msgid = msgid + 1
        else 
            call MPI_RECV(i, 1, MPI_INTEGER, dest, msgid, MPI_COMM_WORLD, stat, ie)
            print *, "Hey I'm processor ",id," and I've just recv ",i," from processor ",dest," on msg :", msgid
            msgid = msgid + 1
        end if
    end if
end do

call MPI_FINAlIZE(ie)

stop

end program ping
