


program ring
implicit none

include 'mpif.h'

integer ie, id, np, rid, i
integer msgid, cpu0, cpuf
integer stat(MPI_STATUS_SIZE)
integer, allocatable :: a(:)

call MPI_INIT(ie)
call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)
call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)

cpu0 = 0
cpuf = np-1
msgid = id
allocate(a(np+1))

do i = 1, np
    if(id == i-1) then
        a(1) = i-1
    end if
end do

do i = 1, np
if (id == cpu0) then
    call MPI_SEND(a(i), 1, MPI_INTEGER, id+1, msgid, MPI_COMM_WORLD, ie)
    call MPI_RECV(rid, 1, MPI_INTEGER, cpuf, cpuf, MPI_COMM_WORLD, stat, ie)
    a(i+1) = rid
    !print *, "Sent:",a(i)," to Proc:",id+1," Recv:",rid," from Proc:",cpuf
else if (id == cpuf) then
    call MPI_RECV(rid, 1, MPI_INTEGER, id-1, msgid-1, MPI_COMM_WORLD, stat, ie)
    call MPI_SEND(a(i), 1, MPI_INTEGER, 0, msgid, MPI_COMM_WORLD, ie)
    a(i+1) = rid
    !print *, "Sent:",a(i)," to Proc:",cpu0," Recv:",rid," from Proc:",id-1
else 
    call MPI_RECV(rid, 1, MPI_INTEGER, id-1, msgid-1, MPI_COMM_WORLD, stat, ie)
    call MPI_SEND(a(i), 1, MPI_INTEGER, id+1, msgid, MPI_COMM_WORLD, ie)
    a(i+1) = rid
    !print *, "Sent:",a(i)," to Proc:",id+1," Recv:",rid," from Proc:",id-1
end if
end do 


if (id == cpu0) then 
    print *, ''
    print *, 'CPU 0 Outcome:'
    do i = 1, np
        print *, a(i)
    end do
end if
if (id == 6) then 
    print *, ''
    print *, 'CPU 7 Outcome:'
    do i = 1, np
        print *, a(i)
    end do
end if

call MPI_FINALIZE(ie)
stop

end program ring




