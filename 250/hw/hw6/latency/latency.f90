!File: Latency.f90
!Author: Dante Buhl
!Dependencies: MPI

program latency

  implicit none
  include 'mpif.h'

  integer :: i, nk, j, fn=10
  integer :: id, src, dest, msgid, ie, np
  integer :: stat(MPI_STATUS_SIZE)
  integer, parameter :: kr=kind(dble(1.0)), num_send=10000, num_repeat = 100
  real(kind=kr), dimension(num_send) :: data, data2, timedata, tot_time_data
  real(kind=kr) :: start, finish

  call MPI_INIT(ie)
  call MPI_COMM_RANK(MPI_COMM_WORLD, id, ie)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ie)

  i = 0
  src = 0
  dest = 1
  msgid = 1
  timedata = 0. 

  if(id.eq.0) then
    do i = 1, num_send
      data(i) = exp(real(i, kind=kr)/num_send)
    end do 
    data2 = 0.
  else
    do i = 1, num_send
      data2(i) = exp(real(i, kind=kr)/num_send)
    end do 
    data = 0.
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ie)

  ! Sending Messages Between Processors
  if (id.eq.src) then
    ! If processor 0, start by sending and then receiving
    do i = 1, num_send
      ! in this loop the number of data being sent changes
      start = MPI_WTIME()
      do j = 1, num_repeat
        ! in this loop data is sent back and forth 50 times.
        call MPI_SEND(data(1:i),i,MPI_REAL8,1,i*1000+2*j-1,MPI_COMM_WORLD,ie)
        call MPI_RECV(data2(1:i),i,MPI_REAL8,1,i*1000+2*j,MPI_COMM_WORLD,stat,ie)
      end do 
      finish = MPI_WTIME()
      timedata(i) = timedata(i) + finish-start
    end do 
  else
    ! If processor 1, start by receiving and then sending
    do i = 1, num_send
      ! in this loop the number of data being sent changes
      start = MPI_WTIME()
      do j = 1, num_repeat
        ! in this loop data is sent back and forth 50 times.
        call MPI_RECV(data(1:i),i,MPI_REAL8,0,i*1000+2*j-1,MPI_COMM_WORLD,stat,ie)
        call MPI_SEND(data2(1:i),i,MPI_REAL8,0,i*1000+2*j,MPI_COMM_WORLD,ie)
      end do 
      finish = MPI_WTIME()
      timedata(i) = timedata(i) + finish-start
    end do 
  end if
  
  call MPI_REDUCE(timedata, tot_time_data, num_send, MPI_REAL8, MPI_SUM, 0,&
    MPI_COMM_WORLD, ie)

  ! Write data to outfile to produce a plot
  if (id.eq.0) then
    tot_time_data = tot_time_data/(num_repeat*4)
    open(fn, file='latency.dat')
      do i = 1, num_send
        write(fn, '(I6, F24.15)') i, tot_time_data(i)
      end do 
    close(fn)
  end if

  call MPI_FINAlIZE(ie)
end program latency
