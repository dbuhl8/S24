
program stuff

  implicit none
  integer :: m,n,i
  integer, allocatable :: A(:,:)

  
  open(10, file='IC')
  read(10, "(2I6)") m,n 
  allocate(A(m,n))
  do i = 1, n
    read(10, '('//trim(str(m))//'I1)') A(:,i)
  end do 
  close(10)

  open(11, file='ICp')
  do i = 1, m
    write(11, '('//trim(str(n))//'I1)') A(i,:)
  end do 
  close(11)

  contains 

    character(len=20) function str(k)
      ! "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
    end function str

end program stuff
