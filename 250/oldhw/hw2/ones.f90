program ones

implicit none


integer :: i, j
integer :: n, seed
integer, allocatable :: a(:, :), b(:, :)
print *, "What size matrices would you like to create? Please type an integer &
 greater than 3."
read (*, *) n
print *, "Please choose a seed for the random generator!: Any integer will do:"
read (*, *) seed

allocate(a(n, n))
allocate(b(n, n))
!read n, size of matrix, from use,
call srand(seed)
do j = 1, n, 1
    do i = 1, n, 1
        a(i, j) = anint(RAND(0))
    end do 
end do

b = transfMat(a, n)
! generate random distrbution of ones in a

open (10, file='ones_output.txt', status='unknown')
! loop through a and assign values to b accordingly
print *, "The Generated Matrix A: "
write (10, *) "The Generated Matrix A: "
do i=1, n
      print *, a(i,:)
      write (10, *) a(i,:)
end do

print *, "The Resulting Matrix B: "
write (10,*) "The Generated Matrix B: "
do i=1, n
      print *, b(i,:)
      write (10, *) b(i,:)
end do




contains
function transfMat(A, n) result(B)
    implicit none
    integer, intent(in) :: n
    integer :: s=0
    integer, intent(in), dimension(n, n) :: A
    integer, dimension(n, n) :: B
    integer :: i, j    
    
    ! Instantializes B    
    do j = 1, n, 1
        do i = 1, n, 1
            B(i, j) = 0
        end do
    end do

    !Down, completed
    B(2:n, :) = B(2:n, :) +  A(1:(n-1), :)
    !Up, completed
    B(1:(n-1), :) = B(1:(n-1), :) + A(2:n, :) 
    !Left, completed
    B(:, 1:(n-1)) = B(:, 1:(n-1)) + A(:, 2:n) 
    !Right, completed
    B(:, 2:n) = B(:, 2:n) + A(:, 1:(n-1)) 
    !Up-right, completed
    B(1:(n-1), 2:n) = B(1:(n-1), 2:n) + A(2:n, 1:(n-1)) 
    !Up-left, completed
    B(1:(n-1), 1:(n-1)) = B(1:(n-1), 1:(n-1)) + A(2:n, 2:n) 
    !Down right, completed
    B(2:n, 2:n) = B(2:n, 2:n) + A(1:(n-1), 1:(n-1)) 
    !Down left, completed
    B(2:n, 1:(n-1)) = B(2:n, 1:(n-1)) + A(1:(n-1), 2:n) 
    
    do j = 1, n, 1
        do i = 1, n, 1
            if (B(i, j) == 3) then
                B(i, j) = 1
            else
                B(i, j) = 0
            end if
        end do
    end do
       
end function transfMat

end program ones

