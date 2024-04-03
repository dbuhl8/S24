
module hw2mod

    implicit none

    contains

!   subroutine ones(args)

!       integer :: i, j
!       integer :: n, seed
!       integer, allocatable :: a(:, :), b(:, :)
!       print *, "What size matrices would you like to create? Please type an integer &
!                greater than 3."
!       read (*, *) n
!       print *, "Please choose a seed for the random generator!: Any integer will do:"
!       read (*, *) seed

!       allocate(a(n, n))
!       allocate(b(n, n))
!       !read n, size of matrix, from use,
!       call srand(seed)
!       do j = 1, n, 1
!           do i = 1, n, 1
!               a(i, j) = anint(RAND(0))
!           end do
!       end do

!       call transfMat(A, B, n)
!       ! generate random distrbution of ones in a

!       open (10, file='ones_output.txt', status='unknown')
!       ! loop through a and assign values to b accordingly
!       print *, "The Generated Matrix A: "
!       write (10, *) "The Generated Matrix A: "
!       do i=1, n
!             print *, a(i,:)
!             write (10, *) a(i,:)
!       end do
!       
!       print *, "The Resulting Matrix B: "
!       write (10,*) "The Generated Matrix B: "
!       do i=1, n
!           print *, b(i,:)
!           write (10, *) b(i,:)
!       end do

!   end subroutine ones

    subroutine transfMat(A, B, n)
        ! A, B are both n by n matrices
        implicit none
        integer, intent(in) :: n
        real, intent(in) :: A(:, :)
        real :: B(:, :)
        integer :: i, j

        B = 0.0
    
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
                if (B(i, j) .eq. 3.0) then
                    B(i, j) = 1.0
                else
                    B(i, j) = 0.0
                end if
            end do
        end do
    
    end subroutine transfMat


    function trapInt(f, a, b, n) result(s)

        implicit none

        integer, intent(in) :: f, n
        real, intent(in) :: a, b
        integer :: val
        real :: s, l
        s = 0
        l = (b-a)/n

        if (f == 1) then
            do val = 1, n, 1
                s = s + l*((l*(val-1))**2 + ((l*val)**2))/2
            end do
        else if (f == 2) then
            do val = 1, (n-1), 1
                s = s + l*(sin(a + l*(val-1)) + (sin(a + l*val)))/2
            end do
        end if

    end function trapInt

!   program trap()

!       real :: a, b, val
!       integer :: n, f
!       print *, "Please enter in the following order separated by commas: interval start, interval end, &
!       number of data points, and the function number"
!       ! This is going to read a value for a
!       read (*, *) a, b, n, f

!       !run specific function
!       !Need to build a function in fortran.
!       val = trapInt(f, a, b, n)

!       print *, "The value of the definite integral is: ", val
!       write (*, *) val

!   end program trap


end module hw2mod
