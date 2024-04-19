

program driver

    use LinAl
    use hw2mod

    implicit none

    !variables
    integer :: i, j, n, seed, f
    real, allocatable :: A(:, :), B(:, :)
    real :: s, end, start


    print *, " "
    print *, "-------------------------------------------------------"
    print *, " "
    print *, " Question 1: Trapezoidal Rule"
    print *, " "

    print *, "Please enter in the following order separated by commas: interval start, interval end, &
    number of data points, and the function number (2 = sin, 1 = x^2)"
    ! This is going to read a value for a
    read (*, *) start, end, n, f

    s = trapint(f, start, end, n)

    print "(A, F8.3, A, I5, A)", " The integral is approximated to ", s, " after ", n, " iterations."


    !repeat to test convergence
    print *, " "
    print *, " Convergence Trials for sin(x)"
    print *, " "
    print *, "Actual Value: 2 "

    do i = 1, 25
        print "(A, I5, A, F8.3)", "Iterations: ", i*10, ", Computed Value: ", trapint(2, 0., pi, i*10)
    end do

    print *, " "
    print *, " Convergence Trials for x^2"
    print *, " "

    print *, "Actual Value: 8/3 ~ 2.6666"

    do i = 1, 25
        print "(A, I5, A, F8.3)", "Iterations: ", i*10, ", Computed Value: ", trapint(1, 0., 2., i*10)
    end do

    print *, " "
    print *, "-------------------------------------------------------"
    print *, " "
    print *, " Question 2: Ones"
    print *, " "

    !generate input
    print *, "What size matrices would you like to create? Please type an integer &
             greater than 3."
    read (*, *) n

    allocate(A(n, n), B(n, n))
    A = 0.0
    B = 0.0

    print *, "Please choose a seed for the random generator!: Any integer will do:"
    read (*, *) seed
    call srand(seed)

    do j = 1, n
        do i = 1, n
            A(i, j) = anint(rand(0))
        end do
    end do
    

    !call subroutines
    print *, " "
    print *, "Original matrix obtained for A"
    print *, " "

    call printmat(A, n, n)
    
    call transfMat(A, n)

    print *, " "
    print *, "Resulting matrix obtained, B"
    print *, " "

    call printmat(A, n, n)

    
    print *, " "
    print *, "-------------------------------------------------------"
    print *, " "

    deallocate(A, B)


end program driver


