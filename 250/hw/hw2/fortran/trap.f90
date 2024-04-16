program trap

implicit none

real :: a, b, val
integer :: n, f
print *, "Please enter in the following order separated by commas: interval start, interval end, &
number of data points, and the function number"
! This is going to read a value for a
read (*, *) a, b, n, f

!run specific function
!Need to build a function in fortran.
val = trapInt(f, a, b, n)

print *, "The value of the definite integral is: ", val
write (*, *) val


contains 
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


end program trap

