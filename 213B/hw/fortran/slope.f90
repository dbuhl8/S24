


program slope

  implicit none

  integer, parameter :: kr = kind(dble(1.0))
  integer, dimension(4) :: n
  real(kind=kr), dimension(4) :: mse
  real(kind=kr) :: mse_slope
  integer :: i

  open(10, file='mse_error.dat')
  do i = 1, 4
    read(10, *) n(i), mse(i)
  end do 
  close(10)

  mse_slope = -(log(mse(4)/mse(1)))/(log(real(n(4))/n(1)))

  print *, " Computed Slope : ", mse_slope

end program slope
