

set terminal png size 1600, 900
# Plot 1
set output "mpe.png"
set title 'Max Pointwise Error between FD and Ch'
set xlabel 'Time'
set ylabel 'Error'
set log y
plot "mpe_error.dat" u 1:2 w l

#Plot 2
reset
set xrange [0.8:1.2]
set output "fdi.png"
set title 'Integral of U over time for FD'
set xlabel 'Time'
set ylabel 'Int of U'
plot "fdint.dat" u 1:2 w l

#Plot 3
reset
set output "mse.png"
set log y
set xrange [20:160]
set title 'Mean Squared Error v N'
set xlabel 'N'
set ylabel 'MSE'
plot "mse_error.dat" u 1:2 w l

