set terminal png size 1600,900

# plot 0
set output "rk3_fail.png"
set title "RK3 Traj for dt = dt_{crit}"
plot "RK3.dat" u 1:2 w l title 'X', "" u 1:3 w l title 'Y', "" u 1:4 w l title 'Z', "" u 1:5 w l title 'W'

# plot 1
set output "rk3_pass.png"
set title "RK3 Traj for dt = dt_{crit}+dt_{pert}"
plot "RK3.dat" u 6:7 w l title 'X', "" u 6:8 w l title 'Y', "" u 6:9 w l title 'Z', "" u 6:10 w l title 'W'

# plot 5
set output "rk3_act.png"
set title "AB3 Traj for dt = 0.001"
plot "prob2_act.dat" u 1:2 w l title 'X', "" u 1:3 w l title 'Y', "" u 1:4 w l title 'Z', "" u 1:5 w l title 'W'

# plot 2
set output "rk4_dy.png"
set title "Shooting Method for DY"
plot "rk4.dat" u 1:3 w l title 'Analytical', "" u 1:5 w l title 'Numerical'

# plot 3
set output "rk4_y.png"
set title "Shooting Method for Y"
plot "rk4.dat" u 1:2 w l title 'Actual', "" u 1:4 w l title 'Numerical'

# plot 4
set output "rk4_error.png"
set title "Shooting Method Error"
plot "rk4.dat" u 1:6 w l title 'Error'

