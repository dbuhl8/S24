set terminal png size 1600,900

# plot 0
set output "rk3_fail.png"
set title "RK3 Traj for dt = dt_crit"
plot "RK3.dat" u 1:2 w l title 'X', "" u 1:3 w l title 'Y', "" u 1:4 w l title 'Z', "" u 1:5 w l title 'W'

# plot 1
set output "rk3_pass.png"
set title "RK3 Traj for dt = dt_crit+dt_pert"
plot "RK3.dat" u 6:7 w l title 'X', "" u 6:8 w l title 'Y', "" u 6:9 w l title 'Z', "" u 6:10 w l title 'W'
