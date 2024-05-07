

set terminal png size 1600,900

# plot 0
set output "1cdt=10-4.png"
set title "AB3 Traj for dt = 10^-4"
plot "AB3Traj.dat" u 1:2 w l title 'X', "" u 1:3 w l title 'Y', "" u 1:4 w l title 'Z'

#plot 1

set output "1cdt<dt_crit.png"
set title "AB3 Traj for dt < dt_crit"
plot "AB3Traj.dat" u 5:6 w l title 'X', "" u 5:7 w l title 'Y', "" u 5:8 w l title 'Z'

#plot 2

set output "1cdt>dt_crit.png"
set title "AB3 Traj for dt > dt_crit"
plot "AB3Traj.dat" u 9:10 w l title 'X', "" u 9:11 w l title 'Y', "" u 9:12 w l title 'Z'

