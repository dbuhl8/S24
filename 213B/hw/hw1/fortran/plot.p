

set terminal png size 1600,900

# plot 0
set output "question2a.png"
plot "df.dat" u 1:2

# plot 1 
set output "BDF3_plot.png"
plot "df.dat" u 1:2, "" u 3:4, "" u 5:6

# plot 3
set output "rk3_am3_v_actual.png"
plot "rk3.dat" u 2:3, "am3.dat" u 2:3, "actual.dat" u 2:3

# plot 4
set log y
set output "error_v_t.png"
plot "error1.dat" u 1:2, "" u 1:3, "error2.dat" u 1:2, "" u 1:3, "error3.dat" u 1:2, "" u 1:3, "error4.dat" u 1:2, "" u 1:3

# plot 2
set log y
set log x
set output "BDF3_error.png"
plot "e.dat" u 1:2, "" u 1:3

# plot 5
set log y
set log x
set output "final_error.png"
plot "finalerror.dat" u 1:2, "" u 1:3 w lp, "" u 1:4, "" u 1:5 w lp



