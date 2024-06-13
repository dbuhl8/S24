
set terminal png size 1600, 900

# Plot 1
reset
set output "spectal.png"
set title 'Spectral Radius v Grid Size'
set xlabel 'Grid Size'
set ylabel 'Spectral Radius'
plot 'dt.dat' u 1:2 with linespoints ps 3


#Plot 2
set output "dt.png"
set title 'Dt_crit v grid_size'
set xlabel 'Grid Size'
set ylabel 'dt_crit'
set log xy
plot "dt.dat" u 1:3 w lp ps 3, x**(-2), x**(-3)

# Plot 3
reset
set output "mpe.png"
set title 'MPE v time'
set xlabel 'Time'
set ylabel 'MPE'
set log y
plot "mpe.dat" u 1:2 title "n=50", "" u 1:3 title "n=100", "" u 1:4 title "n=150", "" u 1:5 title "n=200"


# Plot 4
reset
set output "final_mpe.png"
set title "Final MPE v Grid Size"
set xlabel 'Grid Size'
set ylabel 'E(2)'
set log xy 
set xrange [40:220]
plot "final_mpe.dat" u 1:2 w lp ps 3, x**(-2)

#Plot 5
reset 
set output "r_norm.png"
set title "||R_k|| v k"
set xlabel "k"
set ylabel "||R_k||"
plot "r_norm.dat" u 1:2 w lp ps 3

#Plot 6
reset 
set output "r.png"
set title "Basis Functions"
set xlabel "r"
set ylabel "R_n"
plot "r.dat" u 1:2 title "R1", "" u 1:3 title "R2", "" u 1:4 title "R4", "" u 1:5 title "R6"

