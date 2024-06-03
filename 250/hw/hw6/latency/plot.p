f(x) = m*x + b
g(x) = n*x + c

fit [1:450] f(x) "latency.dat" u 1:2 via m,b
fit [550:5000] g(x) "latency.dat" u 1:2 via n,c

plot "latency.dat" u 1:2, f(x), g(x)
