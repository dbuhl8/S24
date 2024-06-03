f(x) = m*x + b
g(x) = n*x + c

fit [1:900] f(x) "latency.dat" u 1:2 via m,b
fit [1100:9900] g(x) "latency.dat" u 1:2 via n,c

plot "latency.dat" u 1:2, f(x), g(x)
