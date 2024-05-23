clear; close all; clc;

figure(1)

NU1 = readmatrix('../fortran/nu1.dat');
NX1 = readmatrix('../fortran/nx1.dat');
NT1 = readmatrix('../fortran/nt1.dat');


ss = surf(NX1, NT1, NU1, 'FaceColor','interp');
ss.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Finite Difference Evolution of KS IBVP')
 
figure(2)
[row,col] = find(abs(NT1-62) < 0.001);

plot(NX1(row,col), NU1(row,col))
xlabel('X')
ylabel('U')
title('KS IBVP @ T = 62')

