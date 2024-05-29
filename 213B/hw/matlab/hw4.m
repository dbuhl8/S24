close all; clear; clc; 

u = readmatrix('../fortran/u_phys.dat');
x = readmatrix('../fortran/x.dat');
t = readmatrix('../fortran/t.dat');

s = surf(x, t, u, 'FaceColor','interp');
s.EdgeColor = 'none';
