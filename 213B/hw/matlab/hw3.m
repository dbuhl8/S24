clear;
close all;
clc;

format long

U = readmatrix('../fortran/u.dat');
X = readmatrix('../fortran/x.dat');
Y = readmatrix('../fortran/y.dat');
F = readmatrix('../fortran/f.dat');
G = readmatrix('../fortran/g.dat');

figure(1);

surf(X, Y, U)
xlabel("X")
ylabel("Y")
zlabel("U")
title("Q1: Numerical Solution")

figure(2);

surf(X, Y, F)
xlabel("X")
ylabel("Y")
zlabel("U")
title("Q1: Forcing")

figure(3);

surf(X, Y, G)
xlabel("X")
ylabel("Y")
zlabel("U")
title("Q1: Boundary Condition")


figure(4)

I2 = find(Y == 0.2);
I5 = find(Y == 0.5);

hold on;
plot(X(I2), U(I2), 'r');
plot(X(I5), U(I5), 'b');
xlabel("X")
ylabel("U")
legend('Y=0.2', 'Y=0.5')
title("Q1: Plot of U(x)")


figure(5)

UA = readmatrix('../fortran/u2.dat');
XA = readmatrix('../fortran/x2.dat');
TA = readmatrix('../fortran/t2.dat');

sa = surf(XA, TA, UA, 'FaceColor','interp');
sa.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Analytical Evolution of Heat Equation')

figure(6)

UF = readmatrix('../fortran/uf.dat');
XF = readmatrix('../fortran/xf.dat');
TF = readmatrix('../fortran/tf.dat');

sf = surf(XF, TF, UF, 'FaceColor','interp');
sf.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Finite Difference Evolution of Heat Equation')

figure(7)

US = readmatrix('../fortran/us.dat');
XS = readmatrix('../fortran/xs.dat');
TS = readmatrix('../fortran/ts.dat');

ss = surf(XS, TS, US, 'FaceColor','interp');
ss.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Spectral Evolution of Heat Equation')

figure(8)

E = readmatrix('../fortran/q2error.dat');
N = [5; 10; 15; 20; 25; 30; 50; 100; 150; 200];

semilogy(N, E(:,1), 'x', N, E(:,2), 'x')
xlabel('N')
ylabel('Error')
title('Error as a function of grid size')
legend('Spectral', 'Finite Difference')

figure(9) 

NU1 = readmatrix('../fortran/nu1.dat');
NX1 = readmatrix('../fortran/nx1.dat');
NT1 = readmatrix('../fortran/nt1.dat');

%[nur, nuc] = find(NT1 < 0.04);

%s = surf(NX1(nur,nuc), NT1(nur,nuc), NU1(nur,nuc), 'FaceColor','interp');
s = surf(NX1, NT1, NU1, 'FaceColor','interp');
ss.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Finite Difference Evolution of KS IBVP')


