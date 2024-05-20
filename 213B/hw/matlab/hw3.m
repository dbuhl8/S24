clear;
close all;
clc;


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

U2 = readmatrix('../fortran/u2.dat');
X2 = readmatrix('../fortran/x2.dat');
T2 = readmatrix('../fortran/t2.dat');

sa = surf(X2, T2, U2, 'FaceColor','interp');
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

UU = readmatrix('../fortran/us.dat');
XX = readmatrix('../fortran/xs.dat');
TT = readmatrix('../fortran/ts.dat');

ss = surf(XX, TT, UU, 'FaceColor','interp');
ss.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Spectral Evolution of Heat Equation')

