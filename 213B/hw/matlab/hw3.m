clear;
close all;
clc;


U = readmatrix('../fortran/u.dat');
X = readmatrix('../fortran/x.dat');
Y = readmatrix('../fortran/y.dat');
F = readmatrix('../fortran/f.dat');
G = readmatrix('../fortran/g.dat');

%x = linspace(0, 2, 5);
%y = flip(linspace(0, 1, 5));

%[X, Y] = meshgrid(x, y);

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

%x = linspace(-1, 1, 100);
%t = linspace(0, 2, 100);
%[X,T] = meshgrid(x, t);
s = surf(X2, T2, U2, 'FaceColor','interp');
s.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Analytical Evolution of Heat Equation')

figure(6)

UF = readmatrix('../fortran/uf.dat');
XF = readmatrix('../fortran/xf.dat');
TF = readmatrix('../fortran/tf.dat');

s = surf(XF, TF, UF, 'FaceColor','interp');
s.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Finite Difference Evolution of Heat Equation')

figure(7)

[xt,D1t,D2t] = get_GCL_points_and_D_matrices(-1, 1, 10);

%help = D2t

US = readmatrix('../fortran/us.dat');
XS = readmatrix('../fortran/xs.dat');
TS = readmatrix('../fortran/ts.dat');

s = surf(XS, TS, US, 'FaceColor','interp');
s.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Spectral Evolution of Heat Equation')

