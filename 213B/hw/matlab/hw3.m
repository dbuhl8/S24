clear;
close all;
clc;


U = readmatrix('../fortran/u.dat');

x = linspace(0, 2, 81);
y = flip(linspace(0, 1, 51));

[X, Y] = meshgrid(x, y);

figure(1);

surf(X, Y, U)
xlabel("X")
ylabel("Y")
zlabel("U")
title("Q1: Numerical Solution")

figure(2)

I2 = find(Y == 0.2);
I5 = find(Y == 0.5);

hold on;
plot(X(I2), U(I2), 'r');
plot(X(I5), U(I5), 'b');
xlabel("X")
ylabel("U")
legend('Y=0.2', 'Y=0.5')
title("Q1: Plot of U(x)")


figure(3)

U2 = readmatrix('../fortran/u2.dat');

x = linspace(-1, 1, 100);
t = linspace(0, 10, 10001);
[X,T] = meshgrid(x, t);
s = surf(X, T, U2, 'FaceColor','interp');
s.EdgeColor = 'none';
xlabel('X')
ylabel('T')
zlabel('U')
title('Evolution of Heat Equation')
