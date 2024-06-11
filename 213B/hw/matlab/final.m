close all; clear; clc;

u = readmatrix('../fortran/U.dat');
ua = readmatrix('../fortran/UA.dat');

help = size(u)

x = linspace(1,3,202);
t = linspace(0,2,help(1));

[xx,tt] = meshgrid(x,t);

figure(1);
s = surf(xx, tt, u, 'FaceColor','interp');
s.EdgeColor = 'none';
colorbar
%view([0,90])
title('Numerical Solution')
figure(2);
s = surf(xx, tt, ua, 'FaceColor','interp');
s.EdgeColor = 'none';
title('Analytical Solution')
colorbar
%view([0,90])
figure(3);
plot(x,u(1,:),'bx',x,ua(1,:),'r')
xlabel('r')
ylabel('u')
title('Initial Condition')
legend()
