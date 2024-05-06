
clear; close all; clc;
format long

figure(1);
b = [1./6, 2./3, 1./6];
%x = -5:0.01:5;
%y = -5:0.01:5;
x=-6:0.01:0.5;
y=-5:0.1:5;
[xx, yy] = meshgrid(x, y);
z = complex(xx, yy);

s = ((1 + b(1)*z).*(1 + z*(b(2) - 0.25)).*(1 + b(3)*z) + (z.^3)*prod(b) ...
    + (z.^3)*(b(3)*(b(2) - 1.)*(b(1) - 1./4))...
    - (b(3)*b(1)*z.^2).*(1+z*(b(2)-0.25)) - (b(2)*(b(1) - 1./4)*z.^2).*(1 + b(3)*z)...
    - (b(3)*(b(2) - 1.)*z.^2).*(1 + b(1)*z))./(1 - 0.25*z);
 
[C1,h1]=contour(xx, yy, abs(s), [1,1]);
[rw, cl] = find(abs(s) < 1.0);
min_x = min(xx(rw, cl), [], "all")
%clabel(C1, h1);
xlabel("Re(\lambda)")
ylabel("Im(\lambda)")
title('Region of Abs Stability for Implicit RK3')
hold on;

B = [-1.0, 3.0, -5.0, 7.0; 0.0, -2.0, 4.0, -6.0; 0.0, 0.0, -4.0, 6.0; 0.0, 0.0, 0.0, -16.0];
[V, D] = eig(B);
for i = 4:4
    lambda_A=D(i,i);
    xl = real(lambda_A);%/abs(lambda_A);
    yl = imag(lambda_A);%/abs(lambda_A);
    %plot(xl, yl, 'xr');
    %line([xl, 0], [yl, 0], 'Color', 'red', 'LineStyle', '--');
end

figure(2);


i = complex(0, 1);
th = 0:0.001:2*pi;
r = 12*(exp(i*3*th) - (exp(2*i*th) + exp(i*th) + 1)/3)./(23*exp(2*i*th) - 2*exp(i*th) + 3);
x = real(r);
y = imag(r);

plot(x,y);
xlabel("Re(\lambda)")
ylabel("Im(\lambda)")
title("Region of Absolute Stability")


