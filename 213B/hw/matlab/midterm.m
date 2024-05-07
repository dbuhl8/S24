%Set up
clear; close all; clc;
format long

% Find dt
figure(1);
b = [1./6, 2./3, 1./6];
%x = -5:0.01:5;
%y = -5:0.01:5;
x=-5.425:0.000001:-5.415;
y=-0.2:0.1:0.2;
[xx, yy] = meshgrid(x, y);
z = complex(xx, yy);

s = ((1 + b(1)*z).*(1 + z*(b(2) - 0.25)).*(1 + b(3)*z) + (z.^3)*prod(b) ...
    + (z.^3)*(b(3)*(b(2) - 1.)*(b(1) - 1./4))...
    - (b(3)*b(1)*z.^2).*(1+z*(b(2)-0.25))...
    - (b(2)*(b(1) - 0.25)*z.^2).*(1 + b(3)*z)...
    - (b(3)*(b(2) - 1.)*z.^2).*(1 + b(1)*z))./(1 - 0.25*z);

[rw, cl] = find(abs(s) < 1.0);
min_x = (min(xx(rw, cl), [], "all"))/(-16.0)
clear;


% Plotting region of absolute stability
b = [1./6, 2./3, 1./6];
x=-5.5:0.01:0.5;
y=-5:0.1:5;
[xx, yy] = meshgrid(x, y);
z = complex(xx, yy);

s = ((1 + b(1)*z).*(1 + z*(b(2) - 0.25)).*(1 + b(3)*z) + (z.^3)*prod(b) ...
    + (z.^3)*(b(3)*(b(2) - 1.)*(b(1) - 1./4))...
    - (b(3)*b(1)*z.^2).*(1+z*(b(2)-0.25)) - (b(2)*(b(1) - 1./4)*z.^2).*(1 + b(3)*z)...
    - (b(3)*(b(2) - 1.)*z.^2).*(1 + b(1)*z))./(1 - 0.25*z);

[C1,h1]=contour(xx, yy, abs(s), [1,1]);
%clabel(C1, h1);
xlabel("Re(\lambda)")
ylabel("Im(\lambda)")
title('Region of Abs Stability for Implicit RK3')
hold on;

%Double checking eigenvalues are correct
B = [-1.0, 3.0, -5.0, 7.0; 0.0, -2.0, 4.0, -6.0; 0.0, 0.0, -4.0, 6.0; 0.0, 0.0, 0.0, -16.0];
[V, D] = eig(B);
for i = 1:4
    lambda_A=D(i,i);
    xl = real(lambda_A);%/abs(lambda_A);
    yl = imag(lambda_A);%/abs(lambda_A);
    plot(xl, yl, 'xr');
    line([xl, 0], [yl, 0], 'Color', 'red', 'LineStyle', '--');
end

%Problem 2
figure(2);

i = complex(0, 1);
th = 0:0.001:2*pi;
r = 12*(exp(i*3*th) - (exp(2*i*th) + exp(i*th) + 1.)/3.)./(23*exp(2*i*th) - 2*exp(i*th) + 3);
x = real(r);
y = imag(r);

plot(x,y);
xlabel("Re(\lambda)")
ylabel("Im(\lambda)")
title("Region of Absolute Stability")


figure(3);

%shooting method
np4 = 60000;
dt = 1./np4;

v1(1) = 0.;
v2(1) = 0.;
v1(2) = 1.;
v2(2) = 1.;

y0 = [0.;0.;v1(1);v2(1)];
y1 = [0.;0.;v1(2);v2(2)];

[t0, yout0] = ode45(@diffeq, 0.:dt:1., y0);
[t1, yout1] = ode45(@diffeq, 0.:dt:1., y1);

yact = (t0.^6)/360. - (t0.^3)/90. + (t0.^2)/120.;
dyact = (t0.^5)/60. - (t0.^2)/30. + (t0)/60.;

y_error(1) = yact(np4+1) - yout0(np4+1, 1);
y_error(2) = yact(np4+1) -   yout1(np4+1, 1);
dy_error(1) = dyact(np4+1) - yout0(np4+1, 2);
dy_error(2) = dyact(np4+1) - yout1(np4+1, 2);

error_tot = sqrt(dy_error(2)^2+y_error(2)^2);
tol = 1e-10;

  dv1 = v1(2) - v1(1);
  dv2 = v2(2) - v2(1);
  de1 = y_error(2) - y_error(1);
  de2 = dy_error(2) - dy_error(1);
  v1_up = v1(2) - y_error(2)*(dv1/de1)
  v2_up = v2(2) - dy_error(2)*(dv2/de2)

while error_tot > tol
  dv1 = v1(2) - v1(1);
  dv2 = v2(2) - v2(1);
  de1 = y_error(2) - y_error(1);
  de2 = dy_error(2) - dy_error(1);
  v1(3) = v1(2) - y_error(2)*(dv1/de1);
  v2(3) = v2(2) - dy_error(2)*(dv2/de2);

  y0 = [0.;0.;v2(3);v1(3)];
  v1(1:2) = v1(2:3);
  v2(1:2) = v2(2:3);

  [t0, yout0] = ode45(@diffeq, 0.:dt:1., y0);
  y_error(1) = y_error(2);
  y_error(2) = yact(np4+1) - yout0(np4+1, 1);
  dy_error(1) = dy_error(2);
  dy_error(2) = dyact(np4+1) - yout0(np4+1, 2);

  %! update error
  error_tot = sqrt(dy_error(2)^2+y_error(2)^2);
end 

v1_vec = v1;
v2_vec = v2;

plot(t0, yout0(:, 1))
hold on;
plot(t0, yact)

figure(4)

plot(t0, yout0(:, 2))
hold on;
plot(t0, dyact)

clear;
figure(5)

y0 = [1.;1.;1.;1.];
[t0, yout0] = ode45(@diffeq2, 0.:0.001:100., y0);

plot(t0, yout0(:,1))
hold on;
plot(t0, yout0(:,2))
plot(t0, yout0(:,3))
plot(t0, yout0(:,4))

clear; 
%solve a linear system really quickly
