clc; close; clear;
A = [0.0, 10.0, -10.0; -100.0, -1.0, 0.0; 0.0, 10.0, -100.0];

[V, D] = eig(A);
for dt = 0.005:0.00000001:0.01
    format long g
    Current_DT = dt
%    for i = 1:3
        lambda_A = min(real(D),[], "all");
        Current_Lambda = lambda_A;
        C = [0.0, 1.0, 0.0; 0.0, 0.0, 1.0; (5.*dt*lambda_A)/12., (-16.0*dt*lambda_A)/12., 1 + (23.*dt*lambda_A)/12.];
        [C1, C2] = eig(C);
        for j = 1:3
            lambda_C(j) = abs(C2(j, j));
        end 
        if (max(lambda_C) <= 1.)
            Current_DT_Lambda = true;
        else 
            Current_DT_Lambda = false;
            break
        end
        %For_Current_Lambda_DT = lambda_C
%    end 
end

figure(1);
hold on;

c1 = 23./12.;
c2 = -16./12.;
c3 = 5./12.;
th = 0:0.0001:2*pi;
r = (exp(complex(0, 3*th)) - exp(complex(0, 2*th)))./(c1*exp(complex(0,2*th)) + c2*exp(complex(0, th)) + c3);

x = real(r);
y = imag(r);

plot(x, y);
title('Region of Abs Stability for AB3')

for i = 1:3
    lambda_A=D(i,i);
    xl = real(lambda_A)/abs(lambda_A);
    yl = imag(lambda_A)/abs(lambda_A);
    plot(xl, yl, 'xr');
    line([xl, 0], [yl, 0], 'Color', 'red', 'LineStyle', '--');
end

figure(2); 
hold on;

r = (exp(complex(0, 2*th)) - 4.*exp(complex(0, th)) + 3.)/2.;

x = real(r);
y = imag(r);

plot(x, y);
title('Region of Abs Stability for LMM')
line([0 0], [-3 3], 'Color', 'black', 'LineStyle', '--');  %x-axis
line([-1, 5], [0 0], 'Color', 'black', 'LineStyle', '--');  %y-axis

figure(3); 
hold on;

r = (11.*exp(complex(0, 3*th)) - 18.*exp(complex(0, 2*th)) ...
    +9.*exp(complex(0, th)) - 2)./(6.*exp(complex(0,3*th)));

x = real(r);
y = imag(r);

plot(x, y);
title('Region of Abs Stability for BDF3')
line([0 0], [-4.5 4.5], 'Color', 'black', 'LineStyle', '--');  %x-axis
line([-1 7], [0 0], 'Color', 'black', 'LineStyle', '--');  %y-axis

figure(4); 
hold on;

b = [2./9, 1./3, 4./9];
x = -3:0.1:1;
y = -3:0.1:3;
[xx, yy] = meshgrid(x, y);
z = complex(xx, yy);

s = (1 + b(1)*z).*(1 + b(2)*z).*(1 + b(3)*z) + (z.^3)*prod(b) ...
    + (z.^3)*(b(3)*(b(2) - 3./4)*(b(1) - 1./2))...
    - (b(3)*b(1)*z.^2).*(1+b(2)*z) - (b(2)*(b(1) - 1./2)*z.^2).*(1 + b(3)*z)...
    - (b(3)*(b(2) - 3./4)*z.^2).*(1 + b(1)*z);

[C1,h1]=contour(xx, yy, abs(s), [1,1]);
%clabel(C1,h1);
title('Region of Abs Stability for Explicit RK3')
