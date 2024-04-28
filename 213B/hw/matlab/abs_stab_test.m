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

for i = 1:3
    lambda_A=D(i,i);
    xl = real(lambda_A)/abs(lambda_A);
    yl = imag(lambda_A)/abs(lambda_A);
    plot(xl, yl, 'xr');
    line([xl, 0], [yl, 0]);
end

figure(2); 
hold on;

c1 = 23./12.;
c2 = -16./12.;
c3 = 5./12.;
th = 0:0.0001:2*pi;
r = (exp(complex(0, 2*th)) - 4.*exp(complex(0, th)) + 3.)/2.;
%(exp(complex(0, 3*th)) - exp(complex(0, 2*th)))./(c1*exp(complex(0,2*th)) + c2*exp(complex(0, th)) + c3);

x = real(r);
y = imag(r);

plot(x, y);

figure(3); 
hold on;

th = 0:0.0001:2*pi;
r = (exp(complex(0, 3*th)) - 18.*exp(complex(0, 2*th)) ...
    +9.*exp(complex(0, th)) - 2)./(6.*exp(complex(0,3*th)));
%(exp(complex(0, 2*th)) - 4.*exp(complex(0, th)) + 3.)/2.;
%(exp(complex(0, 3*th)) - exp(complex(0, 2*th)))./(c1*exp(complex(0,2*th)) + c2*exp(complex(0, th)) + c3);

x = real(r);
y = imag(r);

plot(x, y);

figure(4); 
hold on;

%th = 0:0.0001:2*pi;
b = [2./9, 1./3, 4./9];
x = -3:0.1:1;
y = -3:0.1:3;
[xx, yy] = meshgrid(x, y);
z = complex(xx, yy); %exp(complex(0., th));
%z2 = z1.^2; %exp(complex(0., 2.*th));
%z3 = z1.^3; %exp(complex(0., 3.*th));
%for i = 1:sz(1)
    %for j = 1:sz(2)
        %z = z1(i,j); 
        %Smat = [1 + b(1)*z, b(2)*z, b(3)*z;...
                %z*(b(1) - 1./2), 1 + b(2)*z, b(3)*z;...
                %b(1)*z, z*(b(2) - 3./4), 1 + b(3)*z];
        %s(i, j) = det(Smat);
        %s = sum(b)*z + (z.^2).*((b(2)^2)/2. + (3*b(3)^2)/4.) ...
                %+ (z.^3)*(prod(b) + (b(3)*(b(2)-3./4)*(b(1)-1./2)) ...
                %+ (b(3)*b(2)*(b(1)-1./2)) + (b(3)*b(1)*(b(2) - 3./4)));
    %end 
%end 
%(exp(complex(0, 3*th)) - 18.*exp(complex(0, 2*th)) ...
%    +9.*exp(complex(0, th)) - 2)./(6.*exp(complex(0,3*th)));
%(exp(complex(0, 2*th)) - 4.*exp(complex(0, th)) + 3.)/2.;
%(exp(complex(0, 3*th)) - exp(complex(0, 2*th)))./(c1*exp(complex(0,2*th)) + c2*exp(complex(0, th)) + c3);

%x = real(cos(th));
%y = imag(sin(th));

bsum = sum(b)
bprod = prod(b)

k = sin(xx).*yy

s = (1 + b(1)*z).*(1 + b(2)*z).*(1 + b(3)*z) + (z.^3)*prod(b) ...
    + (z.^3)*(b(3)*(b(2) - 3./4)*(b(1) - 1./2))...
    - (b(3)*b(1)*z.^2).*(1+b(2)*z) - (b(2)*(b(1) - 1./2)*z.^2).*(1 + b(3)*z)...
    - (b(3)*(b(2) - 3./4)*z.^2).*(1 + b(1)*z);
    

%smax = max(s, [], "all") 

%[row, column] = find(abs(s) >= 1.0);
%s(row, column) = 1.0;

%smax = max(s, [], "all") 

%contour(xx(row, column), yy(row, column), transpose(abs(s(row, column))));
%surf(xx, yy, abs(s));
[C1,h1]=contour(xx, yy, k);
clabel(C1,h1)
