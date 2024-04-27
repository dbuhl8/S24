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
