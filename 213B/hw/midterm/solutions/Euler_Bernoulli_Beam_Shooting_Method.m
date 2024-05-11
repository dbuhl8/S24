close all;
clear all;

EI = 1;                                 % flexural rigidity
L = 1;                                  % length of beam

C = [0 1 0 0 0 0 0 0;                   % Jacobian evolution matrix
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0];
 
j0 =[0; 0; 1; 0; 0; 0; 0; 1];           % Jacobian initial condition
    
j = expm(C*L)*j0;                       % Jacobian at x = L

J = [j(1) j(5);                         % Jacobian components needed for Newton's iteration
    j(2) j(6)];

NSTEPS = 6e4;                           

DT = L/NSTEPS;

IOSTEPS = 500;

v0 = rand(2,1);                         % initial velocity

y0 = [0;0;v0];                          % initial condition for ODE BVP

y_analytical = @(x) (x.^2/120-x.^3/90+x.^6/360)/EI % analytical solution with forcing q(x) = x^2

e = norm(y0);

tol = 1e-3;

f  = @(x,t) [x(2);x(3);x(4);(t^2)/EI];

nIter = 2;

for i = 1:nIter
    
[y,t] = RK4_Method(f,y0,DT,NSTEPS,IOSTEPS);

tt{:,i} = t;
yy{:,i} = y;
  
v1 = v0-J\[y(1,end);y(2,end)];
    
y0 = [0;0;v1];
  
v0 = v1;
end


figure;



plot(tt{:,nIter}(1,:),yy{:,nIter}(1,:),'.r')
hold on
plot(tt{:,nIter}(1,:),y_analytical(tt{:,nIter}(1,:)),'b')

legend('numerical','analytical')


title('Shooting Method for Euler-Bernoulli Beam')
xlabel('$x$','interpreter','latex','fontsize',15)
ylabel('$y(x)$','interpreter','latex','fontsize',15)

print(gcf,'Euler_Bernoulli_Comparison.png','-dpng','-r300'); 
figure;

err = abs(y_analytical(tt{:,nIter}(1,:))-yy{:,nIter}(1,:));
semilogy(tt{:,nIter}(1,:),err);

xlabel('$x_k$','interpreter','latex','fontsize',15)
ylabel('$\log |y(x_k)-u(x_k)| $','interpreter','latex','fontsize',15)
print(gcf,'Euler_Bernoulli_Error.png','-dpng','-r300'); 


