function plot_LMM_stability_region()

alpha_AB2= [1 -1/3 -1/3 -1/3 ];
beta_AB2 = [0 23/12 -1/6 -3/12];

roots(alpha_AB2)

z_AB2 = compute_root(alpha_AB2,beta_AB2);



figure(1)
clf
plot(real(z_AB2),imag(z_AB2),'b','Linewidth',1.5)
hold 
grid
xlabel('{Re}$(z)$','Interpreter','Latex')
ylabel('{Im}$(z)$','Interpreter','Latex')
set(gca,'Fontsize',22)
axis equal
axis([-1.5 0.5 -1 1])

end



function z = compute_root(alpha,beta)

theta=linspace(0,2*pi,500); 

rho   = polyval(alpha,exp(1i*theta));
sigma = polyval(beta,exp(1i.*theta)); 

z=rho./sigma;

end
