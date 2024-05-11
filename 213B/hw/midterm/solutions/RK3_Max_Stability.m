clear all;
close all;

A = [ 0 0 0; 1/4 1/4 0; 0 1 0];
b = [1/6 ; 2/3; 1/6];
h = ones(3,1);
I = eye(length(h));

S = @(z)  abs(det(I- z*A + z*h*b')./det(I- z*A))-1;

zz = fsolve(S,-4);

T = 10;
dt = zz/-16
DT = [3e-2,.3,dt-.01,dt+.01];

NSTEPS  = floor(T./DT)+1;

IOSTEPS = [1,1,1,1,1];

f = @(y,t) [-1 3 -5 7; 0 -2 4 -6; 0 0 -4 6; 0 0 0 -16]*y;

y0 = [1;1;1;1]

A = [-1 3 -5 7; 0 -2 4 -6; 0 0 -4 6; 0 0 0 -16];
I=1;
for i = 1:length(DT)

[yk,tk] = RK3_Method_Implicit_Linear(f,A,y0,DT(i),NSTEPS(i),IOSTEPS(i));

figure(I)
clf
plot(tk,yk(1,:))
hold on
plot(tk,yk(2,:))
plot(tk,yk(3,:))
plot(tk,yk(4,:))

ylabel('$y$','interpreter','latex','fontsize',14)
xlabel('$t$','fontsize',14,'interpreter','latex')
title(sprintf('$\\Delta t  = %2.2g$',DT(i)),'interpreter','latex','fontsize',14)
h=legend('$y_1(t)$','$y_2(t)$','$y_3(t)$','$y_4(t)$','interpreter','latex','fontsize',14);
set(h,'Location','NorthEast')'
print(gcf,sprintf('%RK3_AM2_Comparisong.png',I),'-dpng','-r300'); 
%set(gca,'fontsize',22)
%set(gca,'Xtick',[0 2 4 6 8 10],'Ytick',[1e-50 1e-40 1e-30 1e-20 1e-10 1e0],'fontsize',12)
%axis([0 10 1e-50 1e0])
grid off
I=I+1;
end



