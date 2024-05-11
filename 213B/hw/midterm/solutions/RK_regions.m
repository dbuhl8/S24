function RK_regions()

% Lobatto 3C
A_RK31 = [ 0 0 0; 1/4 1/4 0; 0 1 0];
b_RK31 = [1/6 ; 2/3; 1/6];
h_RK31 = ones(3,1);


z4=region(A_RK31,b_RK31,h_RK31);

figure(1)
clf

plot(z4(1,2:end),z4(2,2:end),'g','Linewidth',1.5)
grid
xlabel('{Re}$(z)$','Interpreter','Latex')
ylabel('{Im}$(z)$','Interpreter','Latex')
%set(gca,'Fontsize',22)
axis equal
axs([-6 2 -5 5])

print(gcf,'RK3_Implicit_Stability_Region.png','-dpng','-r300'); 
end




function z=region(A,b,h)

N=500;
z       = linspace(-6,6,N);
[Z1,Z2] = meshgrid(z,z);

Z=Z1+1i*Z2;

I=eye(length(h));
S=zeros(N,N);

for ii=1:N
    for jj=1:N
        S(ii,jj) = det(I- Z(ii,jj)*A + Z(ii,jj)*h*b')./det(I- Z(ii,jj)*A);
    end
end

[z,hh]=contour(Z1,Z2,abs(S)-1,[0 0]);
end
