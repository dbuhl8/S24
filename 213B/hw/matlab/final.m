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
title('Numerical Solution')
figure(2);
s = surf(xx, tt, ua, 'FaceColor','interp');
s.EdgeColor = 'none';
title('Analytical Solution')
colorbar
figure(3);
plot(x,u(1,:),'bx',x,ua(1,:),'r')
xlabel('r')
ylabel('u')
title('Initial Condition')
legend()

figure(4);
plot(x,u(help(1),:),'bx',x,ua(help(1),:),'r')
xlabel('r')
ylabel('u')
title('Final Solution')
legend()

% Make a movie showing the evolution of the two together on a 2d plot

h = figure();
%v=VideoWriter('u_evol');
%v.FrameRate = 1000;
%open(v)
frames = help(1);
M(frames) = struct('cdata',[],'colormap',[]);
h.Visible = 'off';
xlabel('R');
ylabel('U');
ylim([0,5]);
title('Numerical Solution t = [0,2]');
for i = 1:(frames-1)/100;
    plot(x,u(i*100,:),'bx',x,ua(i*100,:),'r')
    ylim([0,5]);
    drawnow
    M(i) = getframe; 
    %writeVideo(v, M(i));
end
%close(v)
h.Visible = 'on';
movie(M,1,100)

