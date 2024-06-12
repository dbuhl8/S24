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
fn = 4
figure(fn);
fn = fn +1;
plot(x,ua(1,:), 'r', x, ua(10001,:),'b', x, ua(20001,:),'m', x, ua(40001,:), 'c')
xlabel('r')
ylabel('u')
title('Analytical Solution at t = [0, 0.5, 1, 2]')
legend({'t = 0', 't = 0.5', 't = 1', 't = 2'})

%figure(fn);
%fn = fn +1;
%plot(x,ua(20001,:),'b')
%xlabel('r')
%ylabel('u')
%title('Analytical Solution at t = 1')
%legend()

%figure(fn);
%fn = fn +1;
%plot(x,ua(40001,:),'b')
%xlabel('r')
%ylabel('u')
%title('Analytical Solution at t = 2')
%legend()


figure(fn);
plot(x,u(help(1),:),'bx',x,ua(help(1),:),'r')
xlabel('r')
ylabel('u')
title('Final Solution')
legend()



% Make a movie showing the evolution of the two together on a 2d plot

%h = figure();
%v=VideoWriter('u_evol');
%v.FrameRate = 1000;
%open(v)
%frames = help(1);
%M(frames) = struct('cdata',[],'colormap',[]);
%h.Visible = 'off';
%xlabel('R');
%ylabel('U');
%ylim([0,5]);
%title('Numerical Solution t = [0,2]');
%for i = 1:(frames-1)/100;
    %plot(x,u(i*100,:),'bx',x,ua(i*100,:),'r')
    %ylim([0,5]);
    %drawnow
    %M(i) = getframe; 
    %writeVideo(v, M(i));
%end
%close(v)
%h.Visible = 'on';
%movie(M,1,100)

