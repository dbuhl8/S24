close all; clear; clc; 

p = readmatrix('../fortran/params.dat');
u = readmatrix('../fortran/u.dat');
x = readmatrix('../fortran/x.dat');
y = readmatrix('../fortran/y.dat');
t = readmatrix('../fortran/t.dat');
%t = rmmissing(t);
%y = rmmissing(y);
%x = rmmissing(x);
p = rmmissing(p);
%u = rmmissing(u);

U = reshape(transpose(u), p);

fn = 1;
%v=VideoWriter('u_movie');
%v.FrameRate = 400;
%open(v)

%h = figure;
%frames = p(3);
%M(frames) = struct('cdata',[],'colormap',[]);
%h.Visible = 'off';
%xlabel('X');
%ylabel('Y');
%zlabel('U');
%title('Numerical Solution t = [0,2]');
%for i = 1:frames;
    %s = surf(x, y, transpose(U(:,:,i)), 'FaceColor','interp');
    %s.EdgeColor = 'none';
    %drawnow
    %M(i) = getframe; 
    %writeVideo(v, getframe(gcf));
%end
%close(v)
%h.Visible = 'on';

vals = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.];
[lia, ind] = ismember(vals, t(:,3));

% Contour Plots for FD
for i = 1:numel(ind)
    figure(fn);
    fn = fn + 1;
    %s = surf(x, y, transpose(U(:,:,ind(i))),'FaceColor','interp');
    %s.EdgeColor = 'none';
    contour(x, y, transpose(U(:,:,ind(i))),20);
    xlabel('X');
    ylabel('Y');
    view(0,90)
    title(strcat('FD Contour of U(t=',string(t(ind(i),3)),')'));
    colorbar
    int_vals(i) = (sum(U(:,:,ind(i)), 'all')*4*pi^2)/((p(1)-2)^2);
end 

% Contour Plots for Characteristics
char_u = readmatrix('../fortran/char_u.dat');
char_u = reshape(transpose(char_u), [80,80,9]);
cx = readmatrix('../fortran/char_x.dat');
cy = readmatrix('../fortran/char_y.dat');

for i = 1:numel(ind)
    figure(fn);
    fn = fn +1;
    %s = surf(cx, cy, transpose(char_u(:,:,i)),'FaceColor','interp');
    %s.EdgeColor = 'none';
    contour(x, y, transpose(char_u(:,:,i)),20);
    view(0,90)
    xlabel('X');
    ylabel('Y');
    title(strcat('Char Contour of U(t=',string(vals(i)),')'));
    colorbar
end 


% MSE Plots
vals = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.];
[lia, ind] = ismember(vals, t(:,3));

for i = 1:numel(ind)
    mpe_val(i) = max(abs(transpose(U(:,:,ind(i))-transpose(char_u(:,:,i)))),[],'all')
end 

% Max Pointwise Error Plots
figure(fn);
fn=fn+1;
plot(vals,mpe_val)
xlabel('t')
ylabel('MPE')
title('Maximum Pointwise Error v t')

% Integral Plots
figure(fn);
fn = fn +1;
char_int = readmatrix('../fortran/char_int.dat');
plot(vals, char_int(:,2))
xlabel('t')
ylabel('int_u')
title('Volume Integral of U v time, Char')

figure(fn);
fn = fn+1;
plot(vals, int_vals)
xlabel('t')
ylabel('int_u')
title('Volume Integral of U v time, FD')

