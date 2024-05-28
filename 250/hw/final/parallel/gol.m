clear; close all; clc;

G = readmatrix('gol.dat', 'Delimiter', ' ');
p = readmatrix('params.dat');
p = rmmissing(p);

G = reshape(transpose(G), p);

% The routine to make the movie is taken from 
% https://www.mathworks.com/help/matlab/ref/movie.html#d126e1052750

h = figure;
frames = p(3);
M(frames) = struct('cdata',[],'colormap',[]);
h.Visible = 'off';
for i = 1:frames;
    imagesc(G(:,:,i))
    drawnow
    M(i) = getframe; 
end
h.Visible = 'on';
movie(M,1,10);

