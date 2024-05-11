function [y,t] = RK4_Method(f,y0,DT,NSTEPS,IOSTEPS)

% This function computes the numerical solution to 
% the system of ordinary differential equations
% 
%          dy/dt = f(y,t) y(0)=y0
%
% Input f  -> function handle representing f(y,t)
%       y0 -> column vector of inintial conditionn
%       DT -> delta t
%   NSTEPS -> total number of time steps ( T=DT*NSTEPS)
%  IOSTEPS -> input/output steps (one time snapshots is
%             saved into the output matrix every IOSTEPS
%             steps.
%
% Output y -> matrix that collects the time snapshots of 
%             the solution columnwise
%        t -> row vector collecting the times at which the 
%             solution is saved

t=0;
y=y0;


for ii=1:NSTEPS
    
    ts = (ii-1)*DT;
    
    K1 = f(y0,ts);

	K2 = f(y0 + DT * 0.5 *K1, ts + 0.5* DT);


	K3 = f(y0 + DT * 0.5 *K2, ts + 0.5* DT);
    
    K4 = f(y0 + DT * K3, ts+ DT);
    
    y1 = y0 + (DT / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
    
  
    if mod(ii,IOSTEPS)==0
        y =[y y1];
        t = [ t ts+DT];
    end
   y0=y1;
   
end

    
end



