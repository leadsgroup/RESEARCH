%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% velocity_distribution.m - Compute the tangential velocity distribution at the      %
%              midpoint of each panel                                  %
%                                                                      %
%  Input list:                                                         %
%                                                                      %
%  qg      -  Vector of source/sink and vortex strengths               %
%  x       -  Vector of x coordinates of the surface nodes             %
%  y       -  Vector of y coordinates of the surface nodes             %
%  xbar    -  X-coordinate of the midpoint of each panel               %
%  ybar    -  X-coordinate of the midpoint of each panel               %
%  st      -  Sin(theta) for each panel                                %
%  ct      -  Cos(theta) for each panel                                %
%  al      -  Angle of attack in radians                               %
%  npanel  -  Number of panels on the airfoil                          %
%                                                                      %
%  Output list:                                                        %
%                                                                      %
%  vt      -  Vector of tangential velocities                          %
%                                                                      %
%  Written by: Matthew Clarke                                          %
%              Department of Aerospace Engineering                     %
%              University of Illinois, Urbana-Champaign                % 
%              maclarke@illinois.edu                                   %
%                                                                      %
%  Last Modified: Wed July 2023                                        %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vt = velocity_distribution(qg,x,y,xbar,ybar,st,ct,al,npanel);

pi2inv = 1 / (2*pi);
%
% flow tangency boundary condition - source distribution
%

for i=1:npanel
   %<compute the velocity distribution on all panels on the airfoil here>;
   vt(i) = (ct(i)*cos(al) + st(i)*sin(al)); 

   for j = 1:npanel
        r_i_j = sqrt((xbar(i)- x(j))^2+ (ybar(i)- y(j))^2);
        r_i_jplus1 = sqrt((xbar(i)- x(j+1))^2+ (ybar(i)- y(j+1))^2);
        
        if i == j
            beta_i_j = pi;
        else 
             v1 = [xbar(i),ybar(i)] - [x(j), y(j)];
             v2 = [xbar(i),ybar(i)] - [x(j+1), y(j+1)];
             beta_i_j = atan2(v1(1)*v2(2) - v1(2)*v2(1),v1(1)*v2(1) + v1(2)*v2(2));
        end 
        
        vt(i) = vt(i) + qg(j)*pi2inv*((st(i)*ct(j) - ct(i)*st(j))*beta_i_j - (ct(i)*ct(j)...
            + st(i)*st(j))*log(r_i_jplus1/r_i_j)) + qg(npanel+1)*pi2inv*((st(i)*ct(j)  ...
            - ct(i)*st(j))*log(r_i_jplus1/r_i_j) + (ct(i)*ct(j) + st(i)*st(j))*beta_i_j);
       
   end 
   
end

return