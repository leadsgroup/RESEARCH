%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%  panel_geometry.m - Compute airfoil surface panelization parameters  %
%                      for later use in the computation of the matrix  %
%                      of influence coefficients.                      %
%                                                                      %
%                                                                      %
%  Input list:                                                         %
%                                                                      %
%  x       -  Vector of x coordinates of the surface nodes             %
%  y       -  Vector of y coordinates of the surface nodes             %
%  npanel  -  Number of panels on the airfoil                          %
%                                                                      %
%  Output list:                                                        %
%                                                                      %
%  l       -  Panel lenghts                                            %
%  st      -  Sin(theta) for each panel                                %
%  ct      -  Cos(theta) for each panel                                %
%  xbar    -  X-coordinate of the midpoint of each panel               %
%  ybar    -  X-coordinate of the midpoint of each panel               %
%                                                                      %
%  Written by: Matthew Clarke                                          %
%              Department of Aerospace Engineering                     %
%              University of Illinois, Urbana-Champaign                % 
%              maclarke@illinois.edu                                   %
%                                                                      %
%  Last Modified: Wed July 2023                                        %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [l,st,ct,xbar,ybar] = panel_geometry(x,y,npanel);

%
% compute various geometrical quantities
%

for i=1:npanel
    l   (i) = sqrt((x(i+1) -x(i))^2 +(y(i+1) -y(i))^2);
    st  (i) = (y(i+1) -y(i))/l(i);
    ct  (i) = (x(i+1) -x(i))/l(i);
    xbar(i) = (x(i+1) +x(i))/2;
    ybar(i) = (y(i+1) +y(i))/2;
end

return




