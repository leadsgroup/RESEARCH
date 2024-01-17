%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%  hess_smith.m - Main function for the computation of the             %
%                 incompressible, inviscid flow over an airfoil of     %
%                 arbitrary shape using the Hess-Smith panel method.   %
%                                                                      %
%  References: "An introduction to theoretical and computational       %
%               aerodynamics", J. Moran, Wiley, 1984                   %
%                                                                      %
%  Input list:                                                         %
%                                                                      %
%  naca4   -  NACA 4 Series Airfoil Denomination                       %
%  alpha   -  Airfoil angle of attack                                  %
%  npanel  -  Number of panels on the airfoil.  The number of nodes    %
%              is equal to npanel+1, and the ith panel goes from node  %
%              i to node i+1                                           %
%                                                                      %
%  Output list:                                                        %
%                                                                      %
%  cl      -  Airfoil lift coefficient                                 %
%  cd      -  Airfoil drag coefficient                                 %
%  cm      -  Airfoil moment coefficient about the c/4                 %
%  x       -  Vector of x coordinates of the surface nodes             %
%  y       -  Vector of y coordinates of the surface nodes             %
%  cp      -  Vector of coefficients of pressure at the nodes          %
%                                                                      %
%                                                                      %
%  Written by: Matthew Clarke                                          %
%              Department of Aerospace Engineering                     %
%              University of Illinois, Urbana-Champaign                % 
%              maclarke@illinois.edu                                   %
%                                                                      %
%  Last Modified: Wed July 2023                                        %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cl,cd,cm,x,y,cp] = hess_smith(naca4,alpha,npanel)

%
% allocate all necessary arrays
%

x  = zeros(npanel+1,1);
y  = zeros(npanel+1,1);
cp = zeros(npanel+1,1);

%
% generate airfoil surface panelization
%

[x,y] = naca_4series_generator(naca4,npanel);

%
% generate panel geometry data for later use
%

l    = zeros(npanel,1);
st   = zeros(npanel,1);
ct   = zeros(npanel,1);
xbar = zeros(npanel,1);
ybar = zeros(npanel,1);

[l,st,ct,xbar,ybar] = panel_geometry(x,y,npanel);

%
% compute matrix of aerodynamic influence coefficients
%

ainfl = zeros(npanel+1);

ainfl = infl_coeff(x,y,xbar,ybar,st,ct,ainfl,npanel);

%
% compute right hand side vector for the specified angle of attack
%

b  = zeros(npanel+1,1);

al = alpha * pi / 180;

for i=1:npanel
    b(i) = st(i)*cos(al) -sin(al)*ct(i);
end

b(npanel+1) = -(ct(1)     *cos(al) +st(1)     *sin(al)) ...
              -(ct(npanel)*cos(al) +st(npanel)*sin(al));
          
%
% solve matrix system for vector of q_i and gamma
%

qg = inv(ainfl) * b;

%
% compute the tangential velocity distribution at the midpoint of panels
%

vt = velocity_distribution(qg,x,y,xbar,ybar,st,ct,al,npanel);

%
% compute pressure coefficient
%

cp = 1 -vt.^2;

%
% compute force coefficients
%

[cl,cd,cm] = aero_coeff(x,y,cp,al,npanel);

%
% plot the output
%

subplot(2,1,1),plot(xbar,-cp),xlabel('x/c'),ylabel('Cp'),title('Coefficient of Pressure Distribution'),grid

subplot(2,1,2),plot(xbar,ybar,xbar,ybar,'o'),xlabel('x/c'),ylabel('y/c'),title('Airfoil Geometry'),axis('equal'),grid

return
