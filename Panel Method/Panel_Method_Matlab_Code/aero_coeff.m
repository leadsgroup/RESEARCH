%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%  aero_coeff.m - Compute airfoil force and moment coefficients about  %
%                  the quarter chord point                             %
%                                                                      %
%                                                                      %
%  Input list:                                                         %
%                                                                      %
%  x       -  Vector of x coordinates of the surface nodes             %
%  y       -  Vector of y coordinates of the surface nodes             %
%  cp      -  Vector of coefficients of pressure at the nodes          %
%  al      -  Angle of attack in radians                               %
%  npanel  -  Number of panels on the airfoil                          %
%                                                                      %
%  Output list:                                                        %
%                                                                      %
%  cl      -  Airfoil lift coefficient                                 %
%  cd      -  Airfoil drag coefficient                                 %
%  cm      -  Airfoil moment coefficient about the c/4                 %
%                                                                      %
%  Written by: Matthew Clarke                                          %
%              Department of Aerospace Engineering                     %
%              University of Illinois, Urbana-Champaign                % 
%              maclarke@illinois.edu                                   %
%                                                                      %
%  Last Modified: Wed July 2023                                        %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cl,cd,cm] = aero_coeff(x,y,cp,al,npanel);


cl = 0;
cd = 0;
cm = 0;

for i=1:npanel
    dx  = x(i+1) -x(i);
    dy  = y(i+1) -y(i);
    xa  = 0.5*(x(i+1) +x(i)) -0.25;
    ya  = 0.5*(y(i+1) +y(i));
    dcl = -cp(i)*dx;
    dcd = cp(i)*dy;
    cl  = cl +dcl;
    cd  = cd +dcd;
    cm  = cm +dcd*ya -dcl*xa;
end

dcl = cl*cos(al) -cd*sin(al);
cd  = cl*sin(al) +cd*cos(al);
cl  = dcl;


return
