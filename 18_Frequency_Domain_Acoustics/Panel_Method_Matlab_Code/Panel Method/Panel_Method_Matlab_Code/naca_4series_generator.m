%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%  naca_4series_generator.m - Create surface panelization for a NACA4 Series Airfoil    %
%                                                                      %
%                                                                      %
%  Input list:                                                         %
%                                                                      %
%  naca4   -  NACA 4 Series Airfoil Denomination                       %
%  npanel  -  Number of panels on the airfoil.  The number of nodes    %
%              is equal to npanel+1, and the ith panel goes from node  %
%              i to node i+1                                           %
%                                                                      %
%  Output list:                                                        %
%                                                                      %
%  x       -  Vector of x coordinates of the surface nodes             %
%  y       -  Vector of y coordinates of the surface nodes             %
%                                                                      %
%  Written by: Matthew Clarke                                          %
%              Department of Aerospace Engineering                     %
%              University of Illinois, Urbana-Champaign                % 
%              maclarke@illinois.edu                                   %
%                                                                      %
%  Last Modified: Wed July 2023                                        %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x,y] = naca_4series_generator(naca4,npanel)

%
% naca airfoil digits
%
 
n1 = str2num(naca4(4:4));
n2 = str2num(naca4(3:3));
n3 = str2num(naca4(2:2));
n4 = str2num(naca4(1:1));

%
% maximum camber, thickness and location of maximum camber
%

m = n4 / 100;
p = n3 / 10;
t = (n2*10 +n1) /100;

%
% compute thickness and camber distributions
%

if mod(npanel,2) ~= 0
    sprintf('Please choose an even number of panels');
    sprintf('Exiting...');
    exit;
end

nside = npanel / 2 +1;

%
% bunching parameter
%

an    = 1.5;
anp   = an +1;

%
% camber distribution
%

for i=1:nside
    frac  = (i -1)/(nside -1);
    xx(i) = 1 -anp*frac*(1 -frac)^an -(1 -frac)^anp;
    yt(i) = ( 0.29690*sqrt(xx(i)) -0.12600*xx(i)     ...
             -0.35160*xx(i)^2      +0.28430*xx(i)^3  ... 
             -0.10150*xx(i)^4) * t / 0.20;
    if xx(i) < p
        yc(i) = m/p^2 * (2*p*xx(i) -xx(i)^2);
    else
        yc(i) = m/(1 -p)^2 * ((1 -2*p) + 2*p*xx(i)-xx(i)^2);
    end
end

%
% airfoil shape = camber + thickness (you can add the thickness distribution normal to the camberline
%                                     if you wish to be more accurate)
%

for i=1:nside
    x(nside+i-1) = xx(i);
    x(nside-i+1) = xx(i);
    y(nside+i-1) = yc(i) +yt(i);
    y(nside-i+1) = yc(i) -yt(i);

end

return
