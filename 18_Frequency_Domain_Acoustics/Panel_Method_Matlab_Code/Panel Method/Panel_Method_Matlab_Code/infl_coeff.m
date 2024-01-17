%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
%  infl_coeff.m - Compute the matrix of aerodynamic influence          %
%                  coefficients for later use                          %
%                                                                      %
%  Input list:                                                         %
%                                                                      %
%  x       -  Vector of x coordinates of the surface nodes             %
%  y       -  Vector of y coordinates of the surface nodes             %
%  xbar    -  X-coordinate of the midpoint of each panel               %
%  ybar    -  X-coordinate of the midpoint of each panel               %
%  st      -  Sin(theta) for each panel                                %
%  ct      -  Cos(theta) for each panel                                %
%  ainfl   -  Aero influence coefficient matrix                        %
%  npanel  -  Number of panels on the airfoil                          %
%                                                                      %
%  Output list:                                                        %
%                                                                      %
%  ainfl   -  Aero influence coefficient matrix                        %
%                                                                      %
%  Written by: Matthew Clarke                                          %
%              Department of Aerospace Engineering                     %
%              University of Illinois, Urbana-Champaign                % 
%              maclarke@illinois.edu                                   %
%                                                                      %
%  Last Modified: Wed July 2023                                        %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ainfl = infl_coeff(x,y,xbar,ybar,st,ct,ainfl,npanel);

pi2inv = 1 / (2*pi);

%
% Fill the elements of the matrix of aero influence coefficients
%

for i=1:npanel

    %
    % find contribution of the jth panel
    %
  
    %%
    for j=1:npanel
        
        r_i_j = sqrt((xbar(i)- x(j))^2+ (ybar(i)- y(j))^2);
        r_i_jplus1 = sqrt((xbar(i)- x(j+1))^2+ (ybar(i)- y(j+1))^2);
        if i == j
            beta_i_j = pi;
        else 
             v1 = [xbar(i),ybar(i)] - [x(j), y(j)];
             v2 = [xbar(i),ybar(i)] - [x(j+1), y(j+1)];
             beta_i_j = atan2( v1(1)*v2(2) - v1(2)*v2(1),v1(1)*v2(1) + v1(2)*v2(2));
        end
        u_s_star_i_j = -pi2inv*log(r_i_jplus1/r_i_j);
        v_s_star_i_j = pi2inv*beta_i_j;
        ainfl(i,j) = -u_s_star_i_j*(ct(j)*st(i) - st(j)*ct(i)) + v_s_star_i_j*(st(j)*st(i) + ct(j)*ct(i));        
        
        ainfl(i,npanel+1) = ainfl(i,npanel+1) + pi2inv *((ct(i)*ct(j) + st(i)*st(j))*log(r_i_jplus1/r_i_j) - (st(i)*ct(j) - ct(i)*st(j))*beta_i_j); %% correct     
        
        if i == 1 || i == npanel
           ainfl(npanel+1,j) = ainfl(npanel+1,j) +  pi2inv *((st(i)*ct(j) - ct(i)*st(j))*beta_i_j - (ct(i)*ct(j) + st(i)*st(j))*log(r_i_jplus1/r_i_j)); %% correct 
          
           ainfl(npanel+1,npanel+1) = ainfl(npanel+1,npanel+1) + ainfl(i,j);
        end
     end
end
 
if rank(ainfl) ~= npanel+1
    exit
end

return
