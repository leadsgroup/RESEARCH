########################################################################
#                                                                      #
#  infl_coeff.py - Compute the matrix of aerodynamic influence          #
#                  coefficients for later use                          #
#                                                                      #
#  Input list:                                                         #
#                                                                      #
#  x       -  Vector of x coordinates of the surface nodes             #
#  y       -  Vector of y coordinates of the surface nodes             #
#  xbar    -  X-coordinate of the midpoint of each panel               #
#  ybar    -  X-coordinate of the midpoint of each panel               #
#  st      -  Sin(theta) for each panel                                #
#  ct      -  Cos(theta) for each panel                                #
#  A       -  Aero influence coefficient matrix                        #
#  npanel  -  Number of panels on the airfoil                          #
#                                                                      #
#  Output list:                                                        #
#                                                                      #
#  A   -  Aero influence coefficient matrix                            #
#                                                                      #
#  Written by: Matthew Clarke                                          #
#              Department of Aerospace Engineering                     #
#              University of Illinois, Urbana-Champaign                # 
#              maclarke@illinois.edu                                   #
#                                                                      #
#  Last Modified: Wed Aug 30 2023                                      #
#                                                                      #
########################################################################
import numpy as np
from numpy import array, arctan2, cross


def infl_coeff(x, y, xbar, ybar, st, ct, A, npanel):

    # ---------------------------------------------------------------------------  
    # STEP 4.1 precompute common terms 
    # ---------------------------------------------------------------------------      
    pi2inv = 1 / (2 * np.pi)

    # ---------------------------------------------------------------------------     
    # Step 4.2 Define the elements of the matrix of A aero influence coefficients
    # ---------------------------------------------------------------------------     
    # ith panel (major loop)
    for i in range(npanel): 
        # find contribution of the jth panel (inner loop)
        
        for j in range(npanel):
            # compute r_ij
            r_i_j      = np.sqrt((xbar[i] - x[j]) ** 2 + (ybar[i] - y[j]) ** 2)
            
            # compute r_ij+1
            r_i_jplus1 = np.sqrt((xbar[i] - x[j + 1]) ** 2 + (ybar[i] - y[j + 1]) ** 2)
            
            # compute beta_ij
            if i == j:
                beta_i_j = np.pi
            else: 
                v1       = array([[xbar[i]], [ybar[i]]]) - array([[x[j]], [y[j]]])
                v2       = array([[xbar[i]], [ybar[i]]]) - array([[x[j + 1]], [y[j + 1]]])
                beta_i_j = arctan2(cross(v1.T, v2.T), np.dot(v1.T,v2))[0][0]  
                
            # compute u_s_star_i_j
            u_s_star_i_j = -pi2inv * np.log(r_i_jplus1 / r_i_j) 
            
            # compute v_s_star_i_j            
            v_s_star_i_j = pi2inv * beta_i_j 
            # compute all elements of A, the influence coefficient matrix   
            A[i, j]      = -u_s_star_i_j * (ct[j] * st[i] - st[j] * ct[i]) + v_s_star_i_j * (st[j] * st[i] + ct[j] * ct[i])  
            A[i, npanel] = A[i, npanel] + pi2inv * ( (ct[i] * ct[j] + st[i] * st[j]) * np.log( r_i_jplus1 / r_i_j) - (st[i] * ct[j] - ct[i] * st[j]) * beta_i_j)    
            if i == 0 or i == npanel - 1:
                # kutta condition , therefore only applies for the first and last panel 
                A[npanel, j]      = A[npanel, j] + pi2inv * ( (st[i] * ct[j] - ct[i] * st[j]) * beta_i_j - ( ct[i] * ct[j] + st[i] * st[j]) * np.log( r_i_jplus1 / r_i_j))   
                A[npanel, npanel] = A[npanel, npanel] + A[i, j]

    # check to see if matrix is singular 
    if np.linalg.det(A) == 0:
        raise ValueError("Matrix is singular") 
    return A
