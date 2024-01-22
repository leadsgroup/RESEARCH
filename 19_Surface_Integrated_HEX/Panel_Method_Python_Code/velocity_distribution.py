########################################################################
#                                                                      #
# velocity_distribution.py - Compute the tangential velocity 
#              distribution at the midpoint of each panel              #
#                                                                      #
#  Input list:                                                         #
#                                                                      #
#  lambda_gamma      -  Vector of source/sink and vortex strengths     #
#  x       -  Vector of x coordinates of the surface nodes             #
#  y       -  Vector of y coordinates of the surface nodes             #
#  xbar    -  X-coordinate of the midpoint of each panel               #
#  ybar    -  X-coordinate of the midpoint of each panel               #
#  sin_theta      -  Sin(theta) for each panel                         #
#  cos_theta      -  Cos(theta) for each panel                         #
#  alpha      -  Angle of attack in radians                            #
#  npanel  -  Number of panels on the airfoil                          #
#                                                                      #
#  Output list:                                                        #
#                                                                      #
#  vt      -  Vector of tangential velocities                          #
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
from numpy import sin, cos, sqrt, arctan2, array, log, pi, cross 

def velocity_distribution(lambda_gamma, x, y, xbar, ybar, sin_theta, cos_theta, alpha, npanel):

    # ---------------------------------------------------------------------------  
    # STEP 4.1 precompute common terms 
    # ---------------------------------------------------------------------------  
    pi2inv = 1 / (2 * pi)
    
    
    # ---------------------------------------------------------------------------  
    # STEP 7.2 allocate all necessary arrays
    # ---------------------------------------------------------------------------   
    vt     = np.zeros_like(cos_theta)
    
    
    # ---------------------------------------------------------------------------  
    # Step 7.3 flow tangency boundary condition - source distribution
    # ---------------------------------------------------------------------------   
    for i in range(npanel):
        # compute the velocity distribution on all panels on the airfoil here
        vt[i] = cos_theta[i] * cos(alpha) + sin_theta[i] * sin(alpha)

        for j in range(npanel):
            # compute r_ij
            r_i_j      = sqrt((xbar[i] - x[j]) ** 2 + (ybar[i] - y[j]) ** 2) 

            # compute r_ij+1            
            r_i_jplus1 = sqrt((xbar[i] - x[j + 1]) ** 2 + (ybar[i] - y[j + 1]) ** 2)

            # compute beta_ij
            if i == j:
                beta_i_j = pi
            else: 
                v1 = array([[xbar[i]], [ybar[i]]]) - array([[x[j]], [y[j]]])
                v2 = array([[xbar[i]], [ybar[i]]]) - array([[x[j + 1]], [y[j + 1]]])
                beta_i_j = arctan2(cross(v1.T, v2.T), v1.T @ v2)[0][0] 
                
            # compute vt 
            vt[i] = vt[i] + lambda_gamma[j] * pi2inv * ( (sin_theta[i] * cos_theta[j] - cos_theta[i] * sin_theta[j]) * beta_i_j  - (cos_theta[i] * cos_theta[j] + sin_theta[i] * sin_theta[j]) * log( r_i_jplus1 / r_i_j))  \
                + lambda_gamma[npanel] * pi2inv * ( (sin_theta[i] * cos_theta[j] - cos_theta[i] * sin_theta[j]) * log( r_i_jplus1 / r_i_j) + (cos_theta[i] * cos_theta[j] + sin_theta[i] * sin_theta[j]) * beta_i_j)
    return vt
