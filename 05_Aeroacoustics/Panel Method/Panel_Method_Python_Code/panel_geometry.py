########################################################################
#                                                                      #
#  panel_geometry.py - Compute airfoil surface panelization parameters #
#                      for later use in the computation of the matrix  #
#                      of influence coefficients.                      #
#                                                                      #
#                                                                      #
#  Input list:                                                         #
#                                                                      #
#  x       -  Vector of x coordinates of the surface nodes             #
#  y       -  Vector of y coordinates of the surface nodes             #
#  npanel  -  Number of panels on the airfoil                          #
#                                                                      #
#  Output list:                                                        #
#                                                                      #
#  l       -  Panel lenghts                                            #
#  sin_theta      -  Sin(theta) for each panel                                #
#  cos_theta      -  Cos(theta) for each panel                                #
#  xbar    -  X-coordinate of the midpoint of each panel               #
#  ybar    -  X-coordinate of the midpoint of each panel               #
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

def panel_geometry(x, y, npanel):
    # ---------------------------------------------------------------------------  
    # STEP 3.1 allocate all necessary arrays
    # ---------------------------------------------------------------------------  
    l           = np.zeros(npanel)
    sin_theta   = np.zeros(npanel)
    cos_theta   = np.zeros(npanel)
    xbar        = np.zeros(npanel)
    ybar        = np.zeros(npanel)

    # ---------------------------------------------------------------------------  
    # STEP 3.2 compute various geometrical quantities
    # ---------------------------------------------------------------------------  
    for i in range(npanel):
        l[i]    = np.sqrt((x[i + 1] - x[i]) ** 2 + (y[i + 1] - y[i]) ** 2)
        sin_theta[i]   = (y[i + 1] - y[i]) / l[i]
        cos_theta[i]   = (x[i + 1] - x[i]) / l[i]
        xbar[i] = (x[i + 1] + x[i]) / 2
        ybar[i] = (y[i + 1] + y[i]) / 2
     
    return l, sin_theta, cos_theta, xbar, ybar
