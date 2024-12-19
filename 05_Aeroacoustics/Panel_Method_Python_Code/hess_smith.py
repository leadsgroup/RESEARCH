########################################################################
#                                                                      #
#  hess_smith.py - Main function for the computation of the             #
#                 incompressible, inviscid flow over an airfoil of     #
#                 arbitrary shape using the Hess-Smith panel method.   #
#                                                                      #
#  References: "An introduction to theoretical and computational       #
#               aerodynamics", J. Moran, Wiley, 1984                   #
#                                                                      #
#  Input list:                                                         #
#                                                                      #
#  naca4   -  NACA 4 Series Airfoil Denomination                       #
#  alpha   -  Airfoil angle of attack                                  #
#  npanel  -  Number of panels on the airfoil.  The number of nodes    #
#              is equal to npanel+1, and the ith panel goes from node  #
#              i to node i+1                                           #
#                                                                      #
#  Output list:                                                        #
#                                                                      #
#  cl      -  Airfoil lift coefficient                                 #
#  cd      -  Airfoil drag coefficient                                 #
#  cm      -  Airfoil moment coefficient about the c/4                 #
#  x       -  Vector of x coordinates of the surface nodes             #
#  y       -  Vector of y coordinates of the surface nodes             #
#  cp      -  Vector of coefficients of pressure at the nodes          #
#                                                                      #
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
from panel_geometry         import panel_geometry
from velocity_distribution  import velocity_distribution
from infl_coeff             import infl_coeff
from aero_coeff             import aero_coeff
 
def hess_smith(x,y,alpha): 
    # ---------------------------------------------------------------------------  
    # STEP 1: allocate all necessary arrays 
    # ---------------------------------------------------------------------------  
    npanel       = len(x) - 1 
    cp           = np.zeros(npanel + 1) 
    l            = np.zeros(npanel)
    sin_theta    = np.zeros(npanel)
    cos_theta    = np.zeros(npanel)
    xbar         = np.zeros(npanel)
    ybar         = np.zeros(npanel) 

    # ---------------------------------------------------------------------------  
    # STEP 3: generate panel geometry data for later use 
    # ---------------------------------------------------------------------------  
    [l, sin_theta, cos_theta, xbar, ybar] = panel_geometry(x, y, npanel)


    # ---------------------------------------------------------------------------  
    # STEP 4: compute matrix of aerodynamic influence coefficients 
    # ---------------------------------------------------------------------------  
    A = np.zeros((npanel + 1, npanel + 1)) 
    A = infl_coeff(x, y, xbar, ybar, sin_theta, cos_theta, A, npanel) 

    # --------------------------------------------------------------------------- 
    # STEP 5: compute right hand side vector for the specified angle of attack
    # ---------------------------------------------------------------------------  
    b = np.zeros(npanel + 1) 
    al = alpha * np.pi / 180

    for i in range(npanel):
        b[i] = sin_theta[i] * np.cos(al) - np.sin(al) * cos_theta[i] 
    b[npanel] = -(cos_theta[0] * np.cos(al) + sin_theta[0] * np.sin(al)) - (cos_theta[npanel-1] * np.cos(al) + sin_theta[npanel-1] * np.sin(al))


    # ---------------------------------------------------------------------------  
    # STEP 6: solve matrix system for vector of lambda_i and gamma 
    # --------------------------------------------------------------------------- 
    lambda_gamma = np.linalg.inv(A) @ b


    # ---------------------------------------------------------------------------  
    # STEP 7: compute the tangential velocity distribution at the midpoint of panels 
    # --------------------------------------------------------------------------- 
    vt = velocity_distribution(lambda_gamma, x, y, xbar, ybar, sin_theta, cos_theta, al, npanel)


    # ---------------------------------------------------------------------------  
    # STEP 8: compute pressure coefficient  
    # --------------------------------------------------------------------------- 
    cp = 1 - vt ** 2 

    # ---------------------------------------------------------------------------  
    # STEP 9: compute force coefficients 
    # --------------------------------------------------------------------------- 
    cl, cd, cm = aero_coeff(x, y, cp, al, npanel) 

    return cl, cd, cm,cp, xbar, ybar,vt,cos_theta,sin_theta
