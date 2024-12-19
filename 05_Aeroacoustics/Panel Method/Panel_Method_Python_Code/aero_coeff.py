########################################################################
#                                                                      #
#  aero_coeff.py - Compute airfoil force and moment coefficients about  #
#                  the quarter chord point                             #
#                                                                      #
#                                                                      #
#  Input list:                                                         #
#                                                                      #
#  x       -  Vector of x coordinates of the surface nodes             #
#  y       -  Vector of y coordinates of the surface nodes             #
#  cp      -  Vector of coefficients of pressure at the nodes          #
#  al      -  Angle of attack in radians                               #
#  npanel  -  Number of panels on the airfoil                          #
#                                                                      #
#  Output list:                                                        #
#                                                                      #
#  cl      -  Airfoil lift coefficient                                 #
#  cd      -  Airfoil drag coefficient                                 #
#  cm      -  Airfoil moment coefficient about the c/4                 #
#                                                                      #
#  Written by: Matthew Clarke                                          #
#              Department of Aerospace Engineering                     #
#              University of Illinois, Urbana-Champaign                #
#              maclarke@illinois.edu                                   #
#                                                                      #
#  Last Modified: Wed Aug 16 2023                                      #
#                                                                      #
########################################################################

import numpy as np

def aero_coeff(x, y, cp, al, npanel):

    # ---------------------------------------------------------------------------  
    # STEP 3.1 allocate all necessary arrays
    # ---------------------------------------------------------------------------      
    Cn = 0
    Ca = 0
    Cm = 0

    for i in range(npanel):
        # compute dx 
        dx = x[i+1] - x[i]
        
        # compute dy 
        dy = y[i+1] - y[i]

        # compute xa at quarter chord
        xa = 0.5 * (x[i+1] + x[i]) - 0.25 

        # compute ya          
        ya = 0.5 * (y[i+1] + y[i])

        # compute dCn            
        dCn = -cp[i] * dx
        
        # compute dCa    
        dCa = cp[i] * dy
        
        # recursively compute Cn, Ca and Cm
        Cn += dCn
        Ca += dCa
        Cm += - (dCn * xa) + (dCa * ya) 
        
    # compute Cl
    Cl  = Cn * np.cos(al) - Ca * np.sin(al)
    
    # compute Cd
    Cd  = Cn * np.sin(al) + Ca * np.cos(al) 

    return Cl, Cd, Cm