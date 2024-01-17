########################################################################
#
#  naca_4series_generator.py - Create surface panelization for a NACA4 Series Airfoil
#
#
#  Input list:
#
#  naca4   -  NACA 4 Series Airfoil Denomination
#  npanel  -  Number of panels on the airfoil.  The number of nodes
#             is equal to npanel+1, and the ith panel goes from node
#             i to node i+1
#
#  Output list:
#
#  x       -  Vector of x coordinates of the surface nodes
#  y       -  Vector of y coordinates of the surface nodes
#
#  Written by: Matthew Clarke                                          #
#              Department of Aerospace Engineering                     #
#              University of Illinois, Urbana-Champaign                #
#              maclarke@illinois.edu                                   #
#                                                                      #
#  Last Modified: Wed Aug 30 2023
########################################################################

import numpy as np 

def naca_4series_generator(naca4, npanel):
    # ---------------------------------------------------------------------------  
    # STEP 2.1: naca airfoil digits
    # ---------------------------------------------------------------------------  
    x = np.zeros(npanel + 1)
    y = np.zeros(npanel + 1)
    n1 = int(naca4[3])
    n2 = int(naca4[2])
    n3 = int(naca4[1])
    n4 = int(naca4[0])

    # ---------------------------------------------------------------------------  
    # STEP 2.2: maximum camber, thickness, and location of maximum camber
    # ---------------------------------------------------------------------------   
    m = n4 / 100
    p = n3 / 10
    t = (n2 * 10 + n1) / 100

    # ---------------------------------------------------------------------------  
    # STEP 2.3: compute thickness and camber distributions
    # ---------------------------------------------------------------------------   
    if npanel % 2 != 0:
        raise ValueError("Please choose an even number of panels!") 
    nside = int(npanel / 2 + 1)

    # ---------------------------------------------------------------------------  
    # STEP 2.4: bunching parameter for higher resolution near leading edge and 
    #           trailing edge
    # ---------------------------------------------------------------------------  

    an  = 1.5
    anp = an + 1 
    xx  = np.zeros(nside)
    yt  = np.zeros(nside)
    yc  = np.zeros(nside)

    # ---------------------------------------------------------------------------  
    # STEP 2.4: determine camberline 
    # ---------------------------------------------------------------------------      
    for i in range(nside):
        frac = i / (nside - 1)
        xx[i] = 1 - anp * frac * (1 - frac) ** an - (1 - frac) ** anp
        yt[i] = (0.29690 * np.sqrt(xx[i]) - 0.12600 * xx[i]
                 - 0.35160 * xx[i] ** 2 + 0.28430 * xx[i] ** 3
                 - 0.10150 * xx[i] ** 4) * t / 0.2 
        if xx[i] < p:
            yc[i] = m / p ** 2 * (2 * p * xx[i] - xx[i] ** 2)
        else:
            yc[i] = m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * xx[i] - xx[i] ** 2)
        
    # ---------------------------------------------------------------------------  
    # STEP 2.5: determine airfoil shape = camber + thickness 
    #           you can add the thickness distribution normal to the camberline
    #           if you wish to be more accurate
    # ---------------------------------------------------------------------------   
    for i in range(nside):
        x[nside + i - 1] = xx[i]
        x[nside - i - 1] = xx[i]
        y[nside + i - 1] = yc[i] + yt[i]
        y[nside - i - 1] = yc[i] - yt[i]

    return x, y