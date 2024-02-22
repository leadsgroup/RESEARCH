#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 15:05:20 2024

@author: nn18
"""

"""
This will solve the boundary layer over a flat plate having a pressure gradient (acquired from the inviscid panel method)
"""
# RCAIDE imports
#from RCAIDE.Methods.Aerodymanics.Airfoil_Panel_Method import thwaites_method
#from Research.RK4 import RK4

# Python Imports
import numpy as np

def RK4(x0, dx, y0, SlopeFn):
    m1 = SlopeFn(x0, y0)
    m2 = SlopeFn(x0+(dx/2), y0+(m1*dx/2))
    m3 = SlopeFn(x0+(dx/2), y0+(m2*dx/2))
    m4 = SlopeFn(x0+dx, y0+(m3*dx))
    
    change = (dx/6)*(m1 + 2*m2 + 2*m3 + m4)
    return y0 + change

# define functions
def getCf(H, THETA):
    ReTheta = rho*THETA/mu;
    Cf_var = 0.246*np.power(10, -0.678*H)*np.power(ReTheta, -0.268);
    return Cf_var

def getH(H1_var):
    if H1_var<=3.3:
        H_var = 3;
    elif H1_var>3.3 and H1_var<=5.3:
        H_var = 0.6778 + 1.1536*(H1_var - 3.3) - 0.326;
    elif H1_var >5.3:
        H_var = 1.1 + 0.86*(H1_var - 3.3) -0.777;
    else:
        print("Something is wrong, check if H is less than zero")
    return H_var


# initial values #####################################################################################################
n        = 41; # 41 is the upper limit, any further grid points cause floating point errors
rho      = 1.125;
mu       = 1.6*10**-5
######################################################################################################################



# Assume a grid  for the airfoil (x is the number of panels + 1) #####################################################
x        = np.linspace(0, 1, n);
dx       = x[1] - x[0];
######################################################################################################################




# Assume a pressure gradient over the airfoil ########################################################################
P        = 5*np.power(x,2) + 2*x + 1;
######################################################################################################################




# find external velocities ove the airfoil ###########################################################################
v        = np.zeros(n)
v[0]     = 100; # initial value of the velocity
delPterm = np.zeros(n)
for i in range(1,n-1):
    delPterm[i] = 2*(P[0] - P[i])/rho;
    v[i]        = np.sqrt(delPterm[i] + np.power(v[0],2));
######################################################################################################################
  
  
    

# find dVe/dx everywhere using backward difference over the airfoil ##################################################
dv = np.zeros(n-1);
for i in range(1,n-1):
    dv[i]       = (v[i] - v[i-1])/dx;
######################################################################################################################
    


# start importing from here
# initialise for RK4 #################################################################################################
H        = np.zeros(n); H[0] = 1.4;
H1       = np.zeros(n); H1[0] = 3.47485;
Theta    = np.zeros(n); Theta[0] = 0.0003474;
Cf       = np.zeros(n); Cf[0] = 0.0005066;
######################################################################################################################




# RK4 solver #########################################################################################################
for i in range(1,n-1):
    # initialise the variable values at the current grid point using previous grid points (to define the error functions)
    H_er = H[i-1];  Cf_er = Cf[i-1];  H1_er = H[i-1];  Theta_er = Theta[i-1];
    # assign previous grid point values of H and Cf to start RK4
    H[i] = H[i-1]; Cf[i] = Cf[i-1];
    
    #assume some error values
    erH = 0.2; erH1 = 0.2; erTheta = 0.2; erCf = 0.2;
    # iterate to get the variables at the grid point
    while abs(erH)>0.00001 or abs(erH1)>0.00001 or abs(erTheta)>0.00001 or abs(erCf)>0.00001:
        # define RK4 function for Theta
        def dTheta_by_dx(x, THETA):
            return 0.5*Cf[i] - (THETA/v[i])*(2+H[i])*(dv[i])
        # get theta
        Theta[i] = RK4(x[i-1], dx, Theta[i-1], dTheta_by_dx)
        # define RK4 function for H1
        def dH1_by_dx(x, H1_var):
            term1 = 0.0306*(H1_var-3)**(-0.6169)
            term2 = -Theta[i]*H1_var*dv[i]/v[i]
            term3 = -H1_var*Cf[i]/2
            term4 = (H1_var*Theta[i]/v[i])*(2+H[i])*dv[i]
            return (term1 + term2 + term3 +term4)*(1/Theta[i])
        # get H1
        H1[i] = RK4(x[i-1], dx, H1[i-1], dH1_by_dx)
        # get H
        H[i] = getH(H1[i])
        # get skin friction
        Cf[i] = getCf(H[i], Theta[i])
        # define errors
        erH = (H[i]-H_er)/H[i];
        erH1 = (H1[i]-H1_er)/H1[i];
        erTheta = (Theta[i]-Theta_er)/Theta[i];
        erCf = (Cf[i]-Cf_er)/Cf[i];
        # assign current iteration variable values to the Var_er
        H_er = H[i]; H1_er = H1[i]; Theta_er = Theta[i]; Cf_er = Cf[i];
######################################################################################################################





    
    
    
    
    
    
    
    
    
    