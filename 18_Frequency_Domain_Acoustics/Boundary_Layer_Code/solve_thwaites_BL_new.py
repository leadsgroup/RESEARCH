
"""
Created on Wed Apr 17 17:21:40 2024

@author: nn18
"""

import numpy as np

def solve_thwaites_BL_new(nu,l,Re_L,x_i,Ve_i,dVe_i,theta_0):
    
    def dy_by_dx(index, X, Y):
        return 0.45*nu*Ve_i[index]**5    
    
    n = len(x_i)
    rho      = 1.125
    mu       = nu*rho 
    dx_i       = np.diff(x_i)
    
    Theta        = np.zeros(n)
    Theta[0]     = theta_0
    y            = np.zeros(n)
    y[0]         = (theta_0**2)*(Ve_i[0])**6
    
    # determine (Theta**2)*(Ve**6)
    for i in range(1,n):
        y[i] = RK4(i-1, dx_i, x_i, y, dy_by_dx)
    
    
    # determine Theta
    Theta = np.sqrt(y/(Ve_i**6))
    
    # Thwaites separation criteria 
    lambda_val  = (Theta**2)*dVe_i/nu
    
    # Compute H 
    H           = getH(lambda_val) 
    
    # Compute Reynolds numbers based on momentum thickness  
    Re_theta    = Ve_i*Theta/nu
    
    # Compute Reynolds numbers based on distance along airfoil
    Re_x        = Ve_i*x_i/nu
    
    # Compute skin friction 
    cf          = abs(getcf(lambda_val,Re_theta)) 
    
    # Compute displacement thickness
    delta_star    = H*Theta
    
    # Compute boundary layer thickness
    delta       = 5.2*x_i/np.sqrt(Re_x)
    
    return H,delta_star,delta,cf,Theta,Re_x,Re_theta
    


def getH(lambda_val): 
    """ Computes the shape factor, H

    Assumptions:
    None

    Source:
    None

    Inputs: 
    lamdda_val  - thwaites separation criteria [unitless]

    Outputs:  
    H           - shape factor [unitless]

    Properties Used:
    N/A
    """       
    H       = 0.0731/(0.14 + lambda_val ) + 2.088 
    idx1    = (lambda_val>0.0)  
    H[idx1] = 2.61 - 3.75*lambda_val[idx1]  + 5.24*lambda_val[idx1]**2   
    return H 
    

def getcf(lambda_val , Re_theta):
    """ Computes the skin friction coefficient, cf

    Assumptions:
    None

    Source:
    None

    Inputs: 
    lambda_val - thwaites separation criteria                        [unitless]
    Re_theta   - Reynolds Number as a function of momentum thickness [unitless]

    Outputs:  
    cf         - skin friction coefficient [unitless]

    Properties Used:
    N/A 
    """        
    l       = 0.22 + 1.402*lambda_val  + (0.018*lambda_val)/(0.107 + lambda_val ) 
    idx1    = (lambda_val>0.0)   
    l[idx1] = 0.22 + 1.57*lambda_val[idx1] - 1.8*lambda_val[idx1]**2 
    cf      = 2*l/Re_theta  
    return cf   
    
    
def RK4(ind, dx, x, Var1, Slope1):
    m1 = Slope1(ind,  x[ind],  Var1[ind])
    m2 = Slope1(ind,  x[ind] + dx[ind]/2,  Var1[ind] + m1*dx[ind]/2)
    m3 = Slope1(ind,  x[ind] + dx[ind]/2,  Var1[ind] + m2*dx[ind]/2)
    m4 = Slope1(ind,  x[ind] + dx[ind],  Var1[ind] + m3*dx[ind])
    
    change = (dx[ind]/6)*(m1 + 2*m2 + 2*m3 + m4)
    return Var1[ind] + change