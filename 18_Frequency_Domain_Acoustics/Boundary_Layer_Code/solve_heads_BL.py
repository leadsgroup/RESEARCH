 
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import odeint  
# ----------------------------------------------------------------------
# heads_method.py 
# ----------------------------------------------------------------------    
def solve_heads_BL(l,del_0,theta_0,del_star_0,Re_L,x_i,Ve_i,dVe_i):    
    ''' Solved intergral boundary layer equations for turbulent flow  
    '''
    nu           = l/Re_L       
    H_0          = del_star_0 / theta_0                        # page 126 
    H1_0         = getH1(np.atleast_1d(H_0))[0]                # page 127 
    y0           = [theta_0, Ve_i[0]*theta_0*H1_0]   # page 128   
    y            = odeint(odefcn,y0,x_i,args=(Re_L/l, x_i, Ve_i, dVe_i))  
    
    # Compute momentum thickness, theta 
    theta         = y[:,0]    
    Ve_theta_H1   = y[:,1]  
      
    # Compute mass flow shape factor, H1
    H1           = Ve_theta_H1/(theta*Ve_i)
    
    # Compute H 
    H            = getH(np.atleast_1d(H1))    
    
    # Compute Reynolds numbers based on momentum thickness  
    Re_theta     = Re_L/l * Ve_i*theta 
    
    # Compute Reynolds numbers based on distance along airfoil
    Re_x         = Ve_i* x_i / nu
    
    # Compute skin friction 
    cf           = abs( getcf(np.atleast_1d(Re_theta),np.atleast_1d(H))) 
    
    # Compute displacement thickness
    delta_star   = H*theta   
    
    # Compute boundary layer thickness 
    delta        = theta*H1 + delta_star  
                
    return H,delta_star,delta,cf,theta,Re_x,Re_theta  

def getH(H1):
    """ Computes the shape factor, H
    Assumptions:
    None
    Source:
    None
    Inputs: 
    H1       - mass flow shape factor [unitless]
    Outputs:  
    H        - shape factor [unitless]
    Properties Used:
    N/A
    """         
    H       = 0.6778 + 1.1536*(H1-3.3)**-0.326
    idx1    = (H1 < 3.3)
    H[idx1] = 3.0
    idx2    = (H1 > 5.3)
    H[idx2] = 1.1 + 0.86*(H1[idx2] - 3.3)**-0.777 
    return H 

def getH1(H) :    
    """ Computes the mass flow shape factor, H1
    Assumptions:
    None
    Source:
    None
    Inputs: 
    H        - shape factor [unitless]
    Outputs:  
    H1       - mass flow shape factor [unitless]
    Properties Used:
    N/A 
    """
    H1       = 3.3 + 0.8234*(H - 1.1)**-1.287  
    idx1     = (H > 1.6) 
    H1[idx1] = 3.3 + 1.5501*(H[idx1] - 0.6778)**-3.064
    return H1 

def odefcn(y,x,ReL_div_L, x_i, Ve_i, dVe_i): 
    """ Computes boundary layer functions using SciPy ODE solver 
    Assumptions:
    None
    Source:
    None
    Inputs:  
    y           - initial conditions of functions               [unitless]
    x           - new x values at which to solve ODE            [unitless]
    ReL_div_L   - ratio of Reynolds number to length of surface [unitless]
    x_i         - intial array of x values                      [unitless]
    Ve_i        - intial boundary layer velocity                [m/s]
    dVe_i       - initial derivative of bounday layer velocity  [m/s-m]

    Outputs:  
    f           - 2D function of momentum thickness and the product of 
                  the velocity,momentum thickness and the mass flow shape factor
    Properties Used:
    N/A 
    """    
    theta       = y[0]
    Ve_theta_H1 = y[1]  

    if theta == 0:
        H1 = Ve_theta_H1 / (theta + 1e-6) / getVe(x,x_i,Ve_i)
    else:
        H1 = Ve_theta_H1 / theta / getVe(x,x_i,Ve_i)

    H           = getH(np.atleast_1d(H1))
    Re_theta    = ReL_div_L * theta
    cf          = getcf(np.atleast_1d(Re_theta),np.atleast_1d(H))
    dydx_1      = 0.5*cf-(theta/getVe(x,x_i,Ve_i))*(2+H)*getdVe(x, x_i, dVe_i)
    dydx_2      = getVe(x,x_i,Ve_i)*0.0306*(H1 - 3)**-0.6169 
    f           = [dydx_1,dydx_2] 
    return f 

def getVe(x,x_i,Ve_i):
    """ Interpolates the bounday layer velocity over a new dimension of x 
    Assumptions:
    None
    Source:
    None
    Inputs: 
    x         - new x dimension                    [unitless]
    x_i       - old x dimension                    [unitless]
    Ve_i      - old boundary layer velocity values [m/s] 

    Outputs:  
    Ve        - new boundary layer velocity values [m/s]
    Properties Used:
    N/A 
    """    
    Ve_func = interp1d(x_i,Ve_i,fill_value = "extrapolate")
    Ve      = Ve_func(x)
    return Ve 

def getdVe(x,x_i,dVe_i):
    """ Interpolates the derivatives of the bounday layer velocity over a new dimension of x

    Assumptions:
    None
    Source:
    None
    Inputs: 
    x         - new x dimension                                  [unitless]
    x_i       - old x dimension                                  [unitless]
    dVe_i     - old derivative of boundary layer velocity values [m/s-m] 

    Outputs:  
    dVe       - new derivative of boundary layer velocity values [m/s-m]
    Properties Used:
    N/A 
    """        
    dVe_func = interp1d(x_i,dVe_i,fill_value = "extrapolate")
    dVe      = dVe_func(x)
    return dVe  

def getcf(Re_theta,H): 
    """ Computes the skin friction coefficient, cf
    Assumptions:
    None
    Source:
    None
    Inputs: 
    Re_theta - Reynolds Number as a function of momentum thickness [m]
    H        - shape factor                                        [unitless]
    Outputs:  
    cf       - skin friction coefficient  [unitless]
    Properties Used:
    N/A 
    """    
    cf       = 0.246*10**(-0.678*H)*(Re_theta)**-0.268 
    idx1     = (Re_theta == 0) 
    cf[idx1] = 0.246*10**(-0.678*H[idx1])*(1e-3)**-0.268 
    return cf 