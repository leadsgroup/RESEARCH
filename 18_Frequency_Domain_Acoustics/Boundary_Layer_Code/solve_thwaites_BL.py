
import numpy as np  
from scipy.interpolate import interp1d
from scipy.integrate import odeint 
def solve_thwaites_BL(l,Re_L,x_i,Ve_i,dVe_i, theta_0):  
    ''' Solved intergral boundary layer equations for laminar flow  
    '''    
    nu          = l/Re_L     
    y0          = theta_0**2 * getVe(0,x_i,Ve_i)**6   
    theta2_Ve6  = odeint(odefcn, y0,x_i , args=(nu, x_i, Ve_i))  # eqn 25 , Solving Intergral Boundary Layer Lecture 
    
    # Compute momentum thickness, theta 
    theta       = np.sqrt(theta2_Ve6[:,0]/ Ve_i**6) # eqn  25 , Solving Intergral Boundary Layer Lecture 
    
    # Thwaites separation criteria 
    lambda_val  = theta**2 * dVe_i / nu  # 21 , Solving Intergral Boundary Layer Lecture 
    
    # Compute H 
    H           = getH(lambda_val) 
    
    # Compute Reynolds numbers based on momentum thickness  
    Re_theta    = Ve_i * theta / nu
    
    # Compute Reynolds numbers based on distance along airfoil
    Re_x        = Ve_i * x_i/ nu
    
    # Compute skin friction 
    cf          = abs(getcf(lambda_val ,Re_theta)) 
    
    # Compute displacement thickness
    delta_star    = H*theta   
    
    # Compute boundary layer thickness 
    delta       = 5.2*x_i/np.sqrt(Re_x)   # eqn 5 , Solving Intergral Boundary Layer Lecture 
     
    return H,delta_star,delta,cf,theta,Re_x,Re_theta


def getH(lambda_val ): 
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
    
def odefcn(y,x, nu,x_i,Ve_i):
    """ Computes boundary layer functions using SciPy ODE solver 

    Assumptions:
    None

    Source:
    None

    Inputs: 
    y           - initial conditions of functions    [unitless]
    x           - new x values at which to solve ODE [unitless]
    nu          - kinematic viscosity                [m^2/s]
    x_i         - intial array of x values           [unitless]
    Ve_i        - intial boundary layer velocity     [m/s]
    
    Outputs:  
    dydx        - expression for the momentum thickness and velocity (theta**2/Ve**6)

    Properties Used:
    N/A 
    """        
    dydx = 0.45*getVe(x,x_i,Ve_i)**5*nu
    return dydx 
    
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
    Ve_func = interp1d(x_i,Ve_i, axis=0,fill_value = "extrapolate")
    Ve      = Ve_func(x)
    return Ve  

def getdVe(x,x_i,dVe_i):
    """ Interpolates the derivatives of the bounday layer velocity over a new dimension of x 

    Assumptions:
    None

    Source:
    None

    Inputs: 
    x         - new x dimension                                   [unitless]
    x_i       - old x dimension                                   [unitless]
    dVe_i     - old derivative of boundary layer velocity values  [m/s-m]
    
    Outputs:  
    dVe       - new derivative of boundary layer velocity values  [m/s-m]

    Properties Used:
    N/A 
    """
    dVe_func = interp1d(x_i,dVe_i,fill_value = "extrapolate")
    dVe      = dVe_func(x)
    return dVe 

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