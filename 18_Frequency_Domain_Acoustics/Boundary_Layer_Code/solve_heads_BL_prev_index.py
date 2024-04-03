 
import numpy as np 
# ----------------------------------------------------------------------
# heads_method.py 
# ----------------------------------------------------------------------    
def solve_heads_BL_prev_index(nu, l, delta_0, theta_0, delta_star_0, cf_0, H_0_tur, Re_L, x_i, Ve_i, dVe_i):    
    ''' Solved intergral boundary layer equations for turbulent flow  
    '''
    # define functions
    def getCf(H, THETA):
        ReTheta = rho*THETA/mu;
        Cf_var = 0.246*(10**(-0.678*H))*(ReTheta**-0.268);
        return Cf_var

    def getH(H1_var):
        if H1_var<=0:
            print("Something is wrong, check if H is less than zero")        
        elif H1_var < 5.39142:
            H_var = 0.6778 + 1.153793*(H1_var-3.3)**-0.32637;
        elif H1_var >= 5.39142:
            H_var = 1.1 + 0.8598636*(H1_var - 3.3)**-0.777;
        return H_var
    
    # define RK4 function for Theta
    def dTheta_by_dx(index, x, THETA):
        return 0.5*Cf[index] - (THETA/v[index])*(2+H[index])*(dv[index])
    
    def dVeThetaH1_by_dx(index, x, VeThetaH1_var):
        return v[index]*0.0306*(((VeThetaH1_var/(v[index]*Theta[index]))-3)**-0.6169)
    # define RK4 function for H1
    # def dH1_by_dx(index, x, H1_var):
    #     term1 = 0.0306*(H1_var-3)**(-0.6169)
    #     term2 = -Theta[index]*H1_var*dv[index]/v[index]
    #     term3 = -H1_var*Cf[index]/2
    #     term4 = (H1_var*Theta[index]/v[index])*(2+H[index])*dv[index]
    #     return (term1 + term2 + term3 +term4)*(1/Theta[index])
 
    n        = len(x_i)
    rho      = 1.125
    mu       = nu*rho 
    x        = x_i
    dx       = np.diff(x)
    v        = Ve_i 
    dv       = dVe_i

    # start importing from here
    # initialise for RK4 #################################################################################################
    H            = np.zeros(n) 
    H[0]         = H_0_tur
    Theta        = np.zeros(n)
    Theta[0]     = theta_0 
    H1           = np.zeros(n) 
    H1[0]        = (delta_0 - delta_star_0)/theta_0 
    Cf           = np.zeros(n) 
    Cf[0]        = cf_0 
    VeThetaH1    = np.zeros(n)
    VeThetaH1[0] = v[0]*Theta[0]*H1[0]

    # RK4 solver #########################################################################################################
    for i in range(1,n):
        # initialise the variable values at the current grid point using previous grid points (to define the error functions)
        H_er = H[i-1];  Cf_er = Cf[i-1];  H1_er = H[i-1];  Theta_er = Theta[i-1];
        # assign previous grid point values of H and Cf to start RK4
        H[i] = H[i-1]; Cf[i] = Cf[i-1];
        
        #assume some error values
        erH = 0.2; erH1 = 0.2; erTheta = 0.2; erCf = 0.2;
        
        # iterate to get the variables at the grid point
        while abs(erH)>0.00001 or abs(erH1)>0.00001 or abs(erTheta)>0.00001 or abs(erCf)>0.00001:
            
            # get theta
            Theta[i] = RK4(i-1, x, dx, Theta, dTheta_by_dx)
            
            # get VeThetaH1
            VeThetaH1[i] = RK4(i-1, x, dx, VeThetaH1, dVeThetaH1_by_dx)
            
            # get H1
            H1[i] = VeThetaH1[i]/(v[i]*Theta[i])
            
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
            H_er = H[i]
            H1_er = H1[i]
            Theta_er = Theta[i]
            Cf_er = Cf[i]
              
    
    delta_star   = H*Theta    
    delta        = (Theta*H1) + delta_star  
    Re_theta     = Re_L/l * Ve_i*Theta  
    Re_x         = Ve_i* x_i / nu  
    cf           = Cf
    
    return H,delta_star,delta,cf,Theta,Re_x,Re_theta  
 
def RK4(ind, x, dx, y, SlopeFn):
    # ind is the index
    m1 = SlopeFn(ind, x[ind], y[ind])
    m2 = SlopeFn(ind, x[ind]+(dx[ind]/2), y[ind]+(m1*dx[ind]/2))
    m3 = SlopeFn(ind, x[ind]+(dx[ind]/2), y[ind]+(m2*dx[ind]/2))
    m4 = SlopeFn(ind, x[ind]+dx[ind], y[ind]+(m3*dx[ind]))
    
    change = (dx[ind]/6)*(m1 + 2*m2 + 2*m3 + m4)
    return y[ind] + change

# import numpy as np
# from scipy.interpolate import interp1d
# from scipy.integrate import odeint  
# # ----------------------------------------------------------------------
# # heads_method.py 
# # ----------------------------------------------------------------------    
# def solve_heads_BL(nu,l,del_0,theta_0,del_star_0,Re_L,x_i,Ve_i,dVe_i):    
#     ''' Solved intergral boundary layer equations for turbulent flow  
#     '''
#     #nu           = l/Re_L       
#     H_0          = del_star_0 / theta_0                        # page 126 
#     H1_0         = getH1(np.atleast_1d(H_0))[0]                # page 127 
#     y0           = [theta_0, Ve_i[0]*theta_0*H1_0]   # page 128   
#     y            = odeint(odefcn,y0,x_i,args=(Re_L/l, x_i, Ve_i, dVe_i))  
    
#     # Compute momentum thickness, theta 
#     theta         = y[:,0]    
#     Ve_theta_H1   = y[:,1]  
      
#     # Compute mass flow shape factor, H1
#     H1           = Ve_theta_H1/(theta*Ve_i)
    
#     # Compute H 
#     H            = getH(np.atleast_1d(H1))    
    
#     # Compute Reynolds numbers based on momentum thickness  
#     Re_theta     = Re_L/l * Ve_i*theta 
    
#     # Compute Reynolds numbers based on distance along airfoil
#     Re_x         = Ve_i* x_i / nu
    
#     # Compute skin friction 
#     cf           = abs( getcf(np.atleast_1d(Re_theta),np.atleast_1d(H))) 
    
#     # Compute displacement thickness
#     delta_star   = H*theta   
    
#     # Compute boundary layer thickness 
#     delta        = theta*H1 + delta_star  
                
#     return H,delta_star,delta,cf,theta,Re_x,Re_theta  

# def getH(H1):
#     """ Computes the shape factor, H
#     Assumptions:
#     None
#     Source:
#     None
#     Inputs: 
#     H1       - mass flow shape factor [unitless]
#     Outputs:  
#     H        - shape factor [unitless]
#     Properties Used:
#     N/A
#     """         
#     H       = 0.6778 + 1.1536*(H1-3.3)**-0.326
#     idx1    = (H1 < 3.3)
#     H[idx1] = 3.0
#     idx2    = (H1 > 5.3)
#     H[idx2] = 1.1 + 0.86*(H1[idx2] - 3.3)**-0.777 
#     return H 

# def getH1(H) :    
#     """ Computes the mass flow shape factor, H1
#     Assumptions:
#     None
#     Source:
#     None
#     Inputs: 
#     H        - shape factor [unitless]
#     Outputs:  
#     H1       - mass flow shape factor [unitless]
#     Properties Used:
#     N/A 
#     """
#     H1       = 3.3 + 0.8234*(H - 1.1)**-1.287  
#     idx1     = (H > 1.6) 
#     H1[idx1] = 3.3 + 1.5501*(H[idx1] - 0.6778)**-3.064
#     return H1 

# def odefcn(y,x,ReL_div_L, x_i, Ve_i, dVe_i): 
#     """ Computes boundary layer functions using SciPy ODE solver 
#     Assumptions:
#     None
#     Source:
#     None
#     Inputs:  
#     y           - initial conditions of functions               [unitless]
#     x           - new x values at which to solve ODE            [unitless]
#     ReL_div_L   - ratio of Reynolds number to length of surface [unitless]
#     x_i         - intial array of x values                      [unitless]
#     Ve_i        - intial boundary layer velocity                [m/s]
#     dVe_i       - initial derivative of bounday layer velocity  [m/s-m]

#     Outputs:  
#     f           - 2D function of momentum thickness and the product of 
#                   the velocity,momentum thickness and the mass flow shape factor
#     Properties Used:
#     N/A 
#     """    
#     theta       = y[0]
#     Ve_theta_H1 = y[1]  

#     if theta == 0:
#         H1 = Ve_theta_H1 / (theta + 1e-6) / getVe(x,x_i,Ve_i)
#     else:
#         H1 = Ve_theta_H1 / theta / getVe(x,x_i,Ve_i)

#     H           = getH(np.atleast_1d(H1))
#     Re_theta    = ReL_div_L * theta
#     cf          = getcf(np.atleast_1d(Re_theta),np.atleast_1d(H))
#     dydx_1      = 0.5*cf-(theta/getVe(x,x_i,Ve_i))*(2+H)*getdVe(x, x_i, dVe_i)
#     dydx_2      = getVe(x,x_i,Ve_i)*0.0306*(H1 - 3)**-0.6169 
#     f           = [dydx_1,dydx_2] 
#     return f 

# def getVe(x,x_i,Ve_i):
#     """ Interpolates the bounday layer velocity over a new dimension of x 
#     Assumptions:
#     None
#     Source:
#     None
#     Inputs: 
#     x         - new x dimension                    [unitless]
#     x_i       - old x dimension                    [unitless]
#     Ve_i      - old boundary layer velocity values [m/s] 

#     Outputs:  
#     Ve        - new boundary layer velocity values [m/s]
#     Properties Used:
#     N/A 
#     """    
#     Ve_func = interp1d(x_i,Ve_i,fill_value = "extrapolate")
#     Ve      = Ve_func(x)
#     return Ve 

# def getdVe(x,x_i,dVe_i):
#     """ Interpolates the derivatives of the bounday layer velocity over a new dimension of x

#     Assumptions:
#     None
#     Source:
#     None
#     Inputs: 
#     x         - new x dimension                                  [unitless]
#     x_i       - old x dimension                                  [unitless]
#     dVe_i     - old derivative of boundary layer velocity values [m/s-m] 

#     Outputs:  
#     dVe       - new derivative of boundary layer velocity values [m/s-m]
#     Properties Used:
#     N/A 
#     """        
#     dVe_func = interp1d(x_i,dVe_i,fill_value = "extrapolate")
#     dVe      = dVe_func(x)
#     return dVe  

# def getcf(Re_theta,H): 
#     """ Computes the skin friction coefficient, cf
#     Assumptions:
#     None
#     Source:
#     None
#     Inputs: 
#     Re_theta - Reynolds Number as a function of momentum thickness [m]
#     H        - shape factor                                        [unitless]
#     Outputs:  
#     cf       - skin friction coefficient  [unitless]
#     Properties Used:
#     N/A 
#     """    
#     cf       = 0.246*10**(-0.678*H)*(Re_theta)**-0.268 
#     idx1     = (Re_theta == 0) 
#     cf[idx1] = 0.246*10**(-0.678*H[idx1])*(1e-3)**-0.268 
#     return cf 