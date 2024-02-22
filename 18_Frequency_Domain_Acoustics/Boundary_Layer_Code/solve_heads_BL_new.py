 
import numpy as np 
# ----------------------------------------------------------------------
# heads_method.py 
# ----------------------------------------------------------------------    
def solve_heads_BL_new(l,del_0,theta_0,del_star_0,cf_0,H_0_tur,Re_L,x_i,Ve_i,dVe_i):    
    ''' Solved intergral boundary layer equations for turbulent flow  
    '''
    nu       = l/Re_L    
    n        = len(x_i)    # 41 is the upper limit, any further grid points cause floating point errors
    rho      = 1.125
    mu       = nu*rho 
 
    x        = x_i
    dx       = np.diff(x)
    v        = Ve_i 
    dv       = dVe_i

    # start importing from here
    # initialise for RK4 #################################################################################################
    H        = np.zeros(n) 
    H[0]     = 1.4
    Theta    = np.zeros(n)
    Theta[0] = theta_0 
    H1       = np.zeros(n) 
    H1[0]    = getH1(np.array([H[0]])) # 3.47485 
    Cf       = np.zeros(n) 
    Cf[0]    = cf_0 

    # RK4 solver #########################################################################################################
    for i in range(1,n):
        # initialise the variable values at the current grid point using previous grid points (to define the error functions)
        H_er = H[i-1]  
        Cf_er = Cf[i-1] 
        H1_er = H[i-1]  
        Theta_er = Theta[i-1]
        
        # assign previous grid point values of H and Cf to start RK4
        H[i] = H[i-1]
        Cf[i] = Cf[i-1]

        #assume some error values
        erH = 0.2 
        erH1 = 0.2 
        erTheta = 0.2
        erCf = 0.2
        # iterate to get the variables at the grid point
        while abs(erH)>0.00001 or abs(erH1)>0.00001 or abs(erTheta)>0.00001 or abs(erCf)>0.00001:
            # get theta
            Theta[i] = RK4(x[i-1], dx[i-1], Theta[i-1], dTheta_by_dx(x, Theta[i-1],args = (Cf, v, i, H, dv))) # ERROR IN LINE 
            
            # define RK4 function for H1
            # get H1
            H1[i]   = RK4(x[i-1], dx[i-1], H1[i-1], dH1_by_dx(x, H1[i-1], args = (dv, v, H, Cf, Theta, i))) # ERROR IN LINE 
            
            # get H
            H[i]    = getH(H1[i])
            
            # get skin friction
            Cf[i]   = getCf(H[i], Theta[i],rho,mu)
            
            # define errors
            erH     = (H[i]-H_er)/H[i]
            erH1    = (H1[i]-H1_er)/H1[i]
            erTheta = (Theta[i]-Theta_er)/Theta[i]
            erCf    = (Cf[i]-Cf_er)/Cf[i]
            
            # assign current iteration variable values to the Var_er
            H_er     = H[i] 
            H1_er    = H1[i] 
            Theta_er = Theta[i]
            Cf_er    = Cf[i]
              
    
    H1           = getH1(H)
    theta        = Theta 
    delta_star   = H*theta    
    delta        = theta*H1 + delta_star  
    Re_theta     = Re_L/l * Ve_i*theta  
    Re_x         = Ve_i* x_i / nu  
    cf           =  Cf
    
    return H,delta_star,delta,cf,theta,Re_x,Re_theta  


# define RK4 function for Theta
def dTheta_by_dx(x, THETA,Cf,v,i,H,dv):
    return 0.5*Cf[i] - (THETA/v[i])*(2+H[i])*(dv[i])

def dH1_by_dx(x, H1_var,dv,v,H,Cf,Theta,i):
    term1 = 0.0306*(H1_var-3)**(-0.6169)
    term2 = -Theta[i]*H1_var*dv[i]/v[i]
    term3 = -H1_var*Cf[i]/2
    term4 = (H1_var*Theta[i]/v[i])*(2+H[i])*dv[i]
    return (term1 + term2 + term3 +term4)*(1/Theta[i])

def RK4(x0, dx, y0, SlopeFn):
    m1 = SlopeFn(x0, y0)
    m2 = SlopeFn(x0+(dx/2), y0+(m1*dx/2))
    m3 = SlopeFn(x0+(dx/2), y0+(m2*dx/2))
    m4 = SlopeFn(x0+dx, y0+(m3*dx))

    change = (dx/6)*(m1 + 2*m2 + 2*m3 + m4)
    return y0 + change

def getH1(H) :   
    H1          = 3.3 + 0.8234*((H - 1.1)**(-1.287))  
    H1[H > 1.6] = 3.3 + 1.5501*((H[H > 1.6] - 0.6778)**-3.064) 
    return H1 
    

# define functions
def getCf(H, THETA,rho,mu):
    ReTheta = rho*THETA/mu
    Cf_var = 0.246* (10**(-0.678*H)) * (ReTheta**-0.268)
    return Cf_var

def getH(H1_var):
    if H1_var<=3.3:
        H_var = 3
    elif H1_var>3.3 and H1_var<=5.3:
        H_var = 0.6778 + 1.1536*(H1_var - 3.3) - 0.326
    elif H1_var >5.3:
        H_var = 1.1 + 0.86*(H1_var - 3.3) -0.777
    else:
        print("Something is wrong, check if H is less than zero")
    return H_var

 