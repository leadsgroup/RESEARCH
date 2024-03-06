
import numpy as np    
import matplotlib.pyplot as plt 
from  solve_thwaites_BL import solve_thwaites_BL
from  solve_heads_BL_prev_index    import solve_heads_BL_prev_index
from  michel_criterion  import michel_criterion
from solve_heads_BL_curr_index  import solve_heads_BL_curr_index
from solve_heads_BL   import solve_heads_BL 

def main():
     
    l          = 1             # characteristic length of 1
    Re_L       = 10**(7)
    n_points   = 300          # Discretizing domain 
    x          = np.linspace(0,l,n_points)
    eps        = 10**(-12)     # First node location
    x[0]       = eps  
    theta_0    = 1E-12 # initial momentum thickness 
    Ve      = np.ones(n_points)       
    dVe     = np.gradient(Ve,x)
    nu = Ve[0]*l/Re_L
    
    # Thwaites method for laminar flow 
    H_lam,delta_star_lam,delta_lam,cf_lam,theta_lam,Re_x_lam,Re_theta_lam = solve_thwaites_BL(l,Re_L,x,Ve,dVe,theta_0) 
 
    # Find transition point 
    transition_index,transition_flag = michel_criterion(Re_theta_lam,Re_x_lam)
    
    if transition_flag:
        # Get initial conditons for turbulent boundary layer 
        delta_0_tur       = delta_lam[transition_index]
        delta_star_0_tur  = delta_star_lam[transition_index]
        theta_0_tur       = theta_lam[transition_index]
        cf_0_tur          = cf_lam[transition_index]
        H_0_tur           = H_lam[transition_index]
        Re_L_tur          = Re_L 
        x_tur             = x[transition_index:]
        Ve_tur            = Ve[transition_index:]
        dVe_tur           = dVe[transition_index:]
        
        # Turbulent Boundary Layer Calculation
        # old is using previous index
        # new is using current index
        # ode is using odeint function
        H_turb_old,delta_star_turb_old,delta_turb_old,cf_turb_old,theta_turb_old,Re_x_turb_old,Re_theta_turb_old = solve_heads_BL_prev_index(nu, l,delta_0_tur,theta_0_tur,delta_star_0_tur,cf_0_tur,H_0_tur,Re_L_tur,x_tur,Ve_tur,dVe_tur)
        H_turb_new,delta_star_turb_new,delta_turb_new,cf_turb_new,theta_turb_new,Re_x_turb_new,Re_theta_turb_new = solve_heads_BL_curr_index(nu, l,delta_0_tur,theta_0_tur,delta_star_0_tur,cf_0_tur,H_0_tur,Re_L_tur,x_tur,Ve_tur,dVe_tur)
        H_turb_ode,delta_star_turb_ode,delta_turb_ode,cf_turb_ode,theta_turb_ode,Re_x_turb_ode,Re_theta_turb_ode = solve_heads_BL(nu, l,delta_0_tur,theta_0_tur,delta_star_0_tur,Re_L_tur,x_tur,Ve_tur,dVe_tur)
        
        # Concatenate vectors 
        delta_star_old  = np.hstack((delta_star_lam[:transition_index],delta_star_turb_old))  
        delta_old       = np.hstack((delta_lam[:transition_index],delta_turb_old))
        cf_old          = np.hstack((cf_lam[:transition_index],cf_turb_old))
        theta_old       = np.hstack((theta_lam[:transition_index],theta_turb_old))
        Re_x_old        = np.hstack((Re_x_lam[:transition_index],Re_x_turb_old))
        Re_theta_old    = np.hstack((Re_theta_lam[:transition_index],Re_theta_turb_old) )  
        H_old           = np.hstack((H_lam[:transition_index],H_turb_old))
         

        delta_star_new  = np.hstack((delta_star_lam[:transition_index],delta_star_turb_new))  
        delta_new       = np.hstack((delta_lam[:transition_index],delta_turb_new))
        cf_new          = np.hstack((cf_lam[:transition_index],cf_turb_new))
        theta_new       = np.hstack((theta_lam[:transition_index],theta_turb_new))
        Re_x_new        = np.hstack((Re_x_lam[:transition_index],Re_x_turb_new))
        Re_theta_new    = np.hstack((Re_theta_lam[:transition_index],Re_theta_turb_new) )  
        H_new           = np.hstack((H_lam[:transition_index],H_turb_new))
    
        delta_star_ode  = np.hstack((delta_star_lam[:transition_index],delta_star_turb_ode))  
        delta_ode       = np.hstack((delta_lam[:transition_index],delta_turb_ode))
        cf_ode          = np.hstack((cf_lam[:transition_index],cf_turb_ode))
        theta_ode       = np.hstack((theta_lam[:transition_index],theta_turb_ode))
        Re_x_ode        = np.hstack((Re_x_lam[:transition_index],Re_x_turb_ode))
        Re_theta_ode    = np.hstack((Re_theta_lam[:transition_index],Re_theta_turb_ode) )  
        H_ode           = np.hstack((H_lam[:transition_index],H_turb_ode))
        
        
    else:
        delta_star = delta_star_lam
        delta      = delta_lam
        cf         = cf_lam
        theta      = theta_lam 
        Re_x       = Re_x_lam 
        Re_theta   = Re_theta_lam  
        H          = H_lam
        
    
    #print('Transition Point: ', x[transition_index])
    #Cf = np.sum(cf[1:]*np.diff(x))
    #print(Cf)
    # Plots 
    
    # old is blue (previous index) and new is red (current index)
    fig  = plt.figure('Boundary Layer Properties') 
    fig .set_size_inches(12,7)   
    axis1  = fig.add_subplot(3,2,1)  
    axis2  = fig.add_subplot(3,2,2) 
    axis3  = fig.add_subplot(3,2,3) 
    axis4  = fig.add_subplot(3,2,4) 
    axis5  = fig.add_subplot(3,2,5) 
    axis6  = fig.add_subplot(3,2,6)              
    axis1.set_ylabel(r'$\delta$')           
    axis2.set_ylabel(r'$\delta$*')         
    axis3.set_ylabel(r'$\theta$')         
    axis4.set_ylabel(r'$c_f$')         
    axis5.set_ylabel(r'$Re_x$')         
    axis6.set_ylabel(r'H')   
    axis5.set_xlabel('x (m)')    
    axis6.set_xlabel('x (m)')
    axis4.set_ylim([0,0.01])
    axis1.plot(x,delta_old, 'b-')  
    axis2.plot(x,delta_star_old, 'b-')  
    axis3.plot(x,theta_old, 'b-')  
    axis4.plot(x,cf_old, 'b-')  
    axis5.plot(x,Re_x_old, 'b-')  
    axis6.plot(x,H_old, 'b-') 

    axis1.plot(x,delta_new, 'r-')  
    axis2.plot(x,delta_star_new, 'r-')  
    axis3.plot(x,theta_new, 'r-')  
    axis4.plot(x,cf_new, 'r-')  
    axis5.plot(x,Re_x_new, 'r-')  
    axis6.plot(x,H_new, 'r-') 
    
    axis1.plot(x,delta_ode, 'g-')  
    axis2.plot(x,delta_star_ode, 'g-')  
    axis3.plot(x,theta_ode, 'g-')  
    axis4.plot(x,cf_ode, 'g-')  
    axis5.plot(x,Re_x_ode, 'g-')  
    axis6.plot(x,H_ode, 'g-') 
    
    set_axes(axis1)
    set_axes(axis2)
    set_axes(axis3)
    set_axes(axis4)
    set_axes(axis5)
    set_axes(axis6) 
     
    return 
 
## @ingroup Visualization-Performance-Common
def set_axes(axis):
    """This sets the axis parameters for all plots

    Assumptions:
    None

    Source:
    None

    Inputs
    axes

    Outputs:
    axes

    Properties Used:
    N/A
    """

    axis.minorticks_on()
    axis.grid(which='major', linestyle='-', linewidth=0.5, color='grey')
    axis.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')
    axis.grid(True)
    axis.get_yaxis().get_major_formatter().set_scientific(False)
    axis.get_yaxis().get_major_formatter().set_useOffset(False)

    return

if __name__ == '__main__':
    main()
    plt.show()
 