
import numpy as np   
from  solve_thwaites_BL import solve_thwaites_BL
from  solve_heads_BL    import solve_heads_BL
from  michel_criterion  import michel_criterion

def main():
     
    l          = 1             # characteristic length of 1
    Re_L       = 10**(7) 
    n_points   = 401           # Discretizing domain 
    x          = np.linspace(0,l,n_points)
    eps        = 10**(-12) # First node location
    x[0]       = eps  
  
    Ve      = np.ones(n_points)       
    dVe     = np.gradient(Ve,x)  
    
    # Thwaites method for laminar flow 
    theta_0    = 1E-12 # np.sqrt(0.075*l/dVe[0]/Re_L) # initial momentum thickness 
    H_lam,delta_star_lam,delta_lam,cf_lam,theta_lam,Re_x_lam,Re_theta_lam = solve_thwaites_BL(l,Re_L,x,Ve,dVe,theta_0) 
 
    # Find transition point 
    transition_index,transition_flag = michel_criterion(Re_theta_lam,Re_x_lam)
    
    if transition_flag:
        # Get initial conditons for turbulent boundary layer 
        delta_0_tur       = delta_lam[transition_index]
        delta_star_0_tur  = delta_star_lam[transition_index]
        theta_0_tur       = theta_lam[transition_index]
        Re_L_tur          = Re_L 
        x_tur             = x[transition_index:]
        Ve_tur            = Ve[transition_index:]
        dVe_tur           = dVe[transition_index:]
        
        # Turbulent Boundary Layer Calculation
        H_turb,delta_star_turb,delta_turb,cf_turb,theta_turb,Re_x_turb,Re_theta_turb = solve_heads_BL(l,delta_0_tur,theta_0_tur,delta_star_0_tur,Re_L_tur,x_tur,Ve_tur,dVe_tur)
        
        # Concatenate vectors 
        delta_star  = np.hstack((delta_star_lam[:transition_index],delta_star_turb))  
        delta       = np.hstack((delta_lam[:transition_index],delta_turb))
        cf          = np.hstack((cf_lam[:transition_index],cf_turb))
        theta       = np.hstack((theta_lam[:transition_index],theta_turb))
        Re_x        = np.hstack((Re_x_lam[:transition_index],Re_x_turb))
        Re_theta    = np.hstack((Re_theta_lam[:transition_index],Re_theta_turb) )  
        H           = np.hstack((H_lam[:transition_index],H_turb))
    else:
        delta_star = delta_star_lam
        delta      = delta_lam
        cf         = cf_lam
        theta      = theta_lam 
        Re_x       = Re_x_lam 
        Re_theta   = Re_theta_lam  
        H          = H_lam
        
    
    print('Transition Point: ', x[transition_index])
    Cf = np.sum(cf[1:]*np.diff(x))
    print(Cf)
    # Plots 
     
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
    axis1.plot(x,delta)  
    axis2.plot(x,delta_star)  
    axis3.plot(x,theta)  
    axis4.plot(x,cf)  
    axis5.plot(x,Re_x)  
    axis6.plot(x,H) 
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
 