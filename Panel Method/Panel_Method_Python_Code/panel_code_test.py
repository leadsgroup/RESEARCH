########################################################################
#
#  .py - Create surface panelization for a NACA4 Series Airfoil
#
#
#  Input list: 
#  N/A
#
#  Output list:
#  cl       -  Lift Coefficient  
#  cd       -  Drag Coefficient  
#  cm       -  Moment Coefficient  
#  x        -  Vector of x coordinates of the surface nodes
#  y        -  Vector of y coordinates of the surface nodes
#  cp       -  Pressure Coefficient Distribution
#
#  Written by: Matthew Clarke                                          #
#              Department of Aerospace Engineering                     #
#              University of Illinois, Urbana-Champaign                #
#              maclarke@illinois.edu                                   #
#                                                                      #
#  Last Modified: Wed Aug 30 2023
########################################################################

import matplotlib.pyplot as plt 
from os.path import dirname, realpath, sep, pardir
import sys
import numpy as np
sys.path.append(dirname(realpath(__file__)) + sep + pardir + sep + "lib")

from hess_smith             import hess_smith
from naca_4series_generator import naca_4series_generator  
from solve_thwaites_BL      import solve_thwaites_BL
from solve_heads_BL         import solve_heads_BL
from michel_criterion       import michel_criterion  

def main():
    problem_1()
    
    problem_2()
    
    return 
    
    
def problem_1(): 
    naca4  = '2410' 
    alpha  = 4
    npanel = 10      
    
    x, y = naca_4series_generator(naca4, npanel) 
    
    cl,cd,cm,cp, xbar,ybar,vt, ct, st = hess_smith(x,y,alpha)  
     
    # plot the output  
    fig = plt.figure('Problem 1 Part 1')
    fig.set_size_inches(10,7)
    axis_1 = fig.add_subplot(2, 1, 1)
    axis_1.plot(xbar, cp)
    axis_1.set_ylim([1,-2])
    axis_1.set_xlabel('x/c')
    axis_1.set_ylabel('Cp')
    axis_1.set_title('Coefficient of Pressure Distribution')
    axis_1.grid()
    
    axis_2 = fig.add_subplot(2, 1, 2)
    axis_2.plot(xbar, ybar)
    axis_2.plot(xbar, ybar, 'o')
    axis_2.set_xlabel('x/c')
    axis_2.set_ylabel('y/c')
    axis_2.set_title('Airfoil Geometry')
    axis_2.axis('equal')
    axis_2.grid()
     
 
    plt.tight_layout()      
    print('Airfoil           : NACA ' + naca4)
    print('Lift Coefficient  : ',cl)
    print('Drag Coefficient  : ',cd)
    print('Moment Coefficient: ',cm)
    
    alpha_vals = np.linspace(-5,10,16)
    CLs = np.zeros(len(alpha_vals))
    CDs = np.zeros(len(alpha_vals))
    CMs = np.zeros(len(alpha_vals)) 
    for i in range(len(alpha_vals)): 
        cl,cd,cm,cp, xbar,ybar,vt, ct, st = hess_smith(x,y,alpha_vals[i])  
        CLs[i] = cl
        CDs[i] = cd
        CMs[i] = cm
    
    
    
    fig_2 = plt.figure('Problem 1 Part 2')
    fig_2.set_size_inches(10,7)
    axis_1 = fig_2.add_subplot(2, 2, 1)
    axis_1.plot(alpha_vals,CLs)
    axis_1.set_title('Lift Coefficient')
    axis_1.set_xlabel(r'$\alpha$') 
    axis_2.set_ylabel(r'$C_l$')  
    axis_1.grid()
    
    axis_2 = fig_2.add_subplot(2, 2, 2)
    axis_2.plot(alpha_vals,CDs)
    axis_2.set_title('Drag Coefficient')
    axis_2.set_xlabel(r'$\alpha$')
    axis_2.set_ylabel(r'$C_d$')  
    axis_2.grid()

    axis_3 = fig_2.add_subplot(2, 2, 3)
    axis_3.plot(alpha_vals,CMs)
    axis_3.set_title('Moment Coefficient')
    axis_3.set_xlabel(r'$\alpha$')
    axis_3.set_ylabel(r'$C_m$') 
    axis_3.grid()       
    plt.tight_layout()
  
    return 

def problem_2():
    naca4           = '3310'
    alpha           = 3
    npanel          = 100 
    Re_L            = 5E6  

    # ---------------------------------------------------------------------------   
    # STEP 2: generate airfoil surface panelization 
    # ---------------------------------------------------------------------------  
    x, y = naca_4series_generator(naca4, npanel) 
    
    cl,cd,cm,cp, xbar,ybar,vt, ct, st = hess_smith(x,y,alpha)  
    
    # use tangential velocity to upper and lower surfaces 
    stagnation_point = np.where(vt>0)[0][0]  
    x_bot            = x[:stagnation_point][::-1]
    y_bot            = y[:stagnation_point][::-1] 
    x_top            = x[stagnation_point+1:]
    y_top            = y[stagnation_point+1:]  
     
    Ve_bot           = -vt[:stagnation_point][::-1]
    Ve_top           = vt[stagnation_point:]

    x_bl_bot         = np.zeros(len(Ve_bot))
    x_bl_top         = np.zeros(len(Ve_top))
    x_bl_bot[1:]     = np.cumsum(np.sqrt((x_bot[1:] - x_bot[:-1])**2 + (y_bot[1:] - y_bot[:-1])**2))
    x_bl_top[1:]     = np.cumsum(np.sqrt((x_top[1:] - x_top[:-1])**2 + (y_top[1:] - y_top[:-1])**2))
    
    dVe_bot          = np.zeros(len(Ve_bot))
    dVe_bot          = np.gradient(Ve_bot,x_bl_bot)
    dVe_top          = np.zeros(len(Ve_top))
    dVe_top          = np.gradient(Ve_top,x_bl_top)  
     
    # --------------------------------------------------------------------------------------------------------------------------------------------    
    # Top Surface 
    # --------------------------------------------------------------------------------------------------------------------------------------------
    theta_0_top = np.sqrt(0.075*x_bl_top[-1]/dVe_top[0]/Re_L)
    x_lam_top,H_lam_top,delta_star_lam_top,delta_lam_top,cf_lam_top,theta_lam_top,Re_x_lam_top,Re_theta_lam_top = \
        solve_thwaites_BL(x_bl_top[-1],Re_L,x_bl_top,Ve_top,dVe_top,theta_0_top) 
    

    # Find transition points for top and bottom surface 
    transition_index_top,transition_flag_top = michel_criterion(Re_theta_lam_top,Re_x_lam_top) 
 
    if transition_flag_top:
        # Get initial conditons for turbulent boundary layer 
        delta_0_tur_top       = delta_lam_top[transition_index_top]
        delta_star_0_tur_top  = delta_star_lam_top[transition_index_top]
        theta_0_tur_top       = theta_lam_top[transition_index_top]
        Re_L_tur_top          = Re_L  
        x_tur_top             = x_bl_top[transition_index_top:]
        Ve_tur_top            = Ve_top[transition_index_top:]
        dVe_tur_top           = dVe_top[transition_index_top:] 
        
        x_tur_top,H_tur_top,delta_star_tur_top,delta_tur_top,cf_tur_top,theta_tur_top,Re_x_tur_top,Re_theta_tur_top = \
            solve_heads_BL(x_bl_top[-1],delta_0_tur_top,theta_0_tur_top,delta_star_0_tur_top,Re_L_tur_top,x_tur_top,Ve_tur_top,dVe_tur_top)
    
        # Concatenate vectors 
        delta_star_top  = np.hstack((delta_star_lam_top[:transition_index_top],delta_star_tur_top))  
        delta_top       = np.hstack((delta_lam_top[:transition_index_top],delta_tur_top))
        cf_top          = np.hstack((cf_lam_top[:transition_index_top],cf_tur_top))
        theta_top       = np.hstack((theta_lam_top[:transition_index_top],theta_tur_top))
        Re_x_top        = np.hstack((Re_x_lam_top[:transition_index_top],Re_x_tur_top))
        Re_theta_top    = np.hstack((Re_theta_lam_top[:transition_index_top],Re_theta_tur_top) )  
        H_top           = np.hstack((H_lam_top[:transition_index_top],H_tur_top))
    else:
        delta_star_top = delta_star_lam_top
        delta_top      = delta_lam_top
        cf_top         = cf_lam_top
        theta_top      = theta_lam_top 
        Re_x_top       = Re_x_lam_top 
        Re_theta_top   = Re_theta_lam_top  
        H_top          = H_lam_top
       
    # --------------------------------------------------------------------------------------------------------------------------------------------    
    # Bottom Surface 
    # --------------------------------------------------------------------------------------------------------------------------------------------
    theta_0_bot = np.sqrt(0.075*x_bl_bot[-1]/dVe_bot[0]/Re_L)
    x_lam_bot,H_lam_bot,delta_star_lam_bot,delta_lam_bot,cf_lam_bot,theta_lam_bot,Re_x_lam_bot,Re_theta_lam_bot =\
        solve_thwaites_BL(x_bl_bot[-1],Re_L,x_bl_bot,Ve_bot,dVe_bot,theta_0_bot) 


    # Find transition points for bot and botton surface 
    transition_index_bot,transition_flag_bot = michel_criterion(Re_theta_lam_bot,Re_x_lam_bot) 

    if transition_flag_bot:
        # Get initial conditons for turbulent boundary layer 
        delta_0_tur_bot       = delta_lam_bot[transition_index_bot]
        delta_star_0_tur_bot  = delta_star_lam_bot[transition_index_bot]
        theta_0_tur_bot       = theta_lam_bot[transition_index_bot]
        Re_L_tur_bot          = Re_L  
        x_tur_bot             = x_bl_bot[transition_index_bot:]
        Ve_tur_bot            = Ve_bot[transition_index_bot:]
        dVe_tur_bot           = dVe_bot[transition_index_bot:] 

        x_tur_bot,H_tur_bot,delta_star_tur_bot,delta_tur_bot,cf_tur_bot,theta_tur_bot,Re_x_tur_bot,Re_theta_tur_bot = \
            solve_heads_BL(x_bl_bot[-1],delta_0_tur_bot,theta_0_tur_bot,delta_star_0_tur_bot,Re_L_tur_bot,x_tur_bot,Ve_tur_bot,dVe_tur_bot)

        # Concatenate vectors 
        delta_star_bot  = np.hstack((delta_star_lam_bot[:transition_index_bot],delta_star_tur_bot))  
        delta_bot       = np.hstack((delta_lam_bot[:transition_index_bot],delta_tur_bot))
        cf_bot          = np.hstack((cf_lam_bot[:transition_index_bot],cf_tur_bot))
        theta_bot       = np.hstack((theta_lam_bot[:transition_index_bot],theta_tur_bot))
        Re_x_bot        = np.hstack((Re_x_lam_bot[:transition_index_bot],Re_x_tur_bot))
        Re_theta_bot    = np.hstack((Re_theta_lam_bot[:transition_index_bot],Re_theta_tur_bot) )  
        H_bot           = np.hstack((H_lam_bot[:transition_index_bot],H_tur_bot))
    else:
        delta_star_bot = delta_star_lam_bot
        delta_bot      = delta_lam_bot
        cf_bot         = cf_lam_bot
        theta_bot      = theta_lam_bot 
        Re_x_bot       = Re_x_lam_bot 
        Re_theta_bot   = Re_theta_lam_bot  
        H_bot          = H_lam_bot
        
    # ------------------------------------------------------------------------------------------------------
    # Compute effective surface of airfoil with boundary layer  
    # ------------------------------------------------------------------------------------------------------   
    DELTA_STAR     = np.concatenate([delta_star_bot[::-1],delta_star_top])
    DELTA          = np.concatenate([delta_bot[::-1],delta_top])
    CF             = np.concatenate([cf_bot[::-1],cf_top])
    THETA          = np.concatenate([theta_bot[::-1],theta_top])
    RE_X           = np.concatenate([Re_x_bot[::-1],Re_x_top])
    RE_THETA       = np.concatenate([Re_theta_bot[::-1],Re_theta_top])
    H_BOT          = np.concatenate([H_bot[::-1],H_top])
    
    mask           = np.isnan(DELTA)
    DELTA[mask]    = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), DELTA[~mask])
    
    x_bl           = xbar + DELTA*-st
    y_bl           = ybar + DELTA*ct  
    
    cl_bl,cd_bl,cm_bl,cp_bl, xbar_bl,ybar_bl,vt_bl, ct_bl, st_bl = hess_smith(x_bl,y_bl,alpha) 
    
    print('inviscid cl: ', cl )
    print('viscous cl : ', cl_bl)
    print('inviscid cm: ', cm )
    print('viscous cm : ', cm_bl)

    #Plot velocity distributions
    fig  = plt.figure('Extra Credit Part 1')
    fig.set_size_inches(10,7)
    axis_1_1 = fig.add_subplot(3,3,1)     
    axis_1_1.plot(x_top,Ve_top,'r-',x_bot,Ve_bot,'b--') 
    axis_1_1.set_xlabel(r'$x$')
    axis_1_1.set_ylabel(r'$V_e$') 
    axis_1_1.set_xlim([0,1.1])
    axis_1_1.set_ylim([0,1.75])     
    axis_1_1.grid()        
 
    axis_1_2 = fig.add_subplot(3,3,2)     
    axis_1_2.plot(x_top,dVe_top,'r-',x_bot,dVe_bot,'b--') 
    axis_1_2.set_xlabel('$x$')
    axis_1_2.set_ylabel(r'$dV_e$/dx') 
    axis_1_2.set_xlim([0,1.1])
    axis_1_2.set_ylim([-5, 5])   
    axis_1_2.grid()    

    axis_1_3 = fig.add_subplot(3,3,3)     
    axis_1_3.plot(x_top,theta_top,'r-',x_bot,theta_bot,'b--') 
    axis_1_3.set_xlabel('$x$')
    axis_1_3.set_ylabel(r'$\theta$')   
    axis_1_3.set_xlim([0,1.1])
    axis_1_3.grid()   

    axis_1_4 = fig.add_subplot(3,3,4)     
    axis_1_4.plot(x_top,delta_top,'r-',x_bot,delta_bot,'b--') 
    axis_1_4.set_xlabel('$x$')
    axis_1_4.set_ylabel(r'$\delta$')   
    axis_1_4.set_ylim([0,0.05])   
    axis_1_4.set_xlim([0,1.1])
    axis_1_4.grid()   
    

    axis_1_5 = fig.add_subplot(3,3,5)     
    axis_1_5.plot(x_top,delta_star_top,'r-',x_bot,delta_star_bot,'b--') 
    axis_1_5.set_xlabel('$x$')
    axis_1_5.set_ylabel(r'$\delta$*')    
    axis_1_5.set_ylim([0,0.05])    
    axis_1_5.set_xlim([0,1.1])
    axis_1_5.grid()   

    axis_1_6 = fig.add_subplot(3,3,6)     
    axis_1_6.plot(x_top,cf_top,'r-',x_bot,cf_bot,'b--') 
    axis_1_6.set_xlabel('$x$')
    axis_1_6.set_ylabel(r'$c_f$')   
    axis_1_6.set_ylim([0, 0.01])  
    axis_1_6.set_xlim([0,1.1])  
    axis_1_6.grid()     
  
    axis_1_7 = fig.add_subplot(3,3,7)     
    axis_1_7.plot(x_top,H_top,'r-',x_bot,H_bot,'b--') 
    axis_1_7.set_xlabel('$x$')
    axis_1_7.set_ylabel(r'H')     
    axis_1_7.set_ylim([0, 5])  
    axis_1_7.set_xlim([0,1.1])  
    axis_1_7.grid()     

    axis_1_8 = fig.add_subplot(3,3,8)     
    axis_1_8.plot(x_top,Re_x_top,'r-',x_bot,Re_x_bot,'b--') 
    axis_1_8.set_xlabel('$x$')
    axis_1_8.set_ylabel(r'$Re_x$') 
    axis_1_8.grid()     
    axis_1_8.set_xlim([0,1.1]) 

    axis_1_9 = fig.add_subplot(3,3,9)     
    axis_1_9.plot(x_top,Re_theta_top,'r-',x_bot,Re_theta_bot,'b--') 
    axis_1_9.set_xlabel('$x$')
    axis_1_9.set_ylabel(r'$Re_{\theta}$')      
    axis_1_9.set_xlim([0,1.1])  
    axis_1_9.grid()    
    plt.tight_layout()  
    
    fig2  = plt.figure('Extra Credit Part 2')
    axis_2_1 = fig2.add_subplot(1,1,1)   
    axis_2_1.plot(x, y,'r-', label = 'airfoil')
    axis_2_1.plot(x_bl, y_bl,'b--', label = 'boundary layer')
    axis_2_1.set_xlabel(r'airfoil')
    axis_2_1.set_xlabel('x')  
    axis_2_1.set_ylabel('y')  
    axis_2_1.legend()
    axis_2_1.axis('equal')   
     

    # plot the output  
    fig3 = plt.figure('Extra Credit Part 3')
    axis_3_1 = fig3.add_subplot(1, 1, 1)
    axis_3_1.plot(xbar, cp,'r-', label = 'inviscid')
    axis_3_1.plot(xbar_bl, cp_bl,'b--', label = 'viscous') 
    axis_3_1.set_ylim([1,-2])    
    axis_3_1.set_xlabel('x/c')
    axis_3_1.set_ylabel('Cp')
    axis_3_1.set_title('Coefficient of Pressure Distribution')
    axis_3_1.legend()
    axis_3_1.grid()    
    

    plt.tight_layout()        
    return 

if __name__ == '__main__':
    main()
    plt.show()

