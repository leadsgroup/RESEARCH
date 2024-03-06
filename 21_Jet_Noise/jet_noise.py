
import numpy as np 
 

def main():
    LSF     = 0 # the linear scale factor
    X_C     = 0 # empirical length factor 
    theta_f = 0 #  empirical length factor 
    L       = 0 # NEED VALUE the displacement of the start of the mixing region from the fan nozzle exit plane 
    X_s     = (L + (X_C + theta/theta_f)*D)/LDF # eq. 4 source location for a mixing region noise,
    
    R_cor_dif_R_sqrd = 1 + (X_s/R)**2 + 2*(X_s/R)*np.cos(theta)  # 1 
    R_cor_dif_R      = np.sqrt(R_cor_dif_R_sqrd)
    R_cor            = R_cor_dif_R*R
    theta_cor    = np.arccos( np.sqrt(R_cor_dif_R_sqrd)*np.cos(theta)  +  (X_s/R_cor) ) # 1.a 
    delta_theta = theta_cor - theta            # 1.b 
    delta_SPL   = -20*np.log10(R_cor_dif_R)      # 1.c 
    
    

    # -----------------------------------------------------    
    # Mixing noise compoments  
    # -----------------------------------------------------
    C = 0 # NEED VALUE coefficient obtained experimentally 
    N = 0 # NEED VALU Eslope obtained experimentally 
    A  = 0 # NEED VALUE approprite nozzle exia area
    rho = 0 # NEED VALUE  fully exppanded jet density
    k  = 4
    n_c  = 0.6 
    alpha = 0.2 

    V     = 0 #  characteristic velocity 
    V_e   = V*(1 - M_f*(c_amb/V))**(1/2) # eqn 2.a  
    
    M_c   = n_c*((V_e/c_amb) - M_f) # 2.c
    omega  = 3*((V_e/c_amb)**(3.5))/(0.6 + (V_e/c_amb)**(3/5)) - 1 # 2.b
    
    UOL  = C + 10*np.log10( (rho_amb/rho_ISA)**2 *  (c_amb/c_ISA)**4) + 10*np.log10(A/(R_cor**2)) + 10*omega*np.log10(rho/rho_amb) + \
                10*np.log10( ((V_e/c_amb)**N) /(1 + b*(V_e/c_amb))**(N-3) ) - 5*k*np.log10(  (1 + M_c*np.cos(theta_cor))**2 + (alpha**2)*(M_c**2))
    
    theta_prime = theta_cor*((V/c_amb)**0.1)# NEED VALUE effective directivity angle 
    
    # strouhal number 
    D = (4*A/np.pi)**(0.5) # characteristic diameter 
    T = 0  # total temperature 
    S = (f*D/V_e)*((T/_Tamb)**(0.4*(1-np.cos(theta_prime))))
    

    # -----------------------------------------------------    
    # Shock Noise Compoments 
    # ----------------------------------------------------- 
    # not included as of now 
    

    # -----------------------------------------------------    
    # Large Scale Mixing Noise Compoments 
    # ----------------------------------------------------- 

    V_O      = 0 # NEED VALUE 
    V_I      = 0 # NEED VALUE
    m_dot_I  = 0 # NEED VALUE
    m_dot_O  = 0 # NEED VALUE
    

    t_mix    = T_mix - (gamma-1)*(V_mix**2)/(2*gamma*R)
    rho_L    = P_amb/(R*t_mix)
    if dual_stream: 
        n_cL    = 0.4/(1 + 1.5(*(V_O/V_I)**4)) + 0.3
        A_L     = A_I + A_O
        V_mix   = (m_dot_I*V_I + m_dot_O*V_O)/(m_dot_I + m_dot_O )
        T_mix   = (m_dot_I*T_I + m_dot_O*T_O)/(m_dot_I + m_dot_O )
    else: 
        n_cL    = 0.46 
        A_L     = A_I
        V_mix   = V_I  
        T_mix   = T_I   
        
    M_cL    = n_cL*((V_mix/c_amv) - M_f) # 8.c 
    V_eL    = V_mix*(1- M_f*(c_amb/V_mix))**(1/2)  # 8.b 
    omega_L = 3*((V_eL/c_amb)**(3.5))/(0.6 + (V_eL/c_amb)**(3/5)) - 1 # 8.a 
    
    # the normalized baseline (unsuppressed) large scale mixing noise overall level
    UOL_L_norm  = UOL_L - 10*np.log10((rho_amb/rho_ISA)**2 * (c_amb/c_ISA)**4) - 10*np.log10(A_L/(R_corL**2)) -\
        10*omega_L**np.log10(rho_L/rho_amb) + 20*np.log10((1 + M_cL*np.cos(theta_cor_L))**2 + 0.09*(M_cL**2))
    
    
    UOL_L_norm = 140 + 10*log10( ((V_eL/c_amb)**8)/(1 + 0.02*((V_eL/c_amb**5) )      )
    return 