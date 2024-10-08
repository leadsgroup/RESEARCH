import numpy as np
import matplotlib.pyplot as plt

def main():
    
    turboprop_model()
    
    return

def turboprop_model():
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Complete set of equations for the Ideal Case 
    #  Reference: Elements or gas turbine propulsion/Jack D. Mallingly, pag 322
    # ----------------------------------------------------------------------------------------------------------------------    
    
    #eta_P                = 2*V0/(V_jet + V0)
    #F_prop               = m_dot_0*(V_jet - V0)/g_c            
    #eta_prop             = F_prop*V0/W_prop_dot
    #eta_L                = (0.25*m_dot_prop*(V_jet**2 - V0**2))/W_dot_prop
    #eta_prop             = eta_P*eta_L
    #Cc                   = F_c*V_0/(m_dot_0*cp*T0)
    #C_tot                = Cc + C_prop
    #F                    = F_prop + F_c 
    #F                    = (C_tot*m_dot_0*cp*T0)/V0
    #pi_tH                = P_t45/P_t4
    #tau_tH               = T_t45/T_t4
    #pi_tL                = P_t5/P_t45
    #tau_tL               = T_t5/T_t45
    #F_c_m_dot_0          = (a0/g_c)*(V9_a9 - M0)
    #Cc                   = (V0/(cp*T0))*(a0/g_c)*(V9_a0 - M0)
    #V9_a0                = np.sqrt(T9/T0)*M9
    #M9                   = np.sqrt((2/(gamma - 1))*((P_t9/P_9)**((gamma - 1)/gamma) - 1))
    #Pt9/P9               = (P0/P9)*pi_r*pi_d*pi_c*pi_b*pi_tH*pi_tL*pi_n
    #Pt9/P9               = pi_r*pi_c*pi_tH*pi_tL
    #M9                   = np.sqrt((2/(gamma - 1))*(tau_r*tau_c_r*tau_tH*tau_tL - 1))
    #T9_T0                = tau_lambda/tau_r*tau_c_r
    #f                    = ((cp*T0)/h_PR)*(tau_lambda - tau_r*tau_c_r)
    #tau_tH               = 1 - (tau_r/tau_lambda)*(tau_c_r - 1)
    #W_dot_prop           = m_dot_45*cp*(Tt45 - Tt5)
    #C_prop               = (eta_prop*W_dot_prop)/(m_dot_0*cp*T0) 
    #C_prop               = eta_prop*tau_lambda*tau_tH*(1 - tau_tL)
    #W_dot_m_dot_0        = C_tot*cp*T0
    #F_m_dot_0            = (C_tot*cp*T0)/V0
    #S                    = f/(F/m_dot_0)
    #S_P                  = m_dot_f/W_dot
    #S_P                  = f/(W_dot/m_dot_0)
    #S_P                  = f/(C_tot*cp*T0)
    #eta_T                = 1 - 1/(tau_r*tau_c_r)
    #eta_O                = W_dot/(m_dot_f*h_PR)
    #eta_O                = C_tot/(tau_lambda - tau_r*tau_c_r)    
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Input values for ideal turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    M0_i                  = 0.8                       # [-]
    T0_i                  = 240                       # [K]
    gamma                 = 1.4                       # [-]
    cp                    = 1004                      # [J/(kg*K)]
    h_PR                  = 42800000                  # [J/kg]
    Tt4                   = 1370                      # [K]
    pi_c                  = np.linspace(2.5,40,100)   # [-]
    tau_t                 = [0.4, 0.5, 0.6, 0.7, 0.8] # [-]
    eta_prop              = 0.83                      # [-]
     
    # ----------------------------------------------------------------------------------------------------------------------
    # Compiutation for ideal turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    F_mdot0_i             = np.zeros((len(pi_c), len(tau_t)))                    
    tau_c_i               = np.zeros_like(F_mdot0_i)
    f_i                   = np.zeros_like(F_mdot0_i)
    tau_tH_i              = np.zeros_like(F_mdot0_i)
    tau_tL_i              = np.zeros_like(F_mdot0_i)
    V9_a0_i               = np.zeros_like(F_mdot0_i)
    Cc_i                  = np.zeros_like(F_mdot0_i)
    C_prop_i              = np.zeros_like(F_mdot0_i)
    C_tot_i               = np.zeros_like(F_mdot0_i)
    S_i                   = np.zeros_like(F_mdot0_i)
    eta_T_i               = np.zeros_like(F_mdot0_i)
    eta_O_i               = np.zeros_like(F_mdot0_i)
    eta_P_i               = np.zeros_like(F_mdot0_i)
    
    for i in range(len(pi_c)):
        
        for j in range(len(tau_t)):
            
            R                    = ((gamma - 1)/gamma)*cp
            a0                   = np.sqrt(gamma*R*T0_i)
            tau_r                = 1 + ((gamma - 1)/2)*M0_i**2
            tau_lambda           = Tt4/T0_i
            tau_c_i[i,j]         = pi_c[i]**((gamma - 1)/gamma)
            f_i[i,j]             = ((cp*T0_i)/h_PR)*(tau_lambda - tau_r*tau_c_i[i,j])
            tau_tH_i[i,j]        = 1 - (tau_r/tau_lambda)*(tau_c_i[i,j] - 1)
            tau_tL_i[i,j]        = tau_t[j]/tau_tH_i[i,j]
            V9_a0_i[i,j]         = np.sqrt((2/(gamma - 1))*(tau_lambda*tau_t[j] - tau_lambda/(tau_r*tau_c_i[i,j])))
            Cc_i[i,j]            = (gamma - 1)*M0_i*(V9_a0_i[i,j] - M0_i)
            C_prop_i[i,j]        = eta_prop*tau_lambda*tau_tH_i[i,j]*(1 - tau_tL_i[i,j])
            C_tot_i[i,j]         = Cc_i[i,j] + C_prop_i[i,j]
            F_mdot0_i[i,j]       = (C_tot_i[i,j]*cp*T0_i)/(M0_i*a0)
            S_i[i,j]             = (f_i[i,j]/(F_mdot0_i[i,j]))*1000**2
            eta_T_i[i,j]         = 1 - 1/(tau_r*tau_c_i[i,j])
            eta_O_i[i,j]         = C_tot_i[i,j]/(tau_lambda - tau_r*tau_c_i[i,j])
            eta_P_i[i,j]         = eta_O_i[i,j]/eta_T_i[i,j]  
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Complete set of equations for Parametric Cycle Analysis of Real Case
    #  Reference: Elements or gas turbine propulsion/Jack D. Mallingly, pag 433
    # ----------------------------------------------------------------------------------------------------------------------    
    
    #F_c_m_dot_0          = (a0/g_c)*((1 + f)*V9_a0 - M0 + (1 + f)*(Rt/Rc)*((T9/T0)/(V9/a0))*((1 - (P0/P9))/gamma_c))
    #Cc                   = (V0/(cp*T0))*((1 + f)*V9_a0 - M0 + (1 + f)*(Rt/Rc)*((T9/T0)/(V9/a0))*((1 - (P0/P9))/gamma_c))
    #Cc                   = (gamma_c - 1)*M0*((1 + f)*V9_a0 - M0 + (1 + f)*(Rt/Rc)*((T9/T0)/(V9/a0))*((1 - (P0/P9))/gamma_c))
    #(V9_a0)**2           = ((gamma_t*R_t*T_9)/(gamma_c*R_c*T0))*M9**2
    #M9**2                = (2/(gamma_t - 1))*((P_t9/P_9)**((gamma_t - 1)/gamma_t) - 1)
    #Pt9/P9               = (P0/P9)*pi_r*pi_d*pi_c*pi_b*pi_tH*pi_tL*pi_n
    #For unchoked flow:
    # Pt9/P0               <= ((gamma_t +1)/2)**(gamma_t/(gamma_t - 1))
    #For choked flow:
    # Pt9/P0               > ((gamma_t +1)/2)**(gamma_t/(gamma_t - 1))
    # M9                   = 1
    # Pt9/P9               = ((gamma_t +1)/2)**(gamma_t/(gamma_t - 1))
    # P0/P9                = ((Pt9/P9)/(Pt9/P0))
    #T9/T0                = (Tt9/T0)/((Pt9/P9)**((gamma_t - 1)/gamma_t))
    #Tt9/T0               = (Tt4/T0)*tau_tH*tau_tL
    #V9/a0                = np.sqrt(((2*tau_lambda*tau_tH*tau_tL)/(gamma_c - 1))*(1 - (Pt9/P9)**((gamma_t - 1)/gamma_t))))
    #f                    = (tau_lambda - tau_r*tau_c_r)/(eta_b*h_PR/(cp*T0) - tau_lambda)
    #tau_tH               = 1 - (1/(eta_mH*(1 + f)))*(tau_r/tau_lambda)*(tau_c_r - 1)
    #W_dot_prop           = eta_g*eta_mL*m_dot_45*cpt*(Tt45 - Tt5)
    #C_prop               = (eta_prop*W_dot_prop)/(m_dot_0*cp*T0) 
    #C_prop               = eta_prop*eta_g*eta_mL*(1 + f)*tau_lambda*tau_tH*(1 - tau_tL)
    #W_dot_m_dot_0        = C_tot*cpc*T0
    #F_m_dot_0            = (C_tot*cpc*T0)/V0
    #S                    = f/(F/m_dot_0)    
    #S_P                  = f/(C_tot*cpc*T0)
    #eta_T                = (m_dot_0*cpc*T0*C_tot)/(m_dot_f*h_PR)
    #eta_T                = C_tot/(f*h_PR/(cpc*T0))
    #eta_P                = (m_dot_0*cpc*T0*C_tot)/(W_dot_prop*(m_dot_9*V_9**2 - m_dot_0*V_0**2)/(2*g_c)) 
    #eta_P                = C_tot/(C_tot/eta_prop + ((gamma - 1)/2)*((1 + f)*(V9/a0)**2 - M0**2))    
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Input values for real turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    # 1) Real engine
    # 2) Eng library James Aircraft for Validation
    # 3) Integration into RCAIDE
    # 4) Get ATR 72 / Otter Engine data
    
    M0_r                  = 0.8                       # [-]
    T0_r                  = 240                       # [K]
    gamma_c               = 1.4                       # [-]
    cp_c                  = 1004                      # [J/(kg*K)]
    gamma_t               = 1.35                      # [-]
    cp_t                  = 1108                      # [J/(kg*K)]    
    h_PR                  = 42800000                  # [J/kg]
    e_c                   = 0.9                       # [-]
    e_tH                  = 0.89                      # [-]
    e_tL                  = 0.91                      # [-]
    pi_d                  = 0.98                      # [-]
    pi_b                  = 0.96                      # [-] 
    pi_n                  = 0.99                      # [-]
    eta_b                 = 0.99                      # [-]
    eta_mH                = 0.99                      # [-] 
    eta_mL                = 0.99                      # [-] 
    eta_prop              = 0.83                      # [-]
    eta_g                 = 0.99                      # [-]    
    Tt4                   = 1370                      # [K]
    pi_c                  = np.linspace(2.5,40,100)   # [-]
    tau_t                 = [0.4, 0.5, 0.6, 0.7, 0.8] # [-]
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Compiutation for real turboprop
    # ----------------------------------------------------------------------------------------------------------------------
                      
    tau_c_r               = np.zeros((len(pi_c), len(tau_t)))    
    eta_c                 = np.zeros_like(tau_c_r)
    f_r                   = np.zeros_like(tau_c_r)
    tau_tH_r              = np.zeros_like(tau_c_r)
    pi_tH                 = np.zeros_like(tau_c_r)
    eta_tH                = np.zeros_like(tau_c_r)
    tau_tL_r              = np.zeros_like(tau_c_r)
    pi_tL                 = np.zeros_like(tau_c_r)
    eta_tL                = np.zeros_like(tau_c_r) 
    Pt9_P0                = np.zeros_like(tau_c_r)
    M9                    = np.zeros_like(tau_c_r)
    Pt9_P9                = np.zeros_like(tau_c_r)
    P0_P9                 = np.zeros_like(tau_c_r)
    Tt9_T0                = np.zeros_like(tau_c_r)
    T9_T0                 = np.zeros_like(tau_c_r)
    V9_a0_r               = np.zeros_like(tau_c_r)
    C_prop_r              = np.zeros_like(tau_c_r)
    Cc_r                  = np.zeros_like(tau_c_r)  
    C_tot_r               = np.zeros_like(tau_c_r)
    F_mdot0_r             = np.zeros_like(tau_c_r)
    S_r                   = np.zeros_like(tau_c_r)
    W_dot_mdot0           = np.zeros_like(tau_c_r)
    S_P                   = np.zeros_like(tau_c_r)
    eta_T_r               = np.zeros_like(tau_c_r)
    eta_P_r               = np.zeros_like(tau_c_r)
    
    
    for i in range(len(pi_c)):
        for j in range(len(tau_t)):
            
            R_c                  = ((gamma_c - 1)/gamma_c)*cp_c
            R_t                  = ((gamma_t - 1)/gamma_t)*cp_t
            a0                   = np.sqrt(gamma_c*R_c*T0_r)
            V0                   = a0*M0_r
            tau_r                = 1 + ((gamma_c - 1)/2)*M0_r**2
            pi_r                 = tau_r**(gamma_c/(gamma_c - 1))
            if M0_r <= 1:
                eta_r            = 1
            else:
                eta_r            = 1 - 0.0075(M0_r - 1)**1.35
            tau_lambda           = (cp_t*Tt4)/(cp_c*T0_r)
            tau_c_r[i,j]         = pi_c[i]**((gamma_c - 1)/(gamma_c*e_c))
            eta_c[i,j]           = (pi_c[i]**((gamma_c - 1)/(gamma_c)) - 1)/(tau_c_r[i,j] - 1)
            f_r[i,j]             = (tau_lambda - tau_r*tau_c_r[i,j])/((eta_b*h_PR)/(cp_c*T0_r) - tau_lambda)
            tau_tH_r[i,j]        = 1 - ((tau_r*(tau_c_r[i,j] - 1))/(eta_mH*(1 + f_r[i,j])*tau_lambda))
            pi_tH[i,j]           = tau_tH_r[i,j]**(gamma_t/((gamma_t - 1)*e_tH))
            eta_tH[i,j]          = (1 - tau_tH_r[i,j])/(1 - tau_tH_r[i,j]**(1/e_tH))
            tau_tL_r[i,j]        = tau_t[j]/tau_tH_r[i,j]
            pi_tL[i,j]           = tau_tL_r[i,j]**(gamma_t/((gamma_t - 1)*e_tL))
            eta_tL[i,j]          = (1 - tau_tL_r[i,j])/(1 - tau_tL_r[i,j]**(1/e_tL))
            Pt9_P0[i,j]          = pi_r*pi_d*pi_c[i]*pi_b*pi_tH[i,j]*pi_tL[i,j]*pi_n
            
            if Pt9_P0[i,j] > ((gamma_t + 1)/2)**(gamma_t/(gamma_t - 1)):
                M9[i,j]          = 1
                Pt9_P9[i,j]      = ((gamma_t + 1)/2)**(gamma_t/(gamma_t - 1))
                P0_P9[i,j]       = ((Pt9_P9[i,j])/(Pt9_P0[i,j]))
            else:
                Pt9_P9[i,j]      = Pt9_P0[i,j]
                M9[i,j]          = np.sqrt((2/(gamma_t - 1))*((Pt9_P0[i,j])**((gamma_t - 1)/gamma_t)) - 1)
                P0_P9[i,j]       = 1
            Tt9_T0[i,j]          = (Tt4/T0_r)*tau_tH_r[i,j]*tau_tL_r[i,j]
            T9_T0[i,j]           = (Tt9_T0[i,j])/((Pt9_P9[i,j])**((gamma_t - 1)/gamma_t))
            V9_a0_r[i,j]         = np.sqrt(((2*tau_lambda*tau_tH_r[i,j]*tau_tL_r[i,j])/(gamma_c - 1))*(1 - Pt9_P9[i,j]**(-(gamma_t - 1)/(gamma_t))))
            C_prop_r[i,j]        = eta_prop*eta_g*eta_mL*(1 + f_r[i,j])*tau_lambda*tau_tH_r[i,j]*(1 - tau_tL_r[i,j])
            Cc_r[i,j]            = (gamma_c - 1)*M0_r*((1 + f_r[i,j])*V9_a0_r[i,j] - M0_r + (1 + f_r[i,j])*(R_t/R_c)*((T9_T0[i,j])/(V9_a0_r[i,j]))*((1 - (P0_P9[i,j]))/gamma_c))
            C_tot_r[i,j]         = Cc_r[i,j] + C_prop_r[i,j]
            F_mdot0_r[i,j]       = (C_tot_r[i,j]*cp_c*T0_r)/(V0)
            S_r[i,j]             = (f_r[i,j]/(F_mdot0_r[i,j]))*1000**2
            W_dot_mdot0[i,j]     = C_tot_r[i,j]*cp_c*T0_r
            S_P[i,j]             = (f_r[i,j]/(C_tot_r[i,j]*cp_c*T0_r))
            eta_T_r[i,j]         = C_tot_r[i,j]/((f_r[i,j]*h_PR)/(cp_c*T0_r))
            eta_P_r[i,j]         = C_tot_r[i,j]/((C_prop_r[i,j]/eta_prop) + ((gamma_c - 1)/2)*((1 + f_r[i,j])*(V9_a0_r[i,j])**2 - M0_r**2))  
              
    # ----------------------------------------------------------------------------------------------------------------------
    #  Complete set of equations for Parametric Cycle Analysis of Real Case
    #  Reference: Elements or gas turbine propulsion/Jack D. Mallingly, pag 560
    # ----------------------------------------------------------------------------------------------------------------------            
    
    #mdot_0                       = (P_0*pi_r*pi_d*pi_c)/()                         
        
    # ----------------------------------------------------------------------------------------------------------------------
    #  Results for ideal turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    fig1    =  plt.figure('F_mdot0_i')
    fig1.set_size_inches(7, 6)
    axis_1 = fig1.add_subplot(1,1,1)
    axis_1.plot(pi_c,F_mdot0_i[:,0], label='tau_t: 0.4')
    axis_1.plot(pi_c,F_mdot0_i[:,1], label='tau_t: 0.5') 
    axis_1.plot(pi_c,F_mdot0_i[:,2], label='tau_t: 0.6') 
    axis_1.plot(pi_c,F_mdot0_i[:,3], label='tau_t: 0.7') 
    axis_1.plot(pi_c,F_mdot0_i[:,4], label='tau_t: 0.8') 
    axis_1.set_xlabel('pi_c [-]')
    axis_1.set_ylabel('F_m_dot_0 [N/(kg/sec)]')
    axis_1.legend()
    axis_1.set_ylim(800, 1600)
    fig1.tight_layout()
    
    fig2    =  plt.figure('S_i')
    fig2.set_size_inches(7, 6)
    axis_2 = fig2.add_subplot(1,1,1)
    axis_2.plot(pi_c,S_i[:,0], label='tau_t: 0.4')
    axis_2.plot(pi_c,S_i[:,1], label='tau_t: 0.5') 
    axis_2.plot(pi_c,S_i[:,2], label='tau_t: 0.6') 
    axis_2.plot(pi_c,S_i[:,3], label='tau_t: 0.7') 
    axis_2.plot(pi_c,S_i[:,4], label='tau_t: 0.8')     
    axis_2.set_xlabel('pi_c [-]')
    axis_2.set_ylabel('S [mg/(N*sec)]')
    axis_2.legend()
    axis_2.set_ylim(10, 26)
    fig2.tight_layout()   
    
    fig3    =  plt.figure('eta_T')
    fig3.set_size_inches(7, 6)
    axis_3 = fig3.add_subplot(1,1,1)
    axis_3.plot(pi_c,eta_T_i[:,0], label='tau_t: 0.4')
    axis_3.plot(pi_c,eta_T_i[:,1], label='tau_t: 0.5') 
    axis_3.plot(pi_c,eta_T_i[:,2], label='tau_t: 0.6') 
    axis_3.plot(pi_c,eta_T_i[:,3], label='tau_t: 0.7') 
    axis_3.plot(pi_c,eta_T_i[:,4], label='tau_t: 0.8')   
    axis_3.set_xlabel('pi_c [-]')
    axis_3.set_ylabel('eta_T [-]')
    axis_3.legend()
    axis_3.set_ylim(0.20, 0.90)
    fig3.tight_layout()    
    
    fig4    =  plt.figure('eta_O')
    fig4.set_size_inches(7, 6)
    axis_4 = fig4.add_subplot(1,1,1)
    axis_4.plot(pi_c,eta_O_i[:,0], label='tau_t: 0.4')
    axis_4.plot(pi_c,eta_O_i[:,1], label='tau_t: 0.5') 
    axis_4.plot(pi_c,eta_O_i[:,2], label='tau_t: 0.6') 
    axis_4.plot(pi_c,eta_O_i[:,3], label='tau_t: 0.7') 
    axis_4.plot(pi_c,eta_O_i[:,4], label='tau_t: 0.8')     
    axis_4.set_xlabel('pi_c [-]')
    axis_4.set_ylabel('eta_O [-]')
    axis_4.legend()
    axis_4.set_ylim(0.20, 0.90)
    fig4.tight_layout()     
    
    fig5    =  plt.figure('eta_P')
    fig5.set_size_inches(7, 6)
    axis_5 = fig5.add_subplot(1,1,1)
    axis_5.plot(pi_c,eta_P_i[:,0], label='tau_t: 0.4')
    axis_5.plot(pi_c,eta_P_i[:,1], label='tau_t: 0.5') 
    axis_5.plot(pi_c,eta_P_i[:,2], label='tau_t: 0.6') 
    axis_5.plot(pi_c,eta_P_i[:,3], label='tau_t: 0.7') 
    axis_5.plot(pi_c,eta_P_i[:,4], label='tau_t: 0.8')     
    axis_5.set_xlabel('pi_c [-]')
    axis_5.set_ylabel('eta_P [-]')
    axis_5.legend()
    axis_5.set_ylim(0.20, 0.90)
    fig5.tight_layout()          
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Results for real turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    fig6    =  plt.figure('F_mdot0_r')
    fig6.set_size_inches(7, 6)
    axis_6 = fig6.add_subplot(1,1,1)
    axis_6.plot(pi_c,F_mdot0_r[:,0], label='tau_t: 0.4')
    axis_6.plot(pi_c,F_mdot0_r[:,1], label='tau_t: 0.5') 
    axis_6.plot(pi_c,F_mdot0_r[:,2], label='tau_t: 0.6') 
    axis_6.plot(pi_c,F_mdot0_r[:,3], label='tau_t: 0.7') 
    axis_6.plot(pi_c,F_mdot0_r[:,4], label='tau_t: 0.8')     
    axis_6.set_xlabel('pi_c [-]')
    axis_6.set_ylabel('F_m_dot_0 [N/(kg/sec)]')
    axis_6.legend()
    axis_6.set_xlim(0, 40)
    axis_6.set_ylim(500, 1300)
    fig6.tight_layout()
    
    fig7    =  plt.figure('S_r')
    fig7.set_size_inches(7, 6)
    axis_7 = fig7.add_subplot(1,1,1)
    axis_7.plot(pi_c,S_r[:,0], label='tau_t: 0.4')
    axis_7.plot(pi_c,S_r[:,1], label='tau_t: 0.5') 
    axis_7.plot(pi_c,S_r[:,2], label='tau_t: 0.6') 
    axis_7.plot(pi_c,S_r[:,3], label='tau_t: 0.7') 
    axis_7.plot(pi_c,S_r[:,4], label='tau_t: 0.8')     
    axis_7.set_xlabel('pi_c [-]')
    axis_7.set_ylabel('S [mg/(N*sec)]')
    axis_7.legend()
    axis_7.set_xlim(0, 40)
    axis_7.set_ylim(15, 40)
    fig7.tight_layout()    
    
    fig8    =  plt.figure('C_tot')
    fig8.set_size_inches(7, 6)
    axis_8 = fig8.add_subplot(1,1,1)
    axis_8.plot(pi_c,C_tot_r[:,0], label='tau_t: 0.4')
    axis_8.plot(pi_c,C_tot_r[:,1], label='tau_t: 0.5') 
    axis_8.plot(pi_c,C_tot_r[:,2], label='tau_t: 0.6') 
    axis_8.plot(pi_c,C_tot_r[:,3], label='tau_t: 0.7') 
    axis_8.plot(pi_c,C_tot_r[:,4], label='tau_t: 0.8')    
    axis_8.set_xlabel('pi_c [-]')
    axis_8.set_ylabel('C_tot [-]')
    axis_8.legend()
    axis_8.set_ylim(0, 1.4)
    fig8.tight_layout()    
    
    fig9    =  plt.figure('Cc')
    fig9.set_size_inches(7, 6)
    axis_9 = fig9.add_subplot(1,1,1)
    axis_9.plot(pi_c,Cc_r[:,0], label='tau_t: 0.4')
    axis_9.plot(pi_c,Cc_r[:,1], label='tau_t: 0.5') 
    axis_9.plot(pi_c,Cc_r[:,2], label='tau_t: 0.6') 
    axis_9.plot(pi_c,Cc_r[:,3], label='tau_t: 0.7') 
    axis_9.plot(pi_c,Cc_r[:,4], label='tau_t: 0.8')    
    axis_9.set_xlabel('pi_c [-]')
    axis_9.set_ylabel('Cc [-]')
    axis_9.legend()
    axis_9.set_ylim(0, 1.4)
    fig9.tight_layout()     
    
    fig10    =  plt.figure('C_prop')
    fig10.set_size_inches(7, 6)
    axis_10 = fig10.add_subplot(1,1,1)
    axis_10.plot(pi_c,C_prop_r[:,0], label='tau_t: 0.4')
    axis_10.plot(pi_c,C_prop_r[:,1], label='tau_t: 0.5') 
    axis_10.plot(pi_c,C_prop_r[:,2], label='tau_t: 0.6') 
    axis_10.plot(pi_c,C_prop_r[:,3], label='tau_t: 0.7') 
    axis_10.plot(pi_c,C_prop_r[:,4], label='tau_t: 0.8')     
    axis_10.set_xlabel('pi_c [-]')
    axis_10.set_ylabel('C_prop [-]')
    axis_10.legend()
    axis_10.set_ylim(0, 1.4)
    fig10.tight_layout()  
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Results for ideal vs real turboprop
    # ----------------------------------------------------------------------------------------------------------------------    
    
    fig11    =  plt.figure('F_mdot0 ideal vs real')
    fig11.set_size_inches(7, 6)
    axis_11 = fig11.add_subplot(1,1,1)
    axis_11.plot(pi_c,F_mdot0_i[:,0], label='tau_t ideal: 0.4', color = 'b', linestyle = '--')
    axis_11.plot(pi_c,F_mdot0_r[:,0], label='tau_t real: 0.4',  color = 'b', linestyle = '-')
    axis_11.plot(pi_c,F_mdot0_i[:,1], label='tau_t ideal: 0.5', color = 'm', linestyle = '--') 
    axis_11.plot(pi_c,F_mdot0_r[:,1], label='tau_t real: 0.5',  color = 'm', linestyle = '-') 
    axis_11.plot(pi_c,F_mdot0_i[:,2], label='tau_t ideal: 0.6', color = 'c', linestyle = '--') 
    axis_11.plot(pi_c,F_mdot0_r[:,2], label='tau_t real: 0.6',  color = 'c', linestyle = '-') 
    axis_11.plot(pi_c,F_mdot0_i[:,3], label='tau_t ideal: 0.7', color = 'k', linestyle = '--') 
    axis_11.plot(pi_c,F_mdot0_r[:,3], label='tau_t real: 0.7',  color = 'k', linestyle = '-') 
    axis_11.plot(pi_c,F_mdot0_i[:,4], label='tau_t ideal: 0.8', color = 'g', linestyle = '--') 
    axis_11.plot(pi_c,F_mdot0_r[:,4], label='tau_t real: 0.8',  color = 'g', linestyle = '-')     
    axis_11.set_xlabel('pi_c [-]')
    axis_11.set_ylabel('F_m_dot_0 [N/(kg/sec)]')
    axis_11.legend()
    axis_11.set_xlim(0, 40)
    axis_11.set_ylim(500, 1600)
    fig11.tight_layout()   
    
    fig12    =  plt.figure('S ideal vs real')
    fig12.set_size_inches(7, 6)
    axis_12 = fig12.add_subplot(1,1,1)
    axis_12.plot(pi_c,S_i[:,0], label='tau_t ideal: 0.4', color = 'b', linestyle = '--')
    axis_12.plot(pi_c,S_r[:,0], label='tau_t real: 0.4',  color = 'b', linestyle = '-')
    axis_12.plot(pi_c,S_i[:,1], label='tau_t ideal: 0.5', color = 'm', linestyle = '--') 
    axis_12.plot(pi_c,S_r[:,1], label='tau_t real: 0.5',  color = 'm', linestyle = '-') 
    axis_12.plot(pi_c,S_i[:,2], label='tau_t ideal: 0.6', color = 'c', linestyle = '--') 
    axis_12.plot(pi_c,S_r[:,2], label='tau_t real: 0.6',  color = 'c', linestyle = '-') 
    axis_12.plot(pi_c,S_i[:,3], label='tau_t ideal: 0.7', color = 'k', linestyle = '--') 
    axis_12.plot(pi_c,S_r[:,3], label='tau_t real: 0.7',  color = 'k', linestyle = '-') 
    axis_12.plot(pi_c,S_i[:,4], label='tau_t ideal: 0.8', color = 'g', linestyle = '--') 
    axis_12.plot(pi_c,S_r[:,4], label='tau_t real: 0.8',  color = 'g', linestyle = '-')     
    axis_12.set_xlabel('pi_c [-]')
    axis_12.set_ylabel('S [mg/(N*sec)]')
    axis_12.legend()
    axis_12.set_xlim(0, 40)
    axis_12.set_ylim(10, 40)
    fig12.tight_layout() 
    
    fig13    =  plt.figure('eta_T ideal vs real')
    fig13.set_size_inches(7, 6)
    axis_13 = fig13.add_subplot(1,1,1)
    axis_13.plot(pi_c,eta_T_i[:,0], label='tau_t ideal: 0.4', color = 'b', linestyle = '--')
    axis_13.plot(pi_c,eta_T_r[:,0], label='tau_t real: 0.4',  color = 'b', linestyle = '-')
    axis_13.plot(pi_c,eta_T_i[:,1], label='tau_t ideal: 0.5', color = 'm', linestyle = '--') 
    axis_13.plot(pi_c,eta_T_r[:,1], label='tau_t real: 0.5',  color = 'm', linestyle = '-') 
    axis_13.plot(pi_c,eta_T_i[:,2], label='tau_t ideal: 0.6', color = 'c', linestyle = '--') 
    axis_13.plot(pi_c,eta_T_r[:,2], label='tau_t real: 0.6',  color = 'c', linestyle = '-') 
    axis_13.plot(pi_c,eta_T_i[:,3], label='tau_t ideal: 0.7', color = 'k', linestyle = '--') 
    axis_13.plot(pi_c,eta_T_r[:,3], label='tau_t real: 0.7',  color = 'k', linestyle = '-') 
    axis_13.plot(pi_c,eta_T_i[:,4], label='tau_t ideal: 0.8', color = 'g', linestyle = '--') 
    axis_13.plot(pi_c,eta_T_r[:,4], label='tau_t real: 0.8',  color = 'g', linestyle = '-')     
    axis_13.set_xlabel('pi_c [-]')
    axis_13.set_ylabel('eta_T [-]')
    axis_13.legend()
    axis_13.set_xlim(0, 40)
    axis_13.set_ylim(0.20, 0.70)
    fig13.tight_layout()  
    
    fig14    =  plt.figure('eta_P ideal vs real')
    fig14.set_size_inches(7, 6)
    axis_14 = fig14.add_subplot(1,1,1)
    axis_14.plot(pi_c,eta_P_i[:,0], label='tau_t ideal: 0.4', color = 'b', linestyle = '--')
    axis_14.plot(pi_c,eta_P_r[:,0], label='tau_t real: 0.4',  color = 'b', linestyle = '-')
    axis_14.plot(pi_c,eta_P_i[:,1], label='tau_t ideal: 0.5', color = 'm', linestyle = '--') 
    axis_14.plot(pi_c,eta_P_r[:,1], label='tau_t real: 0.5',  color = 'm', linestyle = '-') 
    axis_14.plot(pi_c,eta_P_i[:,2], label='tau_t ideal: 0.6', color = 'c', linestyle = '--') 
    axis_14.plot(pi_c,eta_P_r[:,2], label='tau_t real: 0.6',  color = 'c', linestyle = '-') 
    axis_14.plot(pi_c,eta_P_i[:,3], label='tau_t ideal: 0.7', color = 'k', linestyle = '--') 
    axis_14.plot(pi_c,eta_P_r[:,3], label='tau_t real: 0.7',  color = 'k', linestyle = '-') 
    axis_14.plot(pi_c,eta_P_i[:,4], label='tau_t ideal: 0.8', color = 'g', linestyle = '--') 
    axis_14.plot(pi_c,eta_P_r[:,4], label='tau_t real: 0.8',  color = 'g', linestyle = '-')     
    axis_14.set_xlabel('pi_c [-]')
    axis_14.set_ylabel('eta_P [-]')
    axis_14.legend()
    axis_14.set_xlim(0, 40)
    axis_14.set_ylim(0.20, 0.90)
    fig14.tight_layout()    
    
    return

def turboprop_model_test():
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Complete set of equations for Parametric Cycle Analysis of Real Case
    #  Reference: Elements or gas turbine propulsion/Jack D. Mallingly, pag 433
    # ----------------------------------------------------------------------------------------------------------------------    
    
    #F_c_m_dot_0          = (a0/g_c)*((1 + f)*V9_a0 - M0 + (1 + f)*(Rt/Rc)*((T9/T0)/(V9/a0))*((1 - (P0/P9))/gamma_c))
    #Cc                   = (V0/(cp*T0))*((1 + f)*V9_a0 - M0 + (1 + f)*(Rt/Rc)*((T9/T0)/(V9/a0))*((1 - (P0/P9))/gamma_c))
    #Cc                   = (gamma_c - 1)*M0*((1 + f)*V9_a0 - M0 + (1 + f)*(Rt/Rc)*((T9/T0)/(V9/a0))*((1 - (P0/P9))/gamma_c))
    #(V9_a0)**2           = ((gamma_t*R_t*T_9)/(gamma_c*R_c*T0))*M9**2
    #M9**2                = (2/(gamma_t - 1))*((P_t9/P_9)**((gamma_t - 1)/gamma_t) - 1)
    #Pt9/P9               = (P0/P9)*pi_r*pi_d*pi_c*pi_b*pi_tH*pi_tL*pi_n
    #For unchoked flow:
    # Pt9/P0               <= ((gamma_t +1)/2)**(gamma_t/(gamma_t - 1))
    #For choked flow:
    # Pt9/P0               > ((gamma_t +1)/2)**(gamma_t/(gamma_t - 1))
    # M9                   = 1
    # Pt9/P9               = ((gamma_t +1)/2)**(gamma_t/(gamma_t - 1))
    # P0/P9                = ((Pt9/P9)/(Pt9/P0))
    #T9/T0                = (Tt9/T0)/((Pt9/P9)**((gamma_t - 1)/gamma_t))
    #Tt9/T0               = (Tt4/T0)*tau_tH*tau_tL
    #V9/a0                = np.sqrt(((2*tau_lambda*tau_tH*tau_tL)/(gamma_c - 1))*(1 - (Pt9/P9)**((gamma_t - 1)/gamma_t))))
    #f                    = (tau_lambda - tau_r*tau_c_r)/(eta_b*h_PR/(cp*T0) - tau_lambda)
    #tau_tH               = 1 - (1/(eta_mH*(1 + f)))*(tau_r/tau_lambda)*(tau_c_r - 1)
    #W_dot_prop           = eta_g*eta_mL*m_dot_45*cpt*(Tt45 - Tt5)
    #C_prop               = (eta_prop*W_dot_prop)/(m_dot_0*cp*T0) 
    #C_prop               = eta_prop*eta_g*eta_mL*(1 + f)*tau_lambda*tau_tH*(1 - tau_tL)
    #W_dot_m_dot_0        = C_tot*cpc*T0
    #F_m_dot_0            = (C_tot*cpc*T0)/V0
    #S                    = f/(F/m_dot_0)    
    #S_P                  = f/(C_tot*cpc*T0)
    #eta_T                = (m_dot_0*cpc*T0*C_tot)/(m_dot_f*h_PR)
    #eta_T                = C_tot/(f*h_PR/(cpc*T0))
    #eta_P                = (m_dot_0*cpc*T0*C_tot)/(W_dot_prop*(m_dot_9*V_9**2 - m_dot_0*V_0**2)/(2*g_c)) 
    #eta_P                = C_tot/(C_tot/eta_prop + ((gamma - 1)/2)*((1 + f)*(V9/a0)**2 - M0**2))    
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Input values for real turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    # 1) Real engine
    # 2) Eng library James Aircraft for Validation
    # 3) Integration into RCAIDE
    # 4) Get ATR 72 / Otter Engine data
    
    M0_r                  = 0.8                       # [-]
    T0_r                  = 240                       # [K]
    gamma_c               = 1.4                       # [-]
    cp_c                  = 1004                      # [J/(kg*K)]
    gamma_t               = 1.35                      # [-]
    cp_t                  = 1108                      # [J/(kg*K)]    
    h_PR                  = 42800000                  # [J/kg]
    e_c                   = 0.9                       # [-]
    e_tH                  = 0.89                      # [-]
    e_tL                  = 0.91                      # [-]
    pi_d                  = 0.98                      # [-]
    pi_b                  = 0.96                      # [-] 
    pi_n                  = 0.99                      # [-]
    eta_b                 = 0.99                      # [-]
    eta_mH                = 0.99                      # [-] 
    eta_mL                = 0.99                      # [-] 
    eta_prop              = 0.83                      # [-]
    eta_g                 = 0.99                      # [-]    
    Tt4                   = 1370                      # [K]
    pi_c                  = np.linspace(2.5,40,100)   # [-]
    tau_t                 = [0.4, 0.5, 0.6, 0.7, 0.8] # [-]
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Compiutation for real turboprop
    # ----------------------------------------------------------------------------------------------------------------------
                      
    tau_c_r               = np.zeros((len(pi_c), len(tau_t)))    
    eta_c                 = np.zeros_like(tau_c_r)
    f_r                   = np.zeros_like(tau_c_r)
    tau_tH_r              = np.zeros_like(tau_c_r)
    pi_tH                 = np.zeros_like(tau_c_r)
    eta_tH                = np.zeros_like(tau_c_r)
    tau_tL_r              = np.zeros_like(tau_c_r)
    pi_tL                 = np.zeros_like(tau_c_r)
    eta_tL                = np.zeros_like(tau_c_r) 
    Pt9_P0                = np.zeros_like(tau_c_r)
    M9                    = np.zeros_like(tau_c_r)
    Pt9_P9                = np.zeros_like(tau_c_r)
    P0_P9                 = np.zeros_like(tau_c_r)
    Tt9_T0                = np.zeros_like(tau_c_r)
    T9_T0                 = np.zeros_like(tau_c_r)
    V9_a0_r               = np.zeros_like(tau_c_r)
    C_prop_r              = np.zeros_like(tau_c_r)
    Cc_r                  = np.zeros_like(tau_c_r)  
    C_tot_r               = np.zeros_like(tau_c_r)
    F_mdot0_r             = np.zeros_like(tau_c_r)
    S_r                   = np.zeros_like(tau_c_r)
    W_dot_mdot0           = np.zeros_like(tau_c_r)
    S_P                   = np.zeros_like(tau_c_r)
    eta_T_r               = np.zeros_like(tau_c_r)
    eta_P_r               = np.zeros_like(tau_c_r)
    
    for i in range(len(pi_c)):
        for j in range(len(tau_t)):
            
            R_c                  = ((gamma_c - 1)/gamma_c)*cp_c
            R_t                  = ((gamma_t - 1)/gamma_t)*cp_t
            a0                   = np.sqrt(gamma_c*R_c*T0_r)
            V0                   = a0*M0_r
            tau_r                = 1 + ((gamma_c - 1)/2)*M0_r**2
            pi_r                 = tau_r**(gamma_c/(gamma_c - 1))
            if M0_r <= 1:
                eta_r            = 1
            else:
                eta_r            = 1 - 0.0075(M0_r - 1)**1.35
            tau_lambda           = (cp_t*Tt4)/(cp_c*T0_r)
            tau_c_r[i,j]         = pi_c[i]**((gamma_c - 1)/(gamma_c*e_c))
            eta_c[i,j]           = (pi_c[i]**((gamma_c - 1)/(gamma_c)) - 1)/(tau_c_r[i,j] - 1)
            f_r[i,j]             = (tau_lambda - tau_r*tau_c_r[i,j])/((eta_b*h_PR)/(cp_c*T0_r) - tau_lambda)
            tau_tH_r[i,j]        = 1 - ((tau_r*(tau_c_r[i,j] - 1))/(eta_mH*(1 + f_r[i,j])*tau_lambda))
            pi_tH[i,j]           = tau_tH_r[i,j]**(gamma_t/((gamma_t - 1)*e_tH))
            eta_tH[i,j]          = (1 - tau_tH_r[i,j])/(1 - tau_tH_r[i,j]**(1/e_tH))
            tau_tL_r[i,j]        = tau_t[j]/tau_tH_r[i,j]
            pi_tL[i,j]           = tau_tL_r[i,j]**(gamma_t/((gamma_t - 1)*e_tL))
            eta_tL[i,j]          = (1 - tau_tL_r[i,j])/(1 - tau_tL_r[i,j]**(1/e_tL))
            Pt9_P0[i,j]          = pi_r*pi_d*pi_c[i]*pi_b*pi_tH[i,j]*pi_tL[i,j]*pi_n
            
            if Pt9_P0[i,j] > ((gamma_t + 1)/2)**(gamma_t/(gamma_t - 1)):
                M9[i,j]          = 1
                Pt9_P9[i,j]      = ((gamma_t + 1)/2)**(gamma_t/(gamma_t - 1))
                P0_P9[i,j]       = ((Pt9_P9[i,j])/(Pt9_P0[i,j]))
            else:
                Pt9_P9[i,j]      = Pt9_P0[i,j]
                M9[i,j]          = np.sqrt((2/(gamma_t - 1))*((Pt9_P0[i,j])**((gamma_t - 1)/gamma_t)) - 1)
                P0_P9[i,j]       = 1
            Tt9_T0[i,j]          = (Tt4/T0_r)*tau_tH_r[i,j]*tau_tL_r[i,j]
            T9_T0[i,j]           = (Tt9_T0[i,j])/((Pt9_P9[i,j])**((gamma_t - 1)/gamma_t))
            V9_a0_r[i,j]         = np.sqrt(((2*tau_lambda*tau_tH_r[i,j]*tau_tL_r[i,j])/(gamma_c - 1))*(1 - Pt9_P9[i,j]**(-(gamma_t - 1)/(gamma_t))))
            C_prop_r[i,j]        = eta_prop*eta_g*eta_mL*(1 + f_r[i,j])*tau_lambda*tau_tH_r[i,j]*(1 - tau_tL_r[i,j])
            Cc_r[i,j]            = (gamma_c - 1)*M0_r*((1 + f_r[i,j])*V9_a0_r[i,j] - M0_r + (1 + f_r[i,j])*(R_t/R_c)*((T9_T0[i,j])/(V9_a0_r[i,j]))*((1 - (P0_P9[i,j]))/gamma_c))
            C_tot_r[i,j]         = Cc_r[i,j] + C_prop_r[i,j]
            F_mdot0_r[i,j]       = (C_tot_r[i,j]*cp_c*T0_r)/(V0)
            S_r[i,j]             = (f_r[i,j]/(F_mdot0_r[i,j]))*1000**2
            W_dot_mdot0[i,j]     = C_tot_r[i,j]*cp_c*T0_r
            S_P[i,j]             = (f_r[i,j]/(C_tot_r[i,j]*cp_c*T0_r))
            eta_T_r[i,j]         = C_tot_r[i,j]/((f_r[i,j]*h_PR)/(cp_c*T0_r))
            eta_P_r[i,j]         = C_tot_r[i,j]/((C_prop_r[i,j]/eta_prop) + ((gamma_c - 1)/2)*((1 + f_r[i,j])*(V9_a0_r[i,j])**2 - M0_r**2))  
              
    # ----------------------------------------------------------------------------------------------------------------------
    #  Complete set of equations for Parametric Cycle Analysis of Real Case
    #  Reference: Elements or gas turbine propulsion/Jack D. Mallingly, pag 560
    # ----------------------------------------------------------------------------------------------------------------------            
    
    #mdot_0                       = (P_0*pi_r*pi_d*pi_c)/()                         
        
    # ----------------------------------------------------------------------------------------------------------------------
    #  Results for ideal turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    #fig1    =  plt.figure('F_mdot0_i')
    #fig1.set_size_inches(7, 6)
    #axis_1 = fig1.add_subplot(1,1,1)
    #axis_1.plot(pi_c,F_mdot0_i[:,0], label='tau_t: 0.4')
    #axis_1.plot(pi_c,F_mdot0_i[:,1], label='tau_t: 0.5') 
    #axis_1.plot(pi_c,F_mdot0_i[:,2], label='tau_t: 0.6') 
    #axis_1.plot(pi_c,F_mdot0_i[:,3], label='tau_t: 0.7') 
    #axis_1.plot(pi_c,F_mdot0_i[:,4], label='tau_t: 0.8') 
    #axis_1.set_xlabel('pi_c [-]')
    #axis_1.set_ylabel('F_m_dot_0 [N/(kg/sec)]')
    #axis_1.legend()
    #axis_1.set_ylim(800, 1600)
    #fig1.tight_layout()

    return

if __name__ == '__main__': 
    main()
    plt.show()