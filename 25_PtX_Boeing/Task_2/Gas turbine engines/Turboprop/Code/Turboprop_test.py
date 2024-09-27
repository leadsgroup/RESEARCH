import numpy as np
import matplotlib.pyplot as plt

def main():
    
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Input values for ideal turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    M0         = 0.8                       # [-]
    T0         = 240                       # [K]
    gamma      = 1.4                       # [-]
    cp         = 1004                      # [J/(kg*K)]
    h_PR       = 42800000                  # [J/kg]
    Tt4        = 1370                      # [K]
    pi_c       = np.linspace(2.5,40,100)   # [-]
    tau_t      = [0.4, 0.5, 0.6, 0.7, 0.8] # [-]
    eta_prop   = 0.83                      # [-]
    g_c        = 9.81                      # [m/s**2]
    
    # ----------------------------------------------------------------------------------------------------------------------
    # Compiutation for ideal turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    F_mdot0   = np.zeros((len(pi_c), len(tau_t)))
    
    tau_c     = np.zeros_like(F_mdot0)
    f         = np.zeros_like(F_mdot0)
    tau_tH    = np.zeros_like(F_mdot0)
    tau_tL    = np.zeros_like(F_mdot0)
    V9_a9     = np.zeros_like(F_mdot0)
    Cc        = np.zeros_like(F_mdot0)
    C_prop    = np.zeros_like(F_mdot0)
    C_tot     = np.zeros_like(F_mdot0)
    S         = np.zeros_like(F_mdot0)
    eta_T     = np.zeros_like(F_mdot0)
    eta_O     = np.zeros_like(F_mdot0)
    eta_P     = np.zeros_like(F_mdot0)
    
    for i in range(len(pi_c)):
        for j in range(len(tau_t)):
            
            R                    = ((gamma - 1)/gamma)*cp
            a0                   = np.sqrt(gamma*R*T0)
            tau_r                = 1 + ((gamma - 1)/2)*M0**2
            tau_lambda           = Tt4/T0
            tau_c[i,j]           = pi_c[i]**((gamma - 1)/gamma)
            f[i,j]               = ((cp*T0)/h_PR)*(tau_lambda - tau_r*tau_c[i,j])
            tau_tH[i,j]          = 1 - (tau_r/tau_lambda)*(tau_c[i,j] - 1)
            #tau_t_opt[i,j]       = (tau_lambda/tau_r*tau_c[i,j]) + (((gamma-1)/2)*M0**2)/(tau_lambda * eta_prop**2)
            tau_tL[i,j]          = tau_t[j]/tau_tH[i,j]
            V9_a9[i,j]           = np.sqrt((2/(gamma - 1))*(tau_lambda*tau_t[j] - tau_lambda/(tau_r*tau_c[i,j])))
            Cc[i,j]              = (gamma - 1)*M0*(V9_a9[i,j] - M0)
            C_prop[i,j]          = eta_prop*tau_lambda*tau_tH[i,j]*(1 - tau_tL[i,j])
            C_tot[i,j]           = Cc[i,j] + C_prop[i,j]
            F_mdot0[i,j]         = (C_tot[i,j]*cp*T0)/(M0*a0)
            S[i,j]               = (f[i,j]/(F_mdot0[i,j]))*1000**2
            eta_T[i,j]           = 1 - 1/(tau_r*tau_c[i,j])
            eta_O[i,j]           = C_tot[i,j]/(tau_lambda - tau_r*tau_c[i,j])
            eta_P[i,j]           = eta_O[i,j]/eta_T[i,j]
    
    # ----------------------------------------------------------------------------------------------------------------------
    #  Results for ideal turboprop
    # ----------------------------------------------------------------------------------------------------------------------
    
    fig1    =  plt.figure('F_mdot0')
    fig1.set_size_inches(7, 6)
    axis_1 = fig1.add_subplot(1,1,1)
    axis_1.plot(pi_c,F_mdot0[:,0], label='tau_t: 0.4')
    axis_1.plot(pi_c,F_mdot0[:,1], label='tau_t: 0.5') 
    axis_1.plot(pi_c,F_mdot0[:,2], label='tau_t: 0.6') 
    axis_1.plot(pi_c,F_mdot0[:,3], label='tau_t: 0.7') 
    axis_1.plot(pi_c,F_mdot0[:,4], label='tau_t: 0.8')     
    axis_1.set_xlabel('pi_c [-]')
    axis_1.set_ylabel('F_m_dot_0 [N/(kg/sec)]')
    axis_1.legend()
    axis_1.set_ylim(800, 1600)
    fig1.tight_layout()
    
    fig2    =  plt.figure('S')
    fig2.set_size_inches(7, 6)
    axis_2 = fig2.add_subplot(1,1,1)
    axis_2.plot(pi_c,S[:,0], label='tau_t: 0.4')
    axis_2.plot(pi_c,S[:,1], label='tau_t: 0.5') 
    axis_2.plot(pi_c,S[:,2], label='tau_t: 0.6') 
    axis_2.plot(pi_c,S[:,3], label='tau_t: 0.7') 
    axis_2.plot(pi_c,S[:,4], label='tau_t: 0.8')     
    axis_2.set_xlabel('pi_c [-]')
    axis_2.set_ylabel('S [mg/(N*sec)]')
    axis_2.legend()
    axis_2.set_ylim(10, 26)
    fig2.tight_layout()    
    
    return

if __name__ == '__main__': 
    main()
    plt.show()