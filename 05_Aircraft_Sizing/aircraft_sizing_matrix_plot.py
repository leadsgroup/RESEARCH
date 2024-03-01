# imports  
import numpy as np           
import matplotlib.pyplot as plt            

def main():
    Performance_Requirements()    
    return
    
def Performance_Requirements():
    
    # -------------------------------------------------------------------------------------    
    # User input variables 
    # -------------------------------------------------------------------------------------
    
    MTOW       = 160000                                                    #[lbf]
    n          = 2
    V_stall    = 200                                                       #[ft/s]
    V_2        = 1.2*V_stall                                               #[ft/s]    
    V_LO       = 250                                                       #[ft/s]
    V_LND      = 270                                                       #[ft/s]
    V_cr       = 400                                                       #[ft/s]
    V_turn     = 350                                                       #[ft/s]
    C_L_max    = 2                                                         #[-]
    C_L_TO     = 1.5                                                       #[-]
    C_L_LND    = 1.8                                                       #[-]
    C_L_in_OEI = 1.8                                                       #[-]
    C_L_tr_OEI = 1.4                                                       #[-]
    C_L_es_OEI = 1.2                                                       #[-]
    C_L_cr_OEI = 1.0                                                       #[-]
    C_L_ab_AEO = 1.3                                                       #[-]
    C_L_ab_OEI = 1.2                                                       #[-]
    rho_0      = 2.3769e-3                                                 #[slugs/ft**3]
    rho_stall  = 2.3769e-3                                                 #[slugs/ft**3]
    rho_TO     = 2.3769e-3                                                 #[slugs/ft**3]
    rho_LND    = 2.3769e-3                                                 #[slugs/ft**3]
    rho_in_OEI = 2.3769e-3                                                 #[slugs/ft**3]
    rho_tr_OEI = 2.2e-3                                                    #[slugs/ft**3]
    rho_es_OEI = 2.1e-3                                                    #[slugs/ft**3]
    rho_cr_OEI = 2.1e-3                                                    #[slugs/ft**3]
    rho_ab_AEO = 2.3769e-3                                                 #[slugs/ft**3]
    rho_ab_OEI = 2.3769e-3                                                 #[slugs/ft**3]
    rho_cr     = 2.0e-3                                                    #[slugs/ft**3]
    rho_turn   = 2.0e-3                                                    #[slugs/ft**3]
    h_obs      = 35                                                        #[ft] from FAR 25
    h_tr       = 20                                                        #[ft]
    R          = 50                                                        #[ft]
    gamma      = 0.15                                                      #[rad]
    gamma_in_OEI= 0.075398224                                              #[rad]
    gamma_tr_OEI= 0                                                        #[rad]
    gamma_es_OEI= 0.150796447                                              #[rad]
    gamma_cr_OEI= 0.075398224                                              #[rad]
    gamma_ab_AEO= 0.20106193                                               #[rad]
    gamma_ab_OEI= 0.131946891                                               #[rad]
    T_TO       = 40000                                                     #[lbf]
    C_D_0      = 0.0179                                                    #[-]
    k          = 0.0427                                                    #[-]
    S          = 1341                                                      #[ft**2]  
    mu         = 0.02                                                      #[-]
    N          = 2                                                         #[-]
        
    
    # -------------------------------------------------------------------------------------
    # Sizing calculations
    # ------------------------------------------------------------------------------------- 
    
    # Wing loading vector
    
    WS = np.arange(3,201,1)
    
    
    # Stall
    
    WS_stall = 0.5*rho_stall*(V_stall**2)*C_L_max                       #[psf]
    
    
    # Take off
    
    C_D_TO = C_D_0 + k*C_L_TO**2
    D_TO = 0.5*rho_TO*((V_LO/np.sqrt(2))**2)*S*C_D_TO
    L_TO = 0.5*rho_TO*((V_LO/np.sqrt(2))**2)*S*C_L_TO
    a = (32.174/MTOW)*(T_TO - D_TO - mu*(MTOW - L_TO))
    S_G = (V_LO**2)/(2*a)                                               #[ft]
    
    if h_obs > h_tr:      
        S_TR = R*np.sin(gamma)                                          #[ft]
        S_CL = (h_obs - h_tr)/(np.tan(gamma))                           #[ft]
        S_TO = S_G + S_TR + S_CL                                        #[ft]       
    else:       
        S_TR = np.sqrt(2*R*h_obs - h_obs)                               #[ft]
        S_TO = S_G + S_TR                                               #[ft]
        
    TOP_25 = S_TO/37.5                                                  #[lb/ft**2]
    TW_TO  = np.zeros(len(WS))                                          #[-]
    for i in range(len(WS)):
        TW_TO[i]  = WS[i]*(1/((rho_TO/rho_0)*C_L_TO*TOP_25))            #[-]
        
        
    # Landing
    
    WS_LND = 0.5*rho_LND*(V_LND**2)*C_L_LND                             #[psf] 
    
    
    # Initial Climb OEI
    
    C_D_in_OEI = C_D_0 + k*C_L_in_OEI**2
    D_in_OEI = 0.5*rho_in_OEI*((V_2)**2)*S*C_D_in_OEI
    L_in_OEI = 0.5*rho_in_OEI*((V_2)**2)*S*C_L_in_OEI
    TW_in_OEI = (n/(n - 1))*(np.sin(gamma_in_OEI) + (D_in_OEI/L_in_OEI))
    
    
    # Transitional Climb OEI
    
    C_D_tr_OEI = C_D_0 + k*C_L_tr_OEI**2
    D_tr_OEI = 0.5*rho_tr_OEI*((1.1*V_stall)**2)*S*C_D_tr_OEI
    L_tr_OEI = 0.5*rho_tr_OEI*((1.1*V_stall)**2)*S*C_L_tr_OEI
    TW_tr_OEI = (n/(n - 1))*(np.sin(gamma_tr_OEI) + (D_tr_OEI/L_tr_OEI))
    
    
    # Established Climb OEI
    
    C_D_es_OEI = C_D_0 + k*C_L_es_OEI**2
    D_es_OEI = 0.5*rho_es_OEI*((1.2*V_stall)**2)*S*C_D_es_OEI
    L_es_OEI = 0.5*rho_es_OEI*((1.2*V_stall)**2)*S*C_L_es_OEI
    TW_es_OEI = (n/(n - 1))*(np.sin(gamma_es_OEI) + (D_es_OEI/L_es_OEI))
    
    
    # Cruise Climb OEI
    
    C_D_cr_OEI = C_D_0 + k*C_L_cr_OEI**2
    D_cr_OEI = 0.5*rho_cr_OEI*((1.2*V_stall)**2)*S*C_D_cr_OEI
    L_cr_OEI = 0.5*rho_cr_OEI*((1.2*V_stall)**2)*S*C_L_cr_OEI
    TW_cr_OEI = (n/(n - 1))*(np.sin(gamma_cr_OEI) + (D_cr_OEI/L_cr_OEI))
    
    
    # Aborted landing AEO
    
    C_D_ab_AEO = C_D_0 + k*C_L_ab_AEO**2
    D_ab_AEO = 0.5*rho_ab_AEO*((1.3*V_stall)**2)*S*C_D_ab_AEO
    L_ab_AEO = 0.5*rho_ab_AEO*((1.3*V_stall)**2)*S*C_L_ab_AEO
    TW_ab_AEO = (np.sin(gamma_ab_AEO) + (D_ab_AEO/L_ab_AEO))   
    
    
    # Aborted landing OEI
    
    C_D_ab_OEI = C_D_0 + k*C_L_ab_OEI**2
    D_ab_OEI = 0.5*rho_ab_OEI*((1.5*V_stall)**2)*S*C_D_ab_OEI
    L_ab_OEI = 0.5*rho_ab_OEI*((1.5*V_stall)**2)*S*C_L_ab_OEI
    TW_ab_OEI = (n/(n - 1))*(np.sin(gamma_ab_OEI) + (D_ab_OEI/L_ab_OEI))  
    
    
    # Maximum Cruise Speed
    
    TW_cr  = np.zeros(len(WS))                                          #[-]
    for i in range(len(WS)):
        TW_cr[i]  = (0.5*rho_cr*(V_cr**2)*C_D_0)/WS[i] + (k/(0.5*rho_cr*(V_cr**2)))*WS[i]           #[-]   
        

    # Coordinated turn
    
    TW_turn  = np.zeros(len(WS))                                          #[-]
    for i in range(len(WS)):
        TW_turn[i]  = (0.5*rho_turn*(V_turn**2)*C_D_0)/WS[i] + ((k*N**2)/(0.5*rho_turn*(V_turn**2)))*WS[i]           #[-]   
    
    
    # -------------------------------------------------------------------------------------
    # Sizing Matrix Plot
    # -------------------------------------------------------------------------------------     
    
    plt.figure()
    plt.axvline(WS_stall, color='m', linewidth=2.5)
    plt.axvline(WS_LND, color='c', linewidth=2.5)
    plt.plot(WS, TW_TO, color='b', linewidth=2.5)
    plt.axhline(TW_in_OEI, color='g', linewidth=2.5)
    plt.axhline(TW_tr_OEI, color='g', linestyle='--', linewidth=2.5)
    plt.axhline(TW_es_OEI, color='g', linestyle=':', linewidth=2.5)
    plt.axhline(TW_cr_OEI, color='g', linestyle='-.', linewidth=2.5)
    plt.axhline(TW_ab_AEO, color='k', linestyle=':', linewidth=2.5)
    plt.axhline(TW_ab_OEI, color='k', linestyle='--', linewidth=2.5)    
    plt.plot(WS, TW_cr, color='y', linewidth=2.5)
    plt.plot(WS, TW_turn, color='r', linewidth=2.5)
    plt.legend(['Stall','Landing','Take off','Initial climb OEI',
                'Transitional climb OEI','Established climb OEI',
                'Cruise climb OEI','Aborted landing AEO',
                'Aborted landing OEI', 'Max cruise speed',
                'Coordinated turn'],loc='upper right', fontsize=12)
    plt.box(True)
    plt.title('SMP', fontsize=20)
    plt.xlabel('$(W/S)_{TO}$ [psf]', fontsize=20)
    plt.ylabel('$(T/W)_{TO}$ [-]', fontsize=20)
    plt.grid(True)    
     
if __name__ == '__main__':
    main()
    plt.show()