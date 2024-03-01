# imports  
import numpy as np           
import matplotlib.pyplot as plt            

def main():
    Performance_Requirements()    
    return
    
def Performance_Requirements():
    
    # -------------------------------------------------------------------------------------    
    # user input variables 
    # -------------------------------------------------------------------------------------
    
    #MTOW    = 
    #V_stall = 
    #C_L_max = 
    #rho     = 
    
        
    
    # Lift coefficient
    cL_max = 1.7
    cL_max_TO = 4
    cL_max_LND = 6
    cL = np.linspace(6.3, 6.3, 100)                    # for the drag polar
    
    # cD0 estimation: Roskam
    W_MTO_h = [2200, 2500, 2028.25]                    # only homebuilt
    W_MTO_s = [2500, 3000, 2500, 2350.13, 2450, 2548]  # only single engine
    W_MTO = np.concatenate((W_MTO_s, W_MTO_h))
    c, d = 1.0892, 0.5147                              # single
    c_h, d_h = 1.2362, 0.4319                          # homebuilt
    a, b = -2.2614, 1                                  # interpolation
    L_Swet_s = [c + d * np.log10(W) for W in W_MTO_s]
    L_Swet_h = [c_h + d_h * np.log10(W) for W in W_MTO_h]
    L_Swet = np.concatenate((L_Swet_s, L_Swet_h))
    
    # Parassite area and cD0
    S_h = [167, 180, 134.2]                            # homebuilt
    S_s = [165.6, 231, 131, 131, 174, 170]             # single
    S = np.concatenate((S_s, S_h))
    
    L_f = a + b * L_Swet
    f = 10 ** L_f
    
    cD_0 = 0.02                                        # Initial gues
    
    AR_m = 8.3                                         
    
    e_CLEAN = 1.78 * (1 - 0.045 * AR_m ** 0.68) - 0.64
    e_TO = 0.78                                        # from Roskam
    e_LND = 0.72                                       # from Roskam
    k_CLEAN = 1 / (np.pi * AR_m * e_CLEAN)
    k_TO = 1 / (np.pi * AR_m * e_TO)
    k_LND = 1 / (np.pi * AR_m * e_LND)
    
    # After simulations and verifications
    cD_0_TO = 0.085
    cD_0_LND = 0.105
    cD_0_GEAR = cD_0 + 0.028
    
    cD_CLEAN = cD_0 + k_CLEAN * cL ** 2
    cD_TO = cD_0_TO + k_TO * cL ** 2
    cD_LND = cD_0_LND + k_LND * cL ** 2
        
    # SI unit
    theta_ref = 298.15                                 # by considering ISA + 18 °C
    g = 9.81
    R = 287.05
    lamb = -6.5 * 1e-3
    rho0 = 1.225
    h_1 = 0 * 0.3048
    theta1 = theta_ref + lamb * h_1
    
    sigma1 = (theta1 / theta_ref) ** (-(1 + g / (R * lamb)))
    rho1 = sigma1 * rho0
    
    # sigma 2: 5000 ft
    h_2 = 5000 * 0.3048  # m
    theta2 = theta_ref + lamb * h_2  # K
    
    sigma2 = (theta2 / theta_ref) ** (-(1 + g / (R * lamb)))
    rho2 = sigma2 * rho0
    
    # sigma 3: 10000 ft
    h_3 = 10000 * 0.3048  # m
    theta_ref0 = 288.15  # K
    sigma3 = (1 + lamb * h_3 / theta_ref0) ** (-(1 + g / (R * lamb)))
    rho3 = sigma3 * rho0
    
    # sigma 4: 15000 ft
    h_4 = 15000 * 0.3048  # m
    sigma4 = (1 + lamb * h_4 / theta_ref0) ** (-(1 + g / (R * lamb)))
    rho4 = sigma4 * rho0
    
    # sigma 5: 6000 ft for the second part of the climb
    h_5 = 6000 * 0.3048  # m
    sigma5 = (1 + lamb * h_5 / theta_ref0) ** (-(1 + g / (R * lamb)))
    rho5 = sigma5 * rho0   


    

    # -------------------------------------------------------------------------------------
    # Sizing calculations
    # ------------------------------------------------------------------------------------- 
    
    # Stall
    
    v_s = 28  # kcas
    vs = v_s / 1.94384  # m/s
    
    # TO phase parameters
    S_TO = 300  # ft
    S_TOG = S_TO / 1.66
    mu = [0.025, 0.05]  # concrete and short grass
    PDL = 3             # [hp/ft^2] Single Prop case from Roskam
    K1 = 0.0376
    K2_1 = 5.75 * (sigma1 * PDL) ** (1 / 3)
    K2_2 = 5.75 * (sigma2 * PDL) ** (1 / 3)
    rho1_1 = rho1 * 0.00194  # slug/ft^3
    rho2_1 = rho2 * 0.00194  # slug/ft^3
    
    # TO phase functions
    def y1_1(x1, mu):
        return ((1 / (S_TO * rho1_1 * cL_max_TO * K2_1)) * (
                K1 * 1.66 * x1 + S_TO * rho1_1 * cL_max_TO * mu + S_TO * rho1_1 * 0.72 * cD_0_TO)) ** (-1)
    
    def y1_2(x1, mu):
        return ((1 / (S_TO * rho2_1 * cL_max_TO * K2_2)) * (
                K1 * 1.66 * x1 + S_TO * rho2_1 * cL_max_TO * mu + S_TO * rho2_1 * 0.72 * cD_0_TO)) ** (-1)
    
    x1 = np.arange(0, 40, 0.1)  
    
    # Stall phase calculation
    x2 = ((0.5 * (rho0) * (vs ** 2) * cL_max_LND) / 0.93) * 0.02087825  # 0.975
    
    # Constants
    RC = 1500 * 0.00508  # m/s
    RC_2 = 750 * 0.00508  # m/s
    eta_p = 0.8
    
    # Climb phase equations
    v_Y = lambda x, rho: np.sqrt((2 / rho) * x * np.sqrt(k_CLEAN / (3 * cD_0)))
    v_Y_TO = lambda x, rho: np.sqrt((2 / rho) * x * np.sqrt(k_TO / (3 * cD_0_TO)))
    v_Y_LND = lambda x, rho: np.sqrt((2 / rho) * x * np.sqrt(k_LND / (3 * cD_0_LND)))
    
    # Functions for different scenarios
    y3_1 = lambda x3: 167.58 * (((((1 / eta_p) * (0.5 * rho1 * (v_Y(x3, rho1) ** 3) * cD_0 / x3 +
                                                  2 * k_CLEAN * x3 / (rho1 * v_Y(x3, rho1)) + RC))) ** (-1))) / 0.9997
    y3_1_2 = lambda x3: 167.58 * (((((1 / eta_p) * (0.5 * rho5 * (v_Y(x3, rho5) ** 3) * cD_0 / x3 +
                                                    2 * k_CLEAN * x3 / (rho5 * v_Y(x3, rho1)) + RC_2))) ** (-1))) / 0.9997
    y3_1_c = lambda x3: 167.58 * (((((1 / eta_p) * (0.5 * rho1 * (v_Y(x3, rho1) ** 3) * cD_0 / x3 + 
                                                    2 * k_CLEAN * x3 / (rho1 * v_Y(x3, rho1)) + np.sin(20 / 100) * v_Y(x3, rho1)))) ** (-1))) / 0.9997
    y3_2_OEI = lambda x3: 167.58 * (((((1 / eta_p) * (0.5 * rho2 * (v_Y(x3, rho2) ** 3) * cD_0 / x3 +
                                                      2 * k_CLEAN * x3 / (rho2 * v_Y(x3, rho2)) + np.sin(1.5 / 100) *
                                                    v_Y(x3, rho2)))) ** (-1)))/0.9997
    y3_1_b = lambda x3: 167.58 * (((((1 / eta_p) * (0.5 * rho1 * (v_Y_LND(x3, rho1) ** 3) * cD_0_LND / x3 +
                                                    2 * k_LND * x3 / (rho1 * v_Y_LND(x3, rho1)) + np.sin(3.3 / 100) *
                                                  v_Y_LND(x3, rho1)))) ** (-1))) / 0.9997
    
    x3 = np.arange(0, 40 / 0.02087825, 0.1)
    
    # Service ceiling
    RC_sc = 0.508  # m/s at the condition under examination
    
    y4 = lambda x4: 167.58 * (((((1 / eta_p) * (0.5 * rho4 * v_Y(x4, rho4) ** 3 * cD_0 / x4 +
                                                2 * k_CLEAN * x4 / (rho4 * v_Y(x4, rho4)) + RC_sc))) ** (-1))) / 0.97
    x4 = np.arange(0, 40 / 0.02087825, 0.1)
    
    # Cruise phase
    v_cruise = 170 / 1.94384  # m/s
    v_min_cruise = 150 / 1.94384  # m/s
    
    y5 = lambda x5: 167.58 * (((((0.5 * rho3 * v_cruise ** 3 * cD_0 * x5 ** (-1) / eta_p +
                                  2 * k_CLEAN * x5 / (rho3 * v_cruise * (eta_p)))) ** (-1)) / 0.97))
    y5_min = lambda x5: 167.58 * (((((0.5 * rho3 * v_min_cruise ** 3 * cD_0 * x5 ** (-1) / eta_p +
                                      2 * k_CLEAN * x5 / (rho3 * v_min_cruise * (eta_p)))) ** (-1)) / 0.97))
    x5 = np.arange(0, 40/0.02087825 + 0.1, 0.1)
    
    # Turn: constant velocity
    # 10000 ft
    n = 2  # 60° of turning
    v_A = 100 / 1.94384
    
    x6 = ((((0.5 * rho3 * v_A ** 2 * cL_max) / n)) / 0.97) * 0.02087825
    
    # LND phase
    x7_2 = ((0.5 * rho2 * ((1.3 * vs) ** 2)) * cL_max_LND / 0.93) * 0.02087825  # 0.975
    
    plt.figure()
    plt.plot(x1, y1_2(x1, mu[0]), color=[0.8500, 0.3250, 0.0980], linewidth=2.5)
    plt.plot(x1, y1_2(x1, mu[1]), color=[0.8500, 0.3250, 0.0980], linestyle='--', linewidth=3.5)
    plt.axvline(x2, color='m', linewidth=2.5)
    plt.plot(x3 * 0.02087825, y3_1(x3), color=[0.4660, 0.6740, 0.1880], linewidth=2.5)
    plt.plot(x3 * 0.02087825, y3_1_2(x3), color=[0.4660, 0.6740, 0.1880], linestyle='--', linewidth=2.5)
    plt.plot(x3 * 0.02087825, y3_1_c(x3), color=[0.4660, 0.6740, 0.1880], linestyle=':', linewidth=2.5)
    plt.plot(x3 * 0.02087825, y3_2_OEI(x3), color=[0.4660, 0.6740, 0.1880], linestyle='-.', linewidth=2.5)
    plt.plot(x3 * 0.02087825, y3_1_b(x3), color=[0.4660, 0.6740, 0.1880], linestyle='--', linewidth=2.5)
    plt.plot(x4 * 0.02087825, y4(x4), color=[0.9290, 0.6940, 0.1250], linewidth=2.5)
    plt.plot(x5 * 0.02087825, y5(x5), color=[0.4940, 0.1840, 0.5560], linewidth=2.5)
    plt.plot(x5 * 0.02087825, y5_min(x5), color=[0.4940, 0.1840, 0.5560], linestyle='--', linewidth=2.5)
    plt.axvline(x6, color='b', linestyle='-', linewidth=2.5)
    plt.axvline(x7_2, color='c', linestyle='-', linewidth=2.5)
    plt.xlim(12, 29)
    plt.ylim(0, 40)
    plt.scatter(16.8, 12.57, s=80, color='g', marker='o')  # with 5.9 in LND, 4 in TO and vs = 28 kts; n = 2; no SR 22
    A = [12, 12, 13.53, 17.13, 17.13]
    B = [0, 12.33, 13.66, 13.3, 0]
    plt.fill(A, B, color='green', edgecolor='none', alpha=0.3)
    plt.legend(['$TO_{\mu_1}$ 5000 ft', '$TO_{\mu_2}$ 5000 ft', 'Stall', 'Initial climb, AEO, msl',
                'Climb gradient, AEO, msl', 'Climb, AEO', 'Climb, OEI, 5000 ft', 'Balked LND, msl',
                'Service ceiling', 'Target cruise', 'Min cruise', 'Instantaneous turn',
                'LND 5000 ft', 'Design point'],
               loc='best', fontsize=12)
    plt.box(True)
    plt.title('Final SMP', fontsize=20)
    plt.xlabel('$(W/S)_{TO}$ [psf]', fontsize=20)
    plt.ylabel('$(W/P)_{TO}$ [lbs/Hp]', fontsize=20)
    plt.grid(True)    
     
if __name__ == '__main__':
    main()
    plt.show()