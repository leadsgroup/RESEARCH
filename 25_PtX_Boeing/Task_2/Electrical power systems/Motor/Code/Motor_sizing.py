# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import numpy             as np
from   scipy.integrate   import quad
import matplotlib.pyplot as plt

# References: "Analytical Design and Performance Estimation Methods for Aircraft Permanent Magnet Synchronous Machines" (2023), Thomas F. Tallerico, Aaron D. Anderson, Matthew G. Granger, and Jonathan M. Gutknecht, Glenn Research Center, Cleveland, Ohio
    
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
    B_sign         = 1                                          # average magnitude of the radial flux density produced by the rotor
    D              = 1                                          # stator inner diameter
    k_w            = 0.95                                       # winding factor
    I_tot          = 1                                          # total current that passes through the stator in both axial directions
    D              = 1                                          # stator inner diameter
    L              = 1                                          # motor stack length
    omega          = 1                                          # rotor angular velocity
    
    #compute_basic_motor_sizing(B_sign, k_w, I_tot, D, L, omega)
    example_problem()
    
    return

def compute_basic_motor_sizing(B_sign, k_w, I_tot, D, L, omega):
    
    # -----------------------------------------------------------------------------------------
    # 3 Electromagnetic Motor Sizing
    # -----------------------------------------------------------------------------------------     
    
    # -----------------------------------------------------------------------------------------
    # 3.1 The Basic Motor Sizing Equation
    # -----------------------------------------------------------------------------------------       
    
    A_sign         = (k_w*I_tot)/(np.pi*D)                      # stator electrical loading (Eq.2)   
    tau            = (np.pi/2)*(B_sign*A_sign)*(D**2)*L         # torque (Eq.1)
    P_version_1    = omega*tau                                  # motor power (Eq.1)
    B_g1           = 4*B_sign/np.pi                             # peak magnetic flux density of the fundamental harmonic produced by the rotor (Eq.4)
    K_s1           = A_sign*np.pi/2                             # peak fundamental value of the linear current density (Eq.5)  
    P_version_2    = omega*(np.pi/4)*(B_g1*K_s1)*(D**2)*L       # motor power (Eq.3)
    
    # -----------------------------------------------------------------------------------------
    # 3.1.1 Power and Torque Density
    # -----------------------------------------------------------------------------------------    

    torque_density = (1/2)*(B_sign*A_sign)*D                    # torque per unit volume (Eq.6)  
    power_density  = omega*torque_density                       # power per unit rotor volume (Eq.6)  
    
    # -----------------------------------------------------------------------------------------
    # 3.1.2 Efficiency and Power Factor
    # -----------------------------------------------------------------------------------------     
    
    tau_out        = tau - (LOSS_speed/omega)                   # output torque (Eq.7)  
    P_out          = omega*tau_out                              # output power (Eq.7)  

    def f(i, v):
        return i*v
    
    T              = 2                                          # total time average power is being calculated over
    P              = (1/T)*quad(f, 0, T)                        # real power (Eq.9)  
    I              = 1                                          # amplitude of applied current
    V              = 1                                          # amplitude of applied voltage
    S              = I*V                                        # apparent power (Eq.10)  
    PF             = P/S                                        # power factor (Eq.8)  
    Q              = np.sqrt(S**2 - P**2)                       # reactive power (Eq.11)  
    
    # -----------------------------------------------------------------------------------------
    # 3.2 Calculation of the Average Magnitude of the Radial Flux Density Produced by the Rotor, 𝑩
    # -----------------------------------------------------------------------------------------
    
    # -----------------------------------------------------------------------------------------
    # 3.2.1 Magnetic Reluctance Networks
    # -----------------------------------------------------------------------------------------
    
    MMF            = 1                                          # magnetomotive force applied to the reluctance path
    R              = 1                                          # magnetic reluctance of the path
    phi            = MMF/R                                      # magnetic flux through the reluctance path (Eq.12)
    B              = phi/A                                      # magnetic flux density (Eq.13)
    N              = 1                                          # number of turns
    I              = 1                                          # current in each turn       
    MMF_coil       = N*I                                        # MMF for a coil (Eq.14)
    Br             = 1                                          # remnant flux density
    mu_0           = 1                                          # permeability of free space
    mu_r           = 1                                          # relative permeability of the magnetic material
    l_m            = 1                                          # thickness of the magnet or the size of the element if multiple elements span a magnet
    MMF_magnet     = (Br/(mu_0*mu_r))*l_m                       # MMF for a magnet (Eq.15)
    l              = 1                                          # length of the path
    A              = 1                                          # cross-sectional area of the reluctance path perpendicular to length 𝑙
    R              = l/(A*mu_0*mu_r)                            # reluctance of a given path or given reluctant element (Eq.16)
    
    # -----------------------------------------------------------------------------------------
    # 3.2.1.1 Simple Example of North South Array (see the function "example_problem()")
    # -----------------------------------------------------------------------------------------    
    
    l_g            = 1                                          # airgap size
    A_m            = 1                                          # magnet area
    R_ag           = l_g/(A_m*mu_0)                             # airgap reluctance (Eq.17)
    R_mag          = l_m/(A_m*mu_0*mu_r)                        # magnet reluctance (Eq.17)
    R              = 2*R_ag + 2*R_mag                           # reluctance of the path (Eq.17)
    MMF            = 2*(Br/(mu_0*mu_r))*l_m                     # magnetomotive force in the circuit (Eq.18)
    phi            = MMF/R                                      # flux in the path (Eq.19)
    B_m            = (Br*l_m)/(l_g*mu_r + l_m)                  # airgap field in the gap produced by the magnets (Eq.20)
    B_sign         = B_m                                        # for rotors where magnets span the full
    alpha_m        = 1                                          # magnet’s pole span angle in electrical degrees
    B_g1           = (4/np.pi)*B_m*np.sin(alpha_m/2)            # peak flux density magnitude of the fundamental harmonic, for rotors where the magnets do not span the full magnetic pole  (Eq.21)
    
    # -----------------------------------------------------------------------------------------
    # 3.2.1.2 Stator Inductance Calculation
    # -----------------------------------------------------------------------------------------   
    
    A_tip          = 1                                          # tooth tip area 
    l_tip          = 1                                          # tooth tip gap width
    W_sp           = 1                                          # stator pole width
    R_tip          = l_tip/(A_tip*mu_0)                         # reluctance of the tip (Eq.22)
    R_rotor        = (2/3)*(l_m + l_g)/(W_sp*Stack*mu_0)        # reluctance of the rotor (Eq.22)
    R              = 1/((2/R_tip) + (1/R_rotor))                # reluctance of one coil (Eq.22)
    N              = 1                                          # number of turns
    L_m            = N**2/R                                     # self-inductance of a single coil (Eq.23)
    
    # -----------------------------------------------------------------------------------------
    # 3.2.2 Closed Form Field Solutions
    # ----------------------------------------------------------------------------------------- 
    
    p              = 1                                          # pole count
    Rm             = 1                                          # magnet outer radius
    Rr             = 1                                          # magnet inner radius
    Rs             = 1                                          # stator inner radius
    
    Br             = ((Br*(p/(p + 1))*(1 - (Rr/Rm)**(p + 1)))/(1 - ((Rr**(2*p))/(Rs))))*(((r/Rs)**(p - 1))*((Rm/Rs)**(p + 1)) + ((Rm/r)**(p + 1)))*np.cos(p*theta) # solution for the radial airgap field from iron cored Halbach machines (Eq.24)
    
    r              = Rs
    theta          = 0
    Br_Rs_0        = ((Br*(p/(p + 1))*(1 - (Rr/Rm)**(p + 1)))/(1 - ((Rr**(2*p))/(Rs))))*(((r/Rs)**(p - 1))*((Rm/Rs)**(p + 1)) + ((Rm/r)**(p + 1)))*np.cos(p*theta)
    B_sign         = (2/np.pi)*Br_Rs_0                          # average airgap field (Eq.25)
    
    # -----------------------------------------------------------------------------------------
    # 3.2.3 Magnet Remnant Flux Density Br
    # -----------------------------------------------------------------------------------------
    
    Br_20C         = 1                                          # remnant flux density at 20 °C
    T              = 1                                          # temperature in Celsius
    Br_T           = Br_20C - alpha_mag*(T - 20)*Br_20C         # magnet remnant flux density at temperature (Eq.26)
    alpha_mag      = 1                                          # material property related to the specific magnet grade
    
    # -----------------------------------------------------------------------------------------
    # 3.3 Magnet Losses
    # -----------------------------------------------------------------------------------------    
    
    k                          = 1                              # Steinmetz coefficient
    alpha                      = 1                              # Steinmetz coefficient
    beta                       = 1                              # Steinmetz coefficient
    B                          = 1                              # peak magnetic flux density
    f                          = 1                              # frequency of the magnetic field
    P_v_Steinmetz              = k*(f**alpha)*(B**beta)         # iron loss per volume (Eq.27)
    P_v_Bertotti               = k_h*f*(B**beta) + k_c*(f**2)*(B**2) + k_e*(f**1.5)*(B**1.5) # iron loss per volume (Eq.28)
    B_back                     = (B_sign/t_back)*((2*np.pi*Rs)/(2*p))                        # peak field in the back iron (Eq.29)
    B_tooth_slots_less_than_2p = (B_sign/w_tooth)*((2*np.pi*Rs)/p)                           # peak field in the tooth iron (Eq.30)
    B_tooth_slots_more_than_2p = (B_sign/w_tooth)*((2*np.pi*Rs)/p)*(p/(slots - p))           # peak field in the tooth iron (Eq.30)
    f_tooth                    = f_nom*(Slots/(2*p))            # effective frequency of magnetization of the tooth (Eq.31)
    P_v_tooth                  = (f_nom/f_tooth)*k*(f_tooth**alpha)*(B**beta)                # loss in the stator teeth per unit volume (Eq.32)
    
    # -----------------------------------------------------------------------------------------
    # 3.4 Calculating Current and Resistive Losses
    # -----------------------------------------------------------------------------------------    
        
    I_avg_layer                = I_tot/(Slots*Layers)           # average current per slot per winding layer (Eq.33)
    I_peak_layer               = (np.pi/2)*I_avg_layer          # peak current per winding layer for a 3-phase motor with sinusoidal supply currents (Eq.34)
    I_rms_layer                = (1/np.sqrt(2))*I_peak_layer    # root mean squared current per layer (Eq.35) 
    LOSS_I2R                   = Slots*Layers*rho_copper*(L_layer/(SF*A_layer))*I_rms_layer**2 # resistive losses in the motor (Eq.36) 
    rho_copper_T               = rho_copper_20C*(1 + 0.00393*(T - 20)) # resistivity of copper at a given temperature (Eq.37)
    
    # -----------------------------------------------------------------------------------------
    # 3.4.1 AC Winding Loss
    # -----------------------------------------------------------------------------------------    
    
    d                          = 1                              # diameter of the conductor
    delta                      = 1                              # skin depth in the material at the frequency current is being applied to the conductor
    gamma                      = d/(delta*np.sqrt(2))           # (Eq.39)
    Rac                        = (Rdc/2)-(gamma*(ber_gamma*der_bei_gamma - bei_gamma*der_ber_gamma)/(der_ber_gamma**2 + der_bei_gamma**2)) # AC resistance due to skin effect for round conductors (Eq.38)
    sigma                      = 1                              # material conductivity
    H_e                        = 1                              # peak value of the applied external magnetic field
    P_prox                     = (2*np.pi*gamma/sigma)*((ber_2_gamma*der_ber_gamma + bei_2_gamma*der_ber_gamma)/(ber_gamma**2 + bei_gamma**2))*H_e**2 # proximity loss per unit stack length in a conductor (Eq.40)
    R_ac                       = (Rdc/2)*gamma*(ber_gamma*der_bei_gamma - bei_gamma*der_ber_gamma)/(der_ber_gamma**2 + der_bei_gamma**2) - 2*np.pi*((2*m - 1)**2)*(ber_2_gamma*der_ber_gamma + bei_2_gamma*der_ber_gamma)/(ber_gamma**2 + bei_gamma**2) # AC resistivity of the mth layer of conductors in the slot (Eq.41)
    n                          = 1                              # number of strands 
    Ns                         = 1                              # number of turns 
    b                          = 1                              # winding width
    R_ac                       = Rdc*(1 + ((np.pi*n*Ns)**2)*d**6/(192*(delta**4)*b**2)) # AC resistance of the winding with d<𝛿 (Eq.42)
    P_prox                     = (((np.pi**2)*sigma*d**4)/(32))*(f*B)**2 # proximity loss per unit length generated in a round conductor when d<δ (Eq.43)
    
    # -----------------------------------------------------------------------------------------
    # 3.5 Voltage and Turn Count
    # ----------------------------------------------------------------------------------------- 
    
    EMF_i                      = 1                              # per phase back electromotive force of the motor
    X_d                        = 1                              # d axis reactance of the motor
    X_q                        = 1                              # q axis reactance of the motor
    R_s                        = 1                              # stator resistance
    I_q                        = 1                              # q axis current
    I_d                        = 1                              # d axis current   
    V_q                        = EMF_i + X_d*I_d + R_s*I_q      # q axis voltage of the machine (Eq.45)
    V_d                        = X_q*I_q + R_s*I_d              # d axis voltage of the machine (Eq.46)
    V_ph                       = np.sqrt(V_q**2 + V_d**2)       # peak per phase motor voltage (Eq.44)
    I_d_max_torque             = 0                              # d axis current for max torque      
    I_s                        = 1                              # motor supply current 
    V_q_max_torque             = EMF_i + R_s*I_s                # q axis voltage of the machine (Eq.47)
    V_d_max_torque             = X_q*I_s                        # d axis voltage of the machine (Eq.48)
    V_ph_max_torque            = np.sqrt((EMF_i + R_s*I_s)**2 + (X_q*I_s)**2) # peak per phase motor voltage (Eq.49)    

    # -----------------------------------------------------------------------------------------
    # 3.5.1 Back EMF
    # ----------------------------------------------------------------------------------------- 
    
    P                          = (3/2)*EMF_i*I_s                # back EMF of the motor (Eq.50)    
    P                          = (2/np.pi)*omega*D*L*k_w*B_sign*Nt*Np*I_s # back EMF of the motor (Eq.50)    
    P                          = ((omega*D*L*k_w*B_g1)/2)*Nt*Np*I_s # back EMF of the motor (Eq.50)    
    EMF_i                      = (4/np.pi)*omega*D*L*k_w*B_sign*Nt # per phase back electromotive force of the motor (Eq.51) 
    EMF_i                      = omega*D*L*k_w*B_g1*Nt          # per phase back electromotive force of the motor (Eq.51)  
    E                          = B*l*v                          # flux cutting form of Faraday’s law (Eq.52) 
    
    # -----------------------------------------------------------------------------------------
    # 3.5.2 Reactance
    # -----------------------------------------------------------------------------------------    
    
    X                          = L*2*np.pi*f_elec               # reactance for an inductive load (Eq.53)
    L_aa                       = 1                              # self inductance
    L_bb                       = 1                              # self inductance
    L_cc                       = 1                              # self inductance
    L_ab                       = 1                              # mutual inductance
    L_ac                       = 1                              # mutual inductance  
    L_ba                       = 1                              # mutual inductance
    L_bc                       = 1                              # mutual inductance 
    L                          = [[L_aa, L_ab, L_ac],[L_ab, L_bb, L_bc],[L_ac, L_ba, L_cc]] # inductance matrix of a machine (Eq.54)
    L                          = [[L_l+L_m, -L_m/2, -L_m/2],[-L_m/2, L_l+L_m, -L_m/2],[-L_m/2, -L_m/2, L_l+L_m]] # inductance matrix of a machine (Eq.55)
    
    return


def example_problem():
    
    # -----------------------------------------------------------------------------------------
    # Solution and plot
    # -----------------------------------------------------------------------------------------      
    
    l_m            = np.linspace(0,0.02,100)                    # thickness of the magnet or the size of the element if multiple elements span a magnet
    l_g            = 0.001                                      # airgap size
    Br             = 1.2                                        # remnant flux density
    mu_r           = 1.05                                       # relative permeability of the magnetic material
    B_m            = (Br*l_m)/(l_g*mu_r + l_m)                  # airgap field in the gap produced by the magnets
    
    fig, ax1       = plt.subplots()
    line1,         = ax1.plot(l_m, B_m, 'b-', label="B_m")
    ax1.set_xlabel("Magnet Thickness [m]")
    ax1.set_ylabel("B_m [T]", color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True)
    ax1.set_ylim(0, 1.4)
    line2          = ax1.axhline(y=Br, color='r', linestyle='--', label="Br")
                   
    ax2            = ax1.twinx()
    line3,         = ax2.plot(l_m, [bm / lm for bm, lm in zip(B_m, l_m)], 'k-', label="B_m/l_m")
    ax2.set_ylabel("B_m / l_m [T/m]", color='k')
    ax2.tick_params(axis='y', labelcolor='k')
    ax2.set_ylim(0, 1200)
    
    lines          = [line1, line2, line3]
    labels         = [line.get_label() for line in lines]
    plt.legend(lines, labels, loc="upper left")
    plt.title("Effect of Magnet Thickness on Gap Magnetic Field")
    
    return

if __name__ == '__main__': 
    main()
    plt.show()