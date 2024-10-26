# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import numpy             as np
from   scipy.integrate   import quad
import matplotlib.pyplot as plt

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
    
    compute_basic_motor_sizing(B_sign, k_w, I_tot, D, L, omega)
    
    return

def compute_basic_motor_sizing(B_sign, k_w, I_tot, D, L, omega):
    
    # -----------------------------------------------------------------------------------------
    # The Basic Motor Sizing Equation
    # -----------------------------------------------------------------------------------------       
    
    A_sign         = (k_w*I_tot)/(np.pi*D)                      # stator electrical loading    
    tau            = (np.pi/2)*(B_sign*A_sign)*(D**2)*L         # torque
    P_version_1    = omega*tau                                  # motor power
    B_g1           = 4*B_sign/np.pi                             # peak magnetic flux density of the fundamental harmonic produced by the rotor
    K_s1           = A_sign*np.pi/2                             # peak fundamental value of the linear current density    
    
    # -----------------------------------------------------------------------------------------
    # Power and Torque Density
    # -----------------------------------------------------------------------------------------    
    
    P_version_2    = omega*(np.pi/4)*(B_g1*K_s1)*(D**2)*L       # motor power
    torque_density = (1/2)*(B_sign*A_sign)*D                    # torque per unit volume
    power_density  = omega*torque_density                       # power per unit rotor volume 
    
    # -----------------------------------------------------------------------------------------
    # Efficiency and Power Factor
    # -----------------------------------------------------------------------------------------     
    
    tau_out        = tau - (LOSS_speed/omega)                   # output torque
    P_out          = omega*tau_out                              # output power

    def f(i, v):
        return i*v
    
    T              = 2            
    P              = (1/T)*quad(f, 0, T)                        # real power
    I              = 1                                          # amplitude of applied current
    V              = 1                                          # amplitude of applied voltage
    S              = I*V                                        # apparent power
    PZ             = P/S                                        # power factor
    Q              = np.sqrt(S**2 - P**2)                       # reactive power
    
    # -----------------------------------------------------------------------------------------
    # Calculation of the Average Magnitude of the Radial Flux Density Produced by the Rotor, 𝑩
    # -----------------------------------------------------------------------------------------
    
    # -----------------------------------------------------------------------------------------
    # Magnetic Reluctance Networks
    # -----------------------------------------------------------------------------------------
    
    MMF            = 1                                          # magnetomotive force applied to the reluctance path
    R              = 1                                          # magnetic reluctance of the path
    phi            = MMF/R                                      # magnetic flux through the reluctance path
    B              = phi/A                                      # magnetic flux density
    N              = 1                                          # number of turns
    I              = 1                                          # current in each turn       
    MMF_coil       = N*I                                        # MMF for a coil
    Br             = 1                                          # remnant flux density
    mu_0           = 1                                          # permeability of free space
    mu_r           = 1                                          # relative permeability of the magnetic material
    l_m            = 1                                          # thickness of the magnet or the size of the element if multiple elements span a magnet
    MMF_magnet     = (Br/(mu_0*mu_r))*l_m                       # MMF for a magnet
    l              = 1                                          # length of the path
    A              = 1                                          # cross-sectional area of the reluctance path perpendicular to length 𝑙
    R              = l/(A*mu_0*mu_r)                            # reluctance of a given path or given reluctant element
    
    # -----------------------------------------------------------------------------------------
    # Simple Example of North South Array
    # -----------------------------------------------------------------------------------------    
    
    l_g            = 1                                          # airgap size
    A_m            = 1                                          # magnet area
    R_ag           = l_g/(A_m*mu_0)                             # airgap reluctance
    R_mag          = l_m/(A_m*mu_0*mu_r)                        # magnet reluctance
    R              = 2*R_ag + 2*R_mag                           # reluctance of the path
    MMF            = 2*(Br/(mu_0*mu_r))*l_m                     # magnetomotive force in the circuit
    phi            = MMF/R                                      # flux in the path
    B_m            = (Br*l_m)/(l_g*mu_r + l_m)                  # airgap field in the gap produced by the magnets
    B_sign         = B_m                                        # for rotors where magnets span the full
    alpha_m        = 1                                          # magnet’s pole span angle in electrical degrees
    B_g1           = (4/np.pi)*B_m*np.sin(alpha_m/2)            # peak flux density magnitude of the fundamental harmonic, for rotors where the magnets do not span the full magnetic pole
    
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
    
    # -----------------------------------------------------------------------------------------
    # Stator Inductance Calculation
    # -----------------------------------------------------------------------------------------   
    
    A_tip          = 1                                          # tooth tip area 
    l_tip          = 1                                          # tooth tip gap width
    W_sp           = 1                                          # stator pole width
    R_tip          = l_tip/(A_tip*mu_0)                         # reluctance of the tip
    R_rotor        = (2/3)*(l_m + l_g)/(W_sp*Stack*mu_0)        # reluctance of the rotor
    R              = 1/((2/R_tip) + (1/R_rotor))                # reluctance of one coil
    N              = 1                                          # number of turns
    L_m            = N**2/R                                     # self-inductance of a single coil
    
    # -----------------------------------------------------------------------------------------
    # Closed Form Field Solutions
    # ----------------------------------------------------------------------------------------- 
    
    p              = 1                                          # pole count
    Rm             = 1                                          # magnet outer radius
    Rr             = 1                                          # magnet inner radius
    Rs             = 1                                          # stator inner radius
    
    Br             = ((Br*(p/(p + 1))*(1 - (Rr/Rm)**(p + 1)))/(1 - ((Rr**(2*p))/(Rs))))*(((r/Rs)**(p - 1))*((Rm/Rs)**(p + 1)) + ((Rm/r)**(p + 1)))*np.cos(p*theta)
    
    r              = Rs
    theta          = 0
    Br_Rs_0        = ((Br*(p/(p + 1))*(1 - (Rr/Rm)**(p + 1)))/(1 - ((Rr**(2*p))/(Rs))))*(((r/Rs)**(p - 1))*((Rm/Rs)**(p + 1)) + ((Rm/r)**(p + 1)))*np.cos(p*theta)
    B_sign         = (2/np.pi)*Br_Rs_0                          # average airgap field
    
    # -----------------------------------------------------------------------------------------
    # Magnet Remnant Flux Density Br
    # -----------------------------------------------------------------------------------------
    
    Br_20C         = 1                                          # remnant flux density at 20 °C
    T              = 1                                          # temperature in Celsius
    Br_T           = Br_20C - alpha_mag*(T - 20)*Br_20C         # magnet remnant flux density at temperature
    alpha_mag      = 1                                          # material property related to the specific magnet grade
    
    # -----------------------------------------------------------------------------------------
    # Magnet Losses
    # -----------------------------------------------------------------------------------------    
    
    k                          = 1                                          # Steinmetz coefficient
    alpha                      = 1                                          # Steinmetz coefficient
    beta                       = 1                                          # Steinmetz coefficient
    B                          = 1                                          # peak magnetic flux density
    f                          = 1                                          # frequency of the magnetic field
    P_v_Steinmetz              = k*(f**alpha)*(B**beta)                     # iron loss per volume
    P_v_Bertotti               = k_h*f*(B**beta) + k_c*(f**2)*(B**2) + k_e*(f**1.5)*(B**1.5) # iron loss per volume
    B_back                     = (B_sign/t_back)*((2*np.pi*Rs)/(2*p))
    B_tooth_slots_less_than_2p = (B_sign/w_tooth)*((2*np.pi*Rs)/p)
    B_tooth_slots_more_than_2p = (B_sign/w_tooth)*((2*np.pi*Rs)/p)*(p/(slots - p))
    f_tooth                    = f_nom*(Slots/(2*p))
    P_v_tooth                  = (f_nom/f_tooth)*k*(f_tooth**alpha)*(B**beta)
    
    # -----------------------------------------------------------------------------------------
    # Calculating Current and Resistive Losses
    # -----------------------------------------------------------------------------------------    
        
    I_avg_layer                = I_tot/(Slots*Layers)
    I_peak_layer               = (np.pi/2)*I_avg_layer
    I_rms_layer                = (1/np.sqrt(2))*I_peak_layer
    LOSS_I2R                   = Slots*Layers*rho_copper*(L_layer/(SF*A_layer))*I_rms_layer**2
    rho_copper_T               = rho_copper_20C*(1 + 0.00393*(T - 20))
    
    return

if __name__ == '__main__': 
    main()
    plt.show()