# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import numpy             as np
from   scipy.integrate   import quad
import matplotlib.pyplot as plt
from   scipy.special     import kelvin

# References: "Analytical Design and Performance Estimation Methods for 
#              Aircraft Permanent Magnet Synchronous Machines" (2023), 
#              Tallerico, T. F. et al., Glenn Research Center, 
#              Cleveland, Ohio
    
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
    #compute_basic_motor_sizing(B_sign, k_w, I_tot, D, L, omega)
    #example_problem()
    emrax_348()
    
    return

def emrax_348():
    
    print("--------------------------------------------------------------")
    print("EMRAX 348")
    print("--------------------------------------------------------------")
    
    corr_factor    = 300                                                   # [-]            correction factor
    omega_max      = 4000* (2 * np.pi / 60)                                # [RPM -> rad/s] max rotor angular velocity    
    tau_max        = 1000                                                  # [Nm]           max torque                                                                
    omega          = 2500* (2 * np.pi / 60)                                # [RPM -> rad/s] rotor angular velocity 
    D_in           = 0.16                                                  # [m]            stator inner diameter
    D_out          = 0.348                                                 # [m]            stator outer diameter
    k_w            = 0.95                                                  # [-]            winding factor
    I_tot          = 375                                                   # [A]            total current that passes through the stator in both axial directions
    L_motor_stack  = 0.112                                                 # [m]            motor stack length    
    T              = 10                                                    # [s]            total time average power is being calculated over
    I              = 375                                                   # [A]            current
    V              = 830                                                   # [V]            voltage   
    N              = 100                                                   # [-]            number of turns   
    Br             = 1.2                                                   # [Wb/m**2]      remnant flux density
    mu_0           = 1.256637061e-6                                        # [N/A**2]       permeability of free space
    mu_r           = 1005                                                  # [N/A**2]       relative permeability of the magnetic material
    l_m            = 0.002                                                 # [m]            thickness of the magnet or the size of the element if multiple elements span a magnet
    l              = 0.02                                                  # [m]            length of the path    
    l_g            = 0.001                                                 # [m]            airgap size
    Br             = 1.2                                                   # [Wb/m**2]      remnant flux density
    alpha_m        = 180*(np.pi/180)                                       # [deg -> rad]   magnet’s pole span angle in electrical degrees
    mu_r           = 1.05                                                  # [N/A**2]       relative permeability of the magnetic material
    A_tip          = 0.000114                                              # [m²]           tooth tip area 
    l_tip          = 0.005                                                 # [m]            tooth tip gap width
    W_sp           = 0.018                                                 # [m]            stator pole width
    Stack          = 0.095                                                 # [m]            stack length    
    p              = 20                                                    # [-]            pole count
    Rm             = 0.169                                                 # [m]            magnet outer radius
    Rr             = 0.101                                                 # [m]            magnet inner radius
    Rs             = 0.090                                                 # [m]            stator inner radius
    T              = 30                                                    # [°C]           temperature in Celsius
    alpha_mag      = 0.001                                                 # [1/°C]         temperature coefficient of magnetic material    
    k              = 0.02                                                  # [-]            Steinmetz coefficient
    alpha          = 1.5                                                   # [-]            Steinmetz coefficient
    beta           = 2.0                                                   # [-]            Steinmetz coefficient
    B              = 1.2                                                   # [T]            peak magnetic flux density
    f              = 400                                                   # [Hz]           magnetic field frequency
    k_h            = 0.01                                                  # [W/kg]         hysteresis loss coefficient
    k_c            = 0.0001                                                # [W/kg]         eddy current loss coefficient
    k_e            = 0.001                                                 # [W/kg]         excess loss coefficient
    B_sign         = 1.2                                                   # [T]            sign flux density
    t_back         = 0.005                                                 # [m]            back iron thickness
    w_tooth        = 0.01                                                  # [m]            tooth width
    Rs             = 0.17                                                  # [m]            stator inner radius
    p              = 10                                                    # [-]            number of pole pairs
    slots          = 24                                                    # [-]            number of stator slots
    theta          = 0                                                     # [rad]          angular position
    Layers         = 2                                                     # [-]            number of winding layers per slot    
    rho_copper     = 1.68e-8                                               # [Ω·m]          resistivity of copper                
    L_layer        = 0.9                                                   # [m]            length of a winding layer                
    A_layer        = 1.57e-5                                               # [m²]           cross-sectional area of a winding layer                
    SF             = 0.4                                                   # [-]            slot fill factor
    d              = 0.001                                                 # [m]            diameter of the conductor
    delta          = 0.0001                                                # [m]            skin depth at the applied frequency
    sigma          = 5.8e7                                                 # [S/m]          electrical conductivity of the material
    H_e            = 1000                                                  # [A/m]          peak value of the external magnetic field
    n              = 1                                                     # [-]            number of strands
    b              = 0.02                                                  # [m]            winding width
    Rdc            = 0.1                                                   # [Ω]            DC resistance
    m              = 1                                                     # [-]            number of conductor layers in the slot
    EMF_i          = 900                                                   # [V]            per-phase back electromotive force of the motor
    X_d            = 0.18                                                  # [Ω]            d-axis reactance of the motor
    X_q            = 0.1                                                   # [Ω]            q-axis reactance of the motor
    R_s            = 0.03                                                  # [Ω]            stator resistance
    I_q            = 375                                                   # [A]            q-axis current
    I_d            = 375                                                   # [A]            d-axis current
    I_s            = 375                                                   # [A]            motor supply current
    Nt             = 10                                                    # [-]            number of turns per coil
    Np             = 3                                                     # [-]            number of phases
    L_aa           = 0.12                                                  # [H]            self-inductance for phase A
    L_bb           = 0.12                                                  # [H]            self-inductance for phase B
    L_cc           = 0.12                                                  # [H]            self-inductance for phase C
    L_ab           = 0.06                                                  # [H]            mutual inductance between phases A and B
    L_ac           = 0.06                                                  # [H]            mutual inductance between phases A and C
    L_ba           = 0.06                                                  # [H]            mutual inductance between phases B and A (symmetry)
    L_bc           = 0.06                                                  # [H]            mutual inductance between phases B and C
    f_elec         = 400                                                   # [Hz]           electrical frequency
    L_l            = 0.01                                                  # [H]            leakage inductance
    L_m            = 0.04                                                  # [H]            mutual inductance
    Slots          = 24                                                    # [-]            number of slots in the stator
    Layers         = 2                                                     # [-]            number of winding layers
    N_p            = 10                                                    # [-]            number of pole pairs
    A_tip          = 0.000114                                              # [m^2]          tooth tip area
    l_tip          = 0.005                                                 # [m]            tooth tip gap width
    W_sp           = 0.018                                                 # [m]            stator pole width
    Stack          = 0.095                                                 # [m]            stack length
    l_m            = 0.002                                                 # [m]            magnet thickness
    l_g            = 0.001                                                 # [m]            air gap length
    mu_0           = 1.256637061e-6                                        # [N/A^2]        permeability of free space 
    rho_copper     = 1.68e-8                                               # [Ω·m]          resistivity of copper
    Fill           = 0.6                                                   # [-]            assumed fill factor
    A_Layer        = 1e-6                                                  # [m²]           cross-sectional area of one layer
    Slots          = 24                                                    # [-]            number of slots in stator
    Layers         = 2                                                     # [-]            number of layers
    N_t            = 100                                                   # [-]            number of turns per phase
    N_p            = 1                                                     # [-]            number of parallel paths
    L_turn         = 0.3                                                   # [m]            average length of one turn
    L_Layer        = 0.66                                                  # [m]            average length of one layer
    I_peak_layer   = 50                                                    # [A]            peak current in one layer
    Slots          = 24                                                    # [-]            number of slots
    Layers         = 2                                                     # [-]            number of winding layers
    N_p            = 1                                                     # [-]            number of parallel paths
    N_t            = 100                                                   # [-]            number of turns per phase
    I_tot          = 150                                                   # [A]            total current per phase
    D              = 0.2                                                   # [m]            stator diameter
    L              = 0.05                                                  # [m]            stator length
    k_w            = 0.85                                                  # [-]            winding factor
    B_g1           = 0.8                                                   # [T]            air gap flux density
    N_t            = 100                                                   # [-]            number of turns per phase
    r_copper       = 1.68e-8                                               # [Ω·m]          resistivity of copper
    L_Layer        = 0.66                                                  # [m]            length of winding layer
    Fill           = 0.6                                                   # [-]            fill factor
    A_Layer        = 1e-6                                                  # [m²]           cross-sectional area of winding layer
    I_peak_layer   = 50                                                    # [A]            peak current per layer
    I_tot          = 150                                                   # [A]            total current per phase
    B_sign         = 0.7                                                   # [T]            sign flux density
    R_a            = 0.1                                                   # [Ω]            resistance per winding
    I_rms          = 100                                                   # [A]            RMS current
    I_rms_layer    = 70                                                    # [A]            RMS current per layer
    mu_0           = 4 * np.pi * 1e-7                                      # [H/m]          permeability of free space
    A_tip          = 1e-6                                                  # [m²]           tip cross-sectional area
    l_tip          = 0.01                                                  # [m]            tip length
    W_sp           = 0.01                                                  # [m]            specific weight factor
    Stack          = 1                                                     # [m]            stack height
    l_m            = 0.01                                                  # [m]            magnetic path length
    l_g            = 0.002                                                 # [m]            air gap length
    P              = 5000                                                  # [W]            real power
    S              = 5500                                                  # [VA]           apparent power
    EMF_i          = 230                                                   # [V]            induced EMF
    V_ph           = 240                                                   # [V]            per-phase voltage
    I_s            = 20                                                    # [A]            current
    V_bus          = 400                                                   # [V]            DC bus voltage
    Delta_T        = 50                                                    # [°C or K]      temperature difference
    R              = 0.2                                                   # [K/W]          thermal resistance
    l              = 0.01                                                  # [m]            length of the thermal path
    k              = 200                                                   # [W/m·K]        thermal conductivity (e.g., for copper)
    A              = 0.001                                                 # [m²]           cross-sectional area
    h              = 50                                                    # [W/m²·K]       convective heat transfer coefficient
    A              = 0.01                                                  # [m²]           heat transfer area
    L              = 0.1                                                   # [m]            characteristic length
    k_f            = 0.6                                                   # [W/m·K]        thermal conductivity of fluid
    Re             = 5000                                                  # [-]            Reynolds number
    Pr             = 0.7                                                   # [-]            Prandtl number
    h_fin          = 0.02                                                  # [m]            fin height
    w_channel      = 0.005                                                 # [m]            channel width
    f              = 0.005                                                 # [-]            friction factor
    rho            = 1.225                                                 # [kg/m³]        density of fluid
    v              = 10                                                    # [m/s]          flow velocity
    D_h            = 0.01                                                  # [m]            hydraulic diameter
    L_channel      = 1                                                     # [m]            channel length
    V_dot          = 0.002                                                 # [m³/s]         volumetric flow rate
    Ta             = 50                                                    # [-]            Taylor number
    G              = 0.03                                                  # [-]            Grashof number
    R_1            = 0.1                                                   # [m]            inner radius
    R_2            = 0.15                                                  # [m]            outer radius
    t              = 0.01                                                  # [m]            wall thickness
    omega          = 300                                                   # [rad/s]        rotational speed
    R_cg           = 0.125                                                 # [m]            center of gravity radius
    M_mag          = 5                                                     # [kg]           magnet mass
    M_iron         = 10                                                    # [kg]           iron mass
    L              = 0.05                                                  # [m]            stack length
    k_t            = 1.2                                                   # [-]            tensile stress factor
    k              = 1                                                     # [-]            windage loss factor
    C_f            = 0.005                                                 # [-]            friction coefficient (initial guess)
    rho            = 1.2                                                   # [kg/m³]        air density
    r              = 0.15                                                  # [m]            rotor radius
    Stack          = 0.05                                                  # [m]            stack length
    R              = 0.1                                                   # [m]            radius for calculations
    ag             = 0.002                                                 # [m]            air gap
    nu             = 1.5e-5                                                # [m²/s]         kinematic viscosity of air
    r_1            = 0.1                                                   # [m]            inner radius
    r_2            = 0.2                                                   # [m]            outer radius
    I_gz           = 0.01                                                  # [kg·m²]        rotational inertia about z-axis
    omega_motor    = 300                                                   # [rad/s]        motor angular velocity
    omega_aircraft = 2                                                     # [rad/s]        aircraft angular velocity
    n              = 1000                                                  # [RPM]          rotational speed
    C              = 5000                                                  # [N]            dynamic load rating
    P              = 2000                                                  # [N]            equivalent load
    k              = 3                                                     # [-]            load-life exponent (3 for ball bearings)
    mu             = 0.01                                                  # [-]            friction coefficient (dimensionless)
    R_m            = 0.05                                                  # [m]            mean radius of the bearing
    EI             = 1e6                                                   # [N·m²]         flexural rigidity
    m              = 0.5                                                   # [kg]           mass of the shaft
    L              = 0.2                                                   # [m]            shaft length
    A_w            = 1e-4                                                  # [m²]           winding cross-sectional area
    E_w            = 1.1e11                                                # [Pa]           modulus of elasticity of winding material (e.g., copper)
    alpha_w        = 1.7e-5                                                # [1/K]          thermal expansion coefficient of winding material
    T_w            = 400                                                   # [K]            winding temperature
    T_0            = 300                                                   # [K]            reference temperature
    alpha_i        = 1.2e-5                                                # [1/K]          thermal expansion coefficient of iron
    T_i            = 350                                                   # [K]            iron temperature
    A_p            = 2e-4                                                  # [m²]           shear stress area
    I_RMS          = 150                                                   # [A]            RMS current
    R_DS_on        = 0.01                                                  # [Ω]            on-resistance
    f_sw           = 20e3                                                  # [Hz]           switching frequency
    V_bus          = 400                                                   # [V]            DC bus voltage
    I_DS_test      = 100                                                   # [A]            test drain-source current
    V_DS_test      = 300                                                   # [V]            test drain-source voltage
    E_on_test      = 1e-3                                                  # [J]            test turn-on energy
    E_off_test     = 1.5e-3                                                # [J]            test turn-off energy
    Q_rr           = 1e-6                                                  # [C]            reverse recovery charge
    V_bus          = 400                                                   # [V]            DC bus voltage
    m_a            = 0.8                                                   # [-]            modulation index
    L              = 1e-3                                                  # [H]            inductance
    f_sw           = 20e3                                                  # [Hz]           switching frequency
    I_0            = 150                                                   # [A]            peak current
    A_w            = 1e-4                                                  # [m²]           winding cross-sectional area
    B_max          = 1.2                                                   # [T]            maximum flux density
    Fill           = 0.6                                                   # [-]            fill factor (dimensionless)
    A_c            = 5e-4                                                  # [m²]           core cross-sectional area
    phi            = 0.5                                                   # [Wb]           magnetic flux
    mu             = 4 * np.pi * 1e-7                                      # [H/m]          magnetic permeability
    r_copper       = 1.68e-8                                               # [Ω·m]          resistivity of copper
    l_w            = 1                                                     # [m]            winding length
    I_ph_rms       = 150                                                   # [A]            phase RMS current
    m_a            = 0.8                                                   # [-]            modulation index (dimensionless)
    PF             = 0.9                                                   # [-]            power factor (dimensionless)
    I_0            = 200                                                   # [A]            peak current
    dV             = 5                                                     # [V]            voltage ripple
    f_sw           = 20e3                                                  # [Hz]           switching frequency
    delta          = np.radians(5)                                         # [radians]      loss angle
    f              = 50                                                    # [Hz]           fundamental frequency
    R_cap          = 0.01                                                  # [Ω]            equivalent series resistance of capacitor
    l_RMS          = 0.01                                                  # [Ω]            RMS length/resistance (assumed context)
               
    # -----------------------------------------------------------------------------------------
    # 3.1 The Basic Motor Sizing Equation
    # -----------------------------------------------------------------------------------------   

    P_max          = omega_max*tau_max/1000                                # [kW]           max power
    B_sign         = 1.2*corr_factor                                       #o[V*s/m**2]     average magnitude of the radial flux density produced by the rotor
    A_sign         = (k_w*I_tot)/(np.pi*D_in)                              # [-]            stator electrical loading (Eq.2)    
    tau            = (np.pi/2)*(B_sign*A_sign)*(D_in**2)*L_motor_stack     # [Nm]           torque (Eq.1) 
    P              = ((omega)*tau)/1000                                    # [kW]           power (Eq.1)                                                
    B_g1           = 4*B_sign/np.pi                                        # [V*s/m**2]     peak magnetic flux density of the fundamental harmonic produced by the rotor (Eq.4)
    K_s1           = A_sign*np.pi/2                                        # [V*s/m**2]     peak fundamental value of the linear current density (Eq.5)      
    P_lipo         = ((omega)*(np.pi/4)*(B_g1*K_s1)*(D_in**2)*L_motor_stack)/1000 # [kW]    power (Eq.3)
     
    print("Power =", "%0.3f" % P, "[kW]") 
    print("Torque =", "%0.3f" % tau, "[Nm]")
    print("Power with Lipo's formula =", "%0.3f" % P_lipo, "[kW]")
         
    # -----------------------------------------------------------------------------------------
    # 3.1.1 Power and Torque Density
    # -----------------------------------------------------------------------------------------    

    torque_density = (1/2)*(B_sign*A_sign)*D_in                            # [N/m**2]       torque per unit volume (Eq.6)  
    power_density  = (omega*torque_density)/1000                           # [kW/m**3]      power per unit rotor volume (Eq.6)   
    
    print("Torque density =", "%0.3f" % torque_density, "[N/m**2]") 
    print("Power density =", "%0.3f" % power_density, "[kW/m**3]")   
    
    # -----------------------------------------------------------------------------------------
    # 3.1.2 Efficiency and Power Factor
    # -----------------------------------------------------------------------------------------     
    
    LOSS_speed     = 0.05*P                                                # [kW]           sum of the magnetic and mechanical power losses at the rotational speed 𝜔
    tau_out        = tau - (LOSS_speed/(omega))                            # [N*m]          output torque (Eq.7)  
    P_out          = ((omega)*tau_out)/1000                                # [W]            output power (Eq.7)  
      
    def i(t):
        I = 375                                                            # [A]            current
        return I * np.sin(2 * np.pi * t)                                                    
                                                                                            
    def volt(t):                                                                               
        V = 830                                                            # [V]            voltage
        return V * np.sin(2 * np.pi * t)                                                    
                                                                                            
    def fun(t):                                                                               
        return i(t) * volt(t)                                                                  
                                                                                            
    P_time, _      = quad(fun, 0, T)                                       # [W*s]          Real power over the time interval
    P_time         = (P_time / T)/1000                                     # [kW]           Real power
    S              = (I*V)/1000                                            # [kW]           apparent power (Eq.10)  
    PF             = P_time/S                                              # [-]            power factor (Eq.8)  
    Q              = np.sqrt(S**2 - P_time**2)                             # [kW]           reactive power (Eq.11)   
    
    print("Output power =", "%0.3f" % P_out, "[kW]") 
    print("Real power =", "%0.3f" % P_time, "[kW]")
    print("Apparent power =", "%0.3f" % S, "[kW]")
    print("Reactive power =", "%0.3f" % Q, "[kW]")
    print("Power Factor =", "%0.3f" % PF)    
    
    # -----------------------------------------------------------------------------------------
    # 3.2.1 Magnetic Reluctance Networks
    # -----------------------------------------------------------------------------------------
       
    A              = np.pi * ((D_out**2 - D_in**2) / 4)                    # [m**2]         cross-sectional area of the reluctance path perpendicular to length 𝑙    
    MMF_coil       = N*I                                                   # [A*turns]      magnetomotive force applied to the reluctance path for a coil (Eq.14)    
    MMF_magnet     = (Br/(mu_0*mu_r))*l_m                                  # [A*turns]      MMF for a magnet (Eq.15)
    R              = l/(A*mu_0*mu_r)                                       # [A*turn/Wb]    reluctance of a given path or given reluctant element (Eq.16) 
    phi            = MMF_coil/R                                            # [Wb]           magnetic flux through the reluctance path (Eq.12)
    B              = phi/A                                                 # [Wb/m**2]      magnetic flux density (Eq.13)   
    
    print("Magnetic flux density =", "%0.3f" % B, "[Wb/m**2]")

    # -----------------------------------------------------------------------------------------
    # 3.2.1.1 Simple Example of North South Array 
    # -----------------------------------------------------------------------------------------  
    
    A_m            = A                                                     # [m**2]         magnet area      
    R_ag           = l_g/(A_m*mu_0)                                        # [A*turn/Wb]    airgap reluctance (Eq.17)
    R_mag          = l_m/(A_m*mu_0*mu_r)                                   # [A*turn/Wb]    magnet reluctance (Eq.17)
    R              = 2*R_ag + 2*R_mag                                      # [A*turn/Wb]    reluctance of the path (Eq.17)
    MMF            = 2*(Br/(mu_0*mu_r))*l_m                                # [A*turns]      magnetomotive force in the circuit (Eq.18)
    phi            = MMF/R                                                 # [Wb]           flux in the path (Eq.19)
    B_m            = (Br*l_m)/(l_g*mu_r + l_m)                             # [Wb/m**2]      airgap field in the gap produced by the magnets (Eq.20)
    B_sign         = B_m                                                   # [Wb/m**2]      for rotors where magnets span the full
    B_g1           = (4/np.pi)*B_m*np.sin(alpha_m/2)                       # [Wb/m**2]      peak flux density magnitude of the fundamental harmonic, for rotors where the magnets do not span the full magnetic pole  (Eq.21)    
     
    print("Airgap field in the gap produced by the magnets =", "%0.3f" % B_m, "[Wb/m**2]")
    print("Airgap field for rotors where magnets span the full =", "%0.3f" % B_sign, "[Wb/m**2]")
    print("Airgap field peak flux density magnitude of the fundamental harmonic, for rotors where the magnets do not span the full magnetic pole =", "%0.3f" % B_g1, "[Wb/m**2]")
    
    # -----------------------------------------------------------------------------------------
    # 3.2.1.2 Stator Inductance Calculation
    # -----------------------------------------------------------------------------------------   
    
    R_tip          = l_tip/(A_tip*mu_0)                                    # reluctance of the tip (Eq.22)
    R_rotor        = (2/3)*(l_m + l_g)/(W_sp*Stack*mu_0)                   # reluctance of the rotor (Eq.22)
    R              = 1/((2/R_tip) + (1/R_rotor))                           # reluctance of one coil (Eq.22)
    L_m            = N**2/R                                                # self-inductance of a single coil (Eq.23)
    
    print("Self-inductance of a single coil  =", "%0.3f" % L_m, "[H]")

    # -----------------------------------------------------------------------------------------
    # 3.2.2 Closed Form Field Solutions
    # ----------------------------------------------------------------------------------------- 
    
    r              = Rs
    
    def B_r(r, theta, Br, p, Rr, Rm, Rs):
        
        B_r_value = ((Br * (p / (p + 1)) * (1 - (Rr / Rm)**(p + 1))) / (1 - ((Rr**(2*p))/(Rs)))) * (((r / Rs)**(p - 1)) * ((Rm / Rs)**(p + 1)) + ((Rm / r)**(p + 1))) * np.cos(p * theta)
        return B_r_value  
    
    Br_result      = B_r(r, theta, Br, p, Rr, Rm, Rs)
    B_sign         = (2/np.pi)*Br_result                                   # average airgap field (Eq.25)    
    
    print("Average airgap field  =", "%0.3f" % B_sign, "[T]")
    
    # -----------------------------------------------------------------------------------------
    # 3.2.3 Magnet Remnant Flux Density Br
    # -----------------------------------------------------------------------------------------
    
    Br_20C         = Br                                                    # remnant flux density at 20 °C
    Br_T           = Br_20C - alpha_mag*(T - 20)*Br_20C                    # magnet remnant flux density at temperature (Eq.26)    
        
    print("Magnet remnant flux density at temperature =", "%0.3f" % Br_T, "[T]") 

    # -----------------------------------------------------------------------------------------
    # 3.3 Magnet Losses
    # -----------------------------------------------------------------------------------------    
     
    P_v_Steinmetz              = k*(f**alpha)*(B**beta)                    # iron loss per volume (Eq.27)                                                          
    P_v_Bertotti               = k_h*f*(B**beta) + k_c*(f**2)*(B**2) + k_e*(f**1.5)*(B**1.5) # iron loss per volume (Eq.28)
    B_back                     = (B_sign/t_back)*((2*np.pi*Rs)/(2*p))      # peak field in the back iron (Eq.29)
    
    if slots > 2*p:
        B_tooth_slots = (B_sign/w_tooth)*((2*np.pi*Rs)/p)*(p/(slots - p))  # peak field in the tooth iron (Eq.30)
    else:
        B_tooth_slots = (B_sign/w_tooth)*((2*np.pi*Rs)/p)                  # peak field in the tooth iron (Eq.30)
        
    f_tooth                    = f*(slots/(2*p))                           # effective frequency of magnetization of the tooth (Eq.31)
    P_v_tooth                  = (f/f_tooth)*k*(f_tooth**alpha)*(B**beta)  # loss in the stator teeth per unit volume (Eq.32)
        
    print("Iron loss per volume with Steinmetz' formulation =", "%0.3f" % P_v_Steinmetz, "[W/m³]") 
    print("Iron loss per volume with Bertotti's formulation =", "%0.3f" % P_v_Bertotti, "[W/m³]") 
    print("Peak field in the back iron =", "%0.3f" % B_back, "[T]")       
    print("Peak field in the tooth iron =", "%0.3f" % B_tooth_slots, "[T]") 
    print("Loss in the stator teeth per unit volume =", "%0.3f" % P_v_tooth, "[W/m³]") 
    
    # -----------------------------------------------------------------------------------------
    # 3.4 Calculating Current and Resistive Losses
    # -----------------------------------------------------------------------------------------    
    
    I_avg_layer                = I_tot/(slots*Layers)                      # average current per slot per winding layer (Eq.33)
    I_peak_layer               = (np.pi/2)*I_avg_layer                     # peak current per winding layer for a 3-phase motor with sinusoidal supply currents (Eq.34)
    I_rms_layer                = (1/np.sqrt(2))*I_peak_layer               # root mean squared current per layer (Eq.35) 
    LOSS_I2R                   = slots*Layers*rho_copper*(L_layer/(SF*A_layer))*I_rms_layer**2 # resistive losses in the motor (Eq.36) 
    rho_copper_T               = rho_copper*(1 + 0.00393*(T - 20))         # resistivity of copper at a given temperature (Eq.37)
    
    print("Resistive losses in the motor =", "%0.3f" % LOSS_I2R, "[W]")  
    print("Resistivity of copper at a given temperature =", "%0.10f" % rho_copper_T, "[Ω·m]")  

    
    # -----------------------------------------------------------------------------------------
    # 3.4.1 AC Winding Loss
    # -----------------------------------------------------------------------------------------    
    
    Ns                         = N                                         # [-]            number of turns
    gamma                      = d/(delta*np.sqrt(2))                      #                (Eq.39)  
    ber_gamma, bei_gamma, der_ber_gamma, der_bei_gamma          = kelvin(gamma)
    ber_2_gamma, bei_2_gamma, _, _                              = kelvin(2 * gamma)
    
    Rac                        = np.real((Rdc/2)-(gamma*(ber_gamma*der_bei_gamma - bei_gamma*der_ber_gamma)/(der_ber_gamma**2 + der_bei_gamma**2))) # AC resistance due to skin effect for round conductors (Eq.38)
    P_prox                     = np.real((2*np.pi*gamma/sigma)*((ber_2_gamma*der_ber_gamma + bei_2_gamma*der_ber_gamma)/(ber_gamma**2 + bei_gamma**2))*H_e**2) # proximity loss per unit stack length in a conductor (Eq.40)
    R_ac_layer                 = -np.real((Rdc/2)*gamma*(ber_gamma*der_bei_gamma - bei_gamma*der_ber_gamma)/(der_ber_gamma**2 + der_bei_gamma**2) - 2*np.pi*((2*m - 1)**2)*(ber_2_gamma*der_ber_gamma + bei_2_gamma*der_ber_gamma)/(ber_gamma**2 + bei_gamma**2)) # AC resistivity of the mth layer of conductors in the slot (Eq.41)
    R_ac_d_less_than_delta     = np.real(Rdc*(1 + ((np.pi*n*Ns)**2)*d**6/(192*(delta**4)*b**2))) # AC resistance of the winding with d<𝛿 (Eq.42)
    P_prox_d_less_than_delta   = np.real((((np.pi**2)*sigma*d**4)/(32))*(f*B)**2) # proximity loss per unit length generated in a round conductor when d<δ (Eq.43)    
    
    print("AC resistance due to skin effect for round conductors =", "%0.3f" % Rac, "[Ω]")  
    print("Proximity loss per unit stack length in a conductor =", "%0.3f" % P_prox, "[W/m]")  
    print("AC resistivity of the mth layer of conductors in the slot =", "%0.3f" % R_ac_layer, "[Ω/m]")  
    print("AC resistance of the winding with d<𝛿 =", "%0.3f" % R_ac_d_less_than_delta, "[Ω]")  
    print("Proximity loss per unit length generated in a round conductor when d<δ =", "%0.3f" % P_prox_d_less_than_delta, "[W/m]")  
  
    
    # -----------------------------------------------------------------------------------------
    # 3.5 Voltage and Turn Count
    # ----------------------------------------------------------------------------------------- 
       
    V_q                        = EMF_i + X_d*I_d + R_s*I_q                 # q axis voltage of the machine (Eq.45)
    V_d                        = X_q*I_q + R_s*I_d                         # d axis voltage of the machine (Eq.46)
    V_ph                       = np.sqrt(V_q**2 + V_d**2)                  # voltage of a machine (Eq.44)
    V_q_max_torque             = EMF_i + R_s*I_s                           # q axis voltage of the machine (Eq.47)
    V_d_max_torque             = X_q*I_s                                   # d axis voltage of the machine (Eq.48)
    V_ph_max_torque            = np.sqrt((EMF_i + R_s*I_s)**2 + (X_q*I_s)**2) # peak per phase motor voltage (Eq.49)   
    
    print("Voltage of the machine =",        "%0.3f" % V_ph, "[V]")  
    print("Peak per phase motor voltage =", "%0.3f" % V_ph_max_torque, "[V]") 
    
    # -----------------------------------------------------------------------------------------
    # 3.5.1 Back EMF
    # ----------------------------------------------------------------------------------------- 
    
    leng           = 2 * L * Nt                                            # [m] length of the winding (calculated from stack length and turns)
    D                          = D_in
    vel                        = omega*D_in/2
    P1                         = (3/2)*EMF_i*I_s                           # back EMF of the motor (Eq.50)    
    P2                         = (2/np.pi)*omega*D*L*k_w*B_sign*Nt*Np*I_s  # back EMF of the motor (Eq.50)    
    P3                         = ((omega*D*L*k_w*B_g1)/2)*Nt*Np*I_s        # back EMF of the motor (Eq.50)    
    EMF_i                      = (4/np.pi)*omega*D*L*k_w*B_sign*Nt         # per phase back electromotive force of the motor (Eq.51) 
    EMF_i1                     = omega*D*L*k_w*B_g1*Nt                     # per phase back electromotive force of the motor (Eq.51)  
    E                          = B*leng*vel                                # flux cutting form of Faraday’s law (Eq.52) 

    print("Back EMF of the motor  =", "%0.3f" % P3, "[W]")   
    print("Per phase back electromotive force of the motor  =", "%0.3f" % EMF_i1, "[V]")  
    print("Flux cutting form of Faraday’s law  =", "%0.3f" % E, "[V]") 
   
    
    # -----------------------------------------------------------------------------------------
    # 3.5.2 Reactance
    # -----------------------------------------------------------------------------------------    

    L_general = [
        [L_aa, L_ab, L_ac],
        [L_ba, L_bb, L_bc],
        [L_ac, L_bc, L_cc]
    ]                                                                      # Inductance matrix (Eq.54)

    L_simplified = [
        [L_l + L_m, -L_m / 2, -L_m / 2],
        [-L_m / 2, L_l + L_m, -L_m / 2],
        [-L_m / 2, -L_m / 2, L_l + L_m]
    ]                                                                      # Simplified inductance matrix for balanced machines (Eq.55)

    L_d_1 = ((Slots*Layers)/(2*N_p))*(3/2)*(N**2)*((2*A_tip*mu_0/l_tip + (3/2)*(W_sp*Stack*mu_0)/(l_m + l_g))) # Direct axis inductance (Eq.57)
    L_d_2 = ((2*N_p)/(Layers*Slots))*(3/2)*(N_t**2)*(2*A_tip*mu_0/l_tip + (3/2)*(W_sp*Stack*mu_0)/(l_m + l_g)) # Direct axis inductance (Eq.57)
    L_d_3 = L_l + (3/2)*L_m                                                # direct axis inductance (Eq.56) 
    N_t   = (Slots*Layers*N/(2*N_p))                                       # number of series connected turns per phase (Eq.58)
    X_d = L*2*np.pi*f_elec                                                 # [Ω] reactance for an inductive load (Eq.53)

    print("Reactance Section Results:")
    print(f"Inductance matrix: {np.array(L_general)}")
    print(f"Simplified Inductance Matrix: {np.array(L_simplified)}")
    print(f"Direct Axis Inductance 1 (L_d): {L_d_1:.6f} [H]")
    print(f"Direct Axis Inductance 2 (L_d): {L_d_2:.6f} [H]")
    print(f"Direct Axis Inductance 3 (L_d): {L_d_3:.6f} [H]")
    print(f"Reactance for an inductive load (X_d): {X_d:.6f} [Ω]")    
    
    ## -----------------------------------------------------------------------------------------
    ## 3.5.3 Resistance
    ## -----------------------------------------------------------------------------------------    
    
    R_s1                        = rho_copper*(L_Layer*2*N_t/(Fill*A_Layer*Slots*Layers)/(2*N_t*N_p)) # Resistance (Eq.59)    
    A_turn1                     = Fill*A_Layer/N                           # Turn Area (Eq.60)    
    A_turn2                     = Fill*A_Layer*Slots*Layers/(2*N_t*N_p)    # Turn Area (Eq.60) 
    R_s2                        = rho_copper*(L_turn*N_t)/A_turn2          # Resistance (Eq.59)   
    
    print(f"Resistance 1 (R_s1): {R_s1:.6f} [Ω]")
    print(f"Turn Area 1 (A_turn1): {A_turn1:.6e} [m²]")
    print(f"Turn Area 2 (A_turn2): {A_turn2:.6e} [m²]")
    print(f"Resistance 2 (R_s2): {R_s2:.6f} [Ω]")    
    
    # -----------------------------------------------------------------------------------------
    # 3.5.4 Current
    # -----------------------------------------------------------------------------------------    
    
    I_s1                        = I_peak_layer/N                           # Current (Eq.61)  
    I_s2                        = Slots*Layers*I_peak_layer/(2*N_p*N_t)    # Current (Eq.61) 
    I_s3                        = (np.pi/2)*I_tot/(2*N_p*N_t)              # Current (Eq.61) 
    
    print(f"Current Method 1 (I_s1): {I_s1:.6f} [A]")
    print(f"Current Method 2 (I_s2): {I_s2:.6f} [A]")
    print(f"Current Method 3 (I_s3): {I_s3:.6f} [A]")    
    
    # -----------------------------------------------------------------------------------------
    # 3.5.5 Turn Count
    # -----------------------------------------------------------------------------------------    
    
    V_ph                        = np.sqrt((omega*D*L*k_w*B_g1*N_t + r_copper*((L_Layer*2*N_t)/(Fill*A_Layer))*I_peak_layer)**2 + ((3/2)*N_t*(2*A_tip*mu_0/l_tip + (3/2)*W_sp*Stack*mu_0/(l_m + l_g))*I_peak_layer)**2) # Peak per phase motor voltage (Eq.62) 
    P1                          = (3/2)*(EMF_i)*I_s                        # Output power (Eq.63)
    P2                          = 0.5*(4/np.pi)*omega*D*L*k_w*B_sign*N_t*N_p*I_s # Output power (Eq.63)
    P3                          = omega*D*L*k_w*B_sign*I_tot/2             # Output power (Eq.63) 
    R_loss1                     = (N_p*R_a*I_rms**2)/(I**2)                # Resistive losses (Eq.64)
    R_loss2                     = (1/(I**2))*(N_p*r_copper*L_Layer*4*N_p*N_t**2)/(Slots*Layers*Fill*A_Layer)*(Slots*Layers*I_peak_layer/(2*N_p*N_t*np.sqrt(2)))**2 # Resistive losses (Eq.64)
    R_loss3                     = (1/(I**2))*Slots*Layers*r_copper*(L_Layer/(Fill*A_Layer))*I_rms_layer**2 # Resistive losses (Eq.64)
    
    print(f"Per-phase Peak Voltage (V_ph): {V_ph:.6f} [V]")
    print(f"Output Power Method 1 (P1): {P1:.6f} [W]")
    print(f"Output Power Method 2 (P2): {P2:.6f} [W]")
    print(f"Output Power Method 3 (P3): {P3:.6f} [W]")
    print(f"Resistive Losses Method 1 (R_loss1): {R_loss1:.6f} [W]")
    print(f"Resistive Losses Method 2 (R_loss2): {R_loss2:.6f} [W]")
    print(f"Resistive Losses Method 3 (R_loss3): {R_loss3:.6f} [W]")    
    
    # -----------------------------------------------------------------------------------------
    # 3.5.6 Power Factor
    # -----------------------------------------------------------------------------------------      
    
    PF1                         = P/S                                      # Power Factor (Eq.65)
    PF2                         = ((3/2)*EMF_i*I_s)/((3/2)*V_ph*I_s)       # Power Factor (Eq.65)
    PF3                         = EMF_i/V_ph                               # Power Factor (Eq.65)
    
    print(f"Power Factor Method 1 (PF1): {PF1:.6f}")
    print(f"Power Factor Method 2 (PF2): {PF2:.6f}")
    print(f"Power Factor Method 3 (PF3): {PF3:.6f}")    
    
    # -----------------------------------------------------------------------------------------
    # 3.5.7 Modulation Index
    # -----------------------------------------------------------------------------------------      
        
    m_a                        = 2*V_ph/V_bus                              # Modulation index (Eq.66)   
    
    print(f"Modulation Index (m_a): {m_a:.6f}")
    
    # -----------------------------------------------------------------------------------------
    # 4.0 Thermal Considerations
    # ----------------------------------------------------------------------------------------- 

    Q                          = Delta_T/R                                 # heat through a thermal path (Eq.67)
    
    print(f"Heat Transfer (Q): {Q:.2f} [W]")

    # -----------------------------------------------------------------------------------------
    # 4.1 Conductive Path Thermal Resistances
    # -----------------------------------------------------------------------------------------    
    
    R                          = l/(k*A)                                   # thermal resistance (Eq.68)  
    
    print(f"Thermal Resistance (R): {R:.6f} [K/W]")

    # -----------------------------------------------------------------------------------------
    # 4.2 Fluid Flow Thermal Resistances
    # -----------------------------------------------------------------------------------------   

    Re_d                       = Re                                        # 
    R                          = 1/(h*A)                                   # thermal resistance (Eq.69)
    Nu                         = h*L/(k_f)                                 # Nusselt number (Eq.70)
    Nu_laminar                 = 0.453*(Re**0.5)*(Pr**(1/3))               # Laminar Nusselt number (Eq.71)
    Nu_turbulent               = 0.0308*(Re**(4/5))*(Pr**(1/3))            # Turbulent Nusselt number (Eq.71)
    
    if Re_d < 3000:
        Nu_Re_d     = 1.051*np.log(h_fin/w_channel) + 2.89                 # Nusselt number for cooling flow in rectangular ducts and Re_d < 3000 (Eq.72)
    else:
        Nu_Re_d     = ((f/8)*(Re_d - 1000)*Pr)/(1 + 12.7*((f/8)**0.5)*(Pr**(2/3) - 1)) # Nusselt number for cooling flow in rectangular ducts and Re_d >= 3000 (Eq.72)
    
    f_laminar                  = 64/Re                                     # Laminar Moody friction factor (Eq.73)
    f_turbulent                = (0.79*np.log(Re) - 1.64)**(-2)            # Turbulent Moody friction factor (Eq.73)
    Delta_P_flow               = ((f*rho*v**2)/(2*D_h))*L_channel          # Flow pressure drop (Eq.74)
    Loss_cooling               = Delta_P_flow*V_dot                        # Flow pressure drop (Eq.75)
    if Ta < 41: 
        Nu_airgap  = 2                                                     # Nusselt number for the airgap convection and Ta < 41 (Eq.76)
    elif Ta> 41 and Ta < 100:
        Nu_airgap = 0.202*(Ta**(0.63))*(Pr**0.27)                          # Nusselt number for the airgap convection and 41 < Ta < 100 (Eq.76)
    else:
        Nu_airgap = 0.386*(Ta**0.5)*(Pr**0.27)                             # Nusselt number for the airgap convection and 100 < Ta (Eq.76)
    
    if G == 0.01:
        Nu_laminar_G          = 7.46*Re**(0.32)                            # Nusselt number for laminar flow and G = 0.01 (Eq.77)
        Nu_turbulent_G        = 0.044*Re**(0.75)                           # Nusselt number for turbulent flow and G = 0.01 (Eq.78)
        
    elif G> 0.02 and G < 0.06:
        Nu_laminar_G     = 0.5*(1 + 5.47*(10**-4)*np.exp(112*G))*(Re**0.5) # Nusselt number for laminar flow and G = 0.02 - 0.06 (Eq.77)
        Nu_turbulent_G  = 0.5*(12.57*np.exp(-33.18*G))*(Re**(0.6 + 25*G**(12/7))) # Nusselt number for turbulent flow and G = 0.02 - 0.06 (Eq.78)
        
    elif G > 0.06:
        Nu_laminar_G = 0.35*(Re**0.5)                                      # Nusselt number for laminar flow and G > 0.06 (Eq.77)
        Nu_turbulent_G = 0.0151*(Re**0.6)                                  # Nusselt number for turbulent flow and G > 0.06 (Eq.78)  
        
    print(f"Thermal Resistance (R): {R:.6f} [K/W]")
    print(f"Nusselt Number (Nu): {Nu:.6f}")
    print(f"Laminar Nusselt Number (Nu_laminar): {Nu_laminar:.6f}")
    print(f"Turbulent Nusselt Number (Nu_turbulent): {Nu_turbulent:.6f}")
    print(f"Nusselt Number for Rectangular Ducts (Nu_Re_d): {Nu_Re_d:.6f}")
    print(f"Laminar Friction Factor (f_laminar): {f_laminar:.6f}")
    print(f"Turbulent Friction Factor (f_turbulent): {f_turbulent:.6f}")
    print(f"Flow Pressure Drop (Delta_P_flow): {Delta_P_flow:.6f} [Pa]")
    print(f"Cooling Loss (Loss_cooling): {Loss_cooling:.6f} [W]")
    print(f"Nusselt Number for Airgap Convection (Nu_airgap): {Nu_airgap:.6f}")
    print(f"Nusselt Number for Laminar Flow (Nu_laminar_G): {Nu_laminar_G:.6f}")
    print(f"Nusselt Number for Turbulent Flow (Nu_turbulent_G): {Nu_turbulent_G:.6f}")     
    
    # -----------------------------------------------------------------------------------------
    # 5.0 Mechanical Considerations
    # -----------------------------------------------------------------------------------------    

    # -----------------------------------------------------------------------------------------
    # 5.1 Magnet Retention
    # -----------------------------------------------------------------------------------------     

    sigma_peak                 = P*((R_2**2 + R_1**2)/(R_1**2 - R_2**2)) # peak hoop stress (Eq. 79)
    sigma_thin                 = P*R_1/t                        # thin-walled hoop stress approximation (Eq. 80)
    P                          = (omega**2)*R_cg*M_mag/(2*np.pi*R_2*L) # pressure(Eq. 81)
    sigma_tensile              = ((omega**2)*R_cg*(M_mag + M_iron)/(t*Stack*P))*k_t # Tensile stress on the iron dovetails holding the magnets (Eq. 82)
    
    print(f"Pressure (P): {P:.6f} [Pa]")
    print(f"Peak Hoop Stress (σ_peak): {sigma_peak:.6f} [Pa]")
    print(f"Thin-Walled Hoop Stress (σ_thin): {sigma_thin:.6f} [Pa]")
    print(f"Tensile Stress on Iron Dovetails (σ_tensile): {sigma_tensile:.6f} [Pa]")    

    # -----------------------------------------------------------------------------------------
    # 5.2 Windage Losses
    # -----------------------------------------------------------------------------------------     

    Loss_wind_1                = k*C_f*np.pi*rho*(omega**3)*(r**4)*Stack # windage power loss (Eq. 83)
    Re_ag                      = omega*R*ag/nu                  # airgap Reynolds number (Eq. 85)            
    C_f_laminar_ag             = 0.515*((ag/R)**0.3)/(Re_ag**0.5) # Laminar skin friction coefficient (Eq. 84)
    C_f_turbulent_ag           = 0.0325*((ag/R)**0.3)/(Re_ag**0.2) # Turbulent skin friction coefficient (Eq. 84)
    Loss_wind_2                = 0.5*C_f*rho*(omega**3)*(r_2**5 - r_1**5) # windage power loss (Eq. 86)
    Re_R                       = omega*(R**2)/nu                # tip speed Reynolds number of the machine (Eq. 88)
    C_f_laminar_tip            = 3.87/(Re_R**0.5)               # Approximated laminar skin friction coefficient (Eq. 87)
    C_f_turbulent_tip          = 0.146/(Re_R**0.2)              # Approximated turbulent skin friction coefficient (Eq. 87)
    
    print(f"Airgap Reynolds Number (Re_ag): {Re_ag:.2f}")
    print(f"Laminar Friction Coefficient (C_f_laminar_ag): {C_f_laminar_ag:.6f}")
    print(f"Turbulent Friction Coefficient (C_f_turbulent_ag): {C_f_turbulent_ag:.6f}")
    print(f"Tip Speed Reynolds Number (Re_R): {Re_R:.2f}")
    print(f"Laminar Friction Coefficient (C_f_laminar_tip): {C_f_laminar_tip:.6f}")
    print(f"Turbulent Friction Coefficient (C_f_turbulent_tip): {C_f_turbulent_tip:.6f}")
    print(f"Windage Power Loss Method 1 (Loss_wind_1): {Loss_wind_1:.6f} [W]")
    print(f"Windage Power Loss Method 2 (Loss_wind_2): {Loss_wind_2:.6f} [W]")    

    # -----------------------------------------------------------------------------------------
    # 5.3 Bearing Sizing
    # -----------------------------------------------------------------------------------------    

    M_gyro                     = I_gz*omega_motor*omega_aircraft # Gyroscopic moment (Eq. 89)
    L_10                       = ((10**6)/(60*n))*(C/P)**k      # Bearing life with 90 percent reliability (Eq. 90)
    M_f                        = mu*R_m*P                       # Friction moment (Eq. 92)
    
    print(f"Gyroscopic Moment (M_gyro): {M_gyro:.6f} [N·m]")
    print(f"Bearing Life (L_10): {L_10:.2f} [hours]")
    print(f"Friction Moment (M_f): {M_f:.6f} [N·m]")    

    # -----------------------------------------------------------------------------------------
    # 5.4 Rotor Dynamics
    # -----------------------------------------------------------------------------------------        

    N_c_disk                   = (60/(2*np.pi))*np.sqrt((192*EI)/(m*L**3)) # Shaft critical speed for a disk rotor at the center of the shaft (Eq. 93)
    N_c_shaft                  = 60*1.57*(np.sqrt(EI/m))/(L**2) # Shaft critical speed for a shaft on its own supported at both ends (Eq. 94)
    
    print(f"Critical Speed for Disk Rotor at Center (N_c_disk): {N_c_disk:.2f} [RPM]")
    print(f"Critical Speed for Shaft Supported at Both Ends (N_c_shaft): {N_c_shaft:.2f} [RPM]")    

    # -----------------------------------------------------------------------------------------
    # 5.5 Coil Thermo-Mechanical Stress
    # -----------------------------------------------------------------------------------------        

    F_shear                    = A_w*E_w*(alpha_w*(T_w - T_0) - alpha_i*(T_i - T_0)) # Force generated at the winding to iron interface (Eq. 96)  
    sigma_shear                = F_shear/A_p                    # shear stress (Eq. 97)  
    
    print(f"Force at Winding-Iron Interface (F_shear): {F_shear:.6f} [N]")
    print(f"Shear Stress (σ_shear): {sigma_shear:.6f} [Pa]")    

    # -----------------------------------------------------------------------------------------
    # 6.0 Inverter Sizing
    # -----------------------------------------------------------------------------------------    

    # -----------------------------------------------------------------------------------------
    # 6.1 MOSFET Losses
    # -----------------------------------------------------------------------------------------      

    P_cond                     = 0.5*(l_RMS**2)*R_DS_on         # Conduction loss for a single switch (Eq. 98)  
    P_on                       = f_sw*((((np.sqrt(2))/(np.pi))*I_RMS*V_bus)/(I_DS_test*V_DS_test))*E_on_test # turn-on loss (Eq. 99)  
    P_off                      = f_sw*((((np.sqrt(2))/(np.pi))*I_RMS*V_bus)/(I_DS_test*V_DS_test))*E_off_test # turn-off loss (Eq. 100)
    P_rr                       = 0.25*Q_rr*V_bus*f_sw           # Reverse recovery loss (Eq. 101)
    P_loss_sw                  = P_cond + P_on + P_off + P_rr   # Sum of all losses for a single switching device (Eq. 102)
    
    print(f"Conduction Loss (P_cond): {P_cond:.6f} [W]")
    print(f"Turn-On Loss (P_on): {P_on:.6f} [W]")
    print(f"Turn-Off Loss (P_off): {P_off:.6f} [W]")
    print(f"Reverse Recovery Loss (P_rr): {P_rr:.6f} [W]")
    print(f"Total Switching Device Loss (P_loss_sw): {P_loss_sw:.6f} [W]")    

    # -----------------------------------------------------------------------------------------
    # 6.2 Ripple Current
    # -----------------------------------------------------------------------------------------   

    I_ripple                   = (V_bus/2)*(m_a/(2*np.sqrt(3)*L))*(1/f_sw) # Ripple current at the output of the inverter (Eq. 103)  
    W_a                        = 2*L*I_0*A_w/(B_max * Fill * A_c) # Minimum winding area within a core (Eq. 104)     
    N_1                        = I_0*L/(B_max*A_c)              # Number of turns (Eq. 105)    
    N_2                        = I_0*L/phi                      # Number of turns (Eq. 105)
    l_path                     = 2*(np.pi**0.5)*(W_a**0.5)      # Length of the flux path in the iron (Eq. 107)
    L_inductor                 = A_c*mu*N/l_path                # Inductor flux path length (Eq. 106)
    Loss_I2R                   = N*r_copper*(l_w/A_w)*(I_0**2)/2 # Resistive loss per inductor in the winding (Eq. 108)
    
    print(f"Ripple Current (I_ripple): {I_ripple:.6f} [A]")
    print(f"Minimum Winding Area (W_a): {W_a:.6e} [m²]")
    print(f"Number of Turns (N_1): {N_1:.2f}")
    print(f"Number of Turns (N_2): {N_2:.2f}")
    print(f"Flux Path Length (l_path): {l_path:.6f} [m]")
    print(f"Inductor Length (L_inductor): {L_inductor:.6f} [H]")
    print(f"Resistive Loss (Loss_I2R): {Loss_I2R:.6f} [W]")    

    # -----------------------------------------------------------------------------------------
    # 6.3 DC Link Capacitor
    # ----------------------------------------------------------------------------------------- 

    I_In_avg                   = (3/4)*I_0*m_a*PF               # Average inverter current (Eq. 111) 
    I_In_rms                   = I_ph_rms*np.sqrt(2*(np.sqrt(3)/np.pi)*m_a*(PF**2 + 0.25)) # Root mean squared inverter input current (Eq. 110) 
    diff                       = max(I_In_rms**2 - I_In_avg**2, 0)  
    I_c_rms                    = np.sqrt(diff)                  # Current in the DC link capacitor (Eq. 109)
    if I_c_rms > 0:
        C                          = I_c_rms/(dV*f_sw)          # Needed capacitance to limit voltage ripple on the supply side to a desired value (Eq. 112)  
        R_cap_f                    = np.tan(delta)/(2*np.pi*f*C)# Equivalent series resistance of the capacitor (Eq. 114)  
        Loss_cap                   = R_cap*I_c_rms**2           # Losses in the capacitor due to the current ripple (Eq. 113)   
    else:
        C = np.nan
        R_cap_f = np.nan
        Loss_cap = np.nan           
    
    print(f"RMS Inverter Input Current (I_In_rms): {I_In_rms:.6f} [A]")
    print(f"Average Inverter Input Current (I_In_avg): {I_In_avg:.6f} [A]")
    print(f"DC Link Capacitor Current (I_c_rms): {I_c_rms:.6f} [A]")
    print(f"Required Capacitance (C): {C:.6e} [F]")
    print(f"Equivalent Series Resistance (R_cap_f): {R_cap_f:.6f} [Ohm]")
    print(f"Capacitor Losses (Loss_cap): {Loss_cap:.6f} [W]")    
    
    return
    

def compute_basic_motor_sizing(B_sign, k_w, I_tot, D, L, omega):
    
    # -----------------------------------------------------------------------------------------
    # 3 Electromagnetic Motor Sizing
    # -----------------------------------------------------------------------------------------     
    
    # -----------------------------------------------------------------------------------------
    # 3.1 The Basic Motor Sizing Equation
    # -----------------------------------------------------------------------------------------   
    
    B_sign         = 1                                          # average magnitude of the radial flux density produced by the rotor
    D              = 1                                          # stator inner diameter
    k_w            = 0.95                                       # winding factor
    I_tot          = 1                                          # total current that passes through the stator in both axial directions
    L              = 1                                          # motor stack length
    omega          = 1                                          # rotor angular velocity    
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
    I              = 1                                          # amplitude of applied current
    V              = 1                                          # amplitude of applied voltage    
    P              = (1/T)*quad(f, 0, T)                        # real power (Eq.9)  
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
    N              = 1                                          # number of turns
    I              = 1                                          # current in each turn      
    Br             = 1                                          # remnant flux density
    mu_0           = 1                                          # permeability of free space
    mu_r           = 1                                          # relative permeability of the magnetic material
    l_m            = 1                                          # thickness of the magnet or the size of the element if multiple elements span a magnet
    l              = 1                                          # length of the path
    A              = 1                                          # cross-sectional area of the reluctance path perpendicular to length 𝑙    
    phi            = MMF/R                                      # magnetic flux through the reluctance path (Eq.12)
    B              = phi/A                                      # magnetic flux density (Eq.13)    
    MMF_coil       = N*I                                        # MMF for a coil (Eq.14)    
    MMF_magnet     = (Br/(mu_0*mu_r))*l_m                       # MMF for a magnet (Eq.15)
    R              = l/(A*mu_0*mu_r)                            # reluctance of a given path or given reluctant element (Eq.16)
    
    # -----------------------------------------------------------------------------------------
    # 3.2.1.1 Simple Example of North South Array (see the function "example_problem()")
    # -----------------------------------------------------------------------------------------    
    
    l_g            = 1                                          # airgap size
    A_m            = 1                                          # magnet area
    alpha_m        = 1                                          # magnet’s pole span angle in electrical degrees
    R_ag           = l_g/(A_m*mu_0)                             # airgap reluctance (Eq.17)
    R_mag          = l_m/(A_m*mu_0*mu_r)                        # magnet reluctance (Eq.17)
    R              = 2*R_ag + 2*R_mag                           # reluctance of the path (Eq.17)
    MMF            = 2*(Br/(mu_0*mu_r))*l_m                     # magnetomotive force in the circuit (Eq.18)
    phi            = MMF/R                                      # flux in the path (Eq.19)
    B_m            = (Br*l_m)/(l_g*mu_r + l_m)                  # airgap field in the gap produced by the magnets (Eq.20)
    B_sign         = B_m                                        # for rotors where magnets span the full
    B_g1           = (4/np.pi)*B_m*np.sin(alpha_m/2)            # peak flux density magnitude of the fundamental harmonic, for rotors where the magnets do not span the full magnetic pole  (Eq.21)
    
    # -----------------------------------------------------------------------------------------
    # 3.2.1.2 Stator Inductance Calculation
    # -----------------------------------------------------------------------------------------   

    N              = 1                                          # number of turns    
    A_tip          = 1                                          # tooth tip area 
    l_tip          = 1                                          # tooth tip gap width
    W_sp           = 1                                          # stator pole width
    R_tip          = l_tip/(A_tip*mu_0)                         # reluctance of the tip (Eq.22)
    R_rotor        = (2/3)*(l_m + l_g)/(W_sp*Stack*mu_0)        # reluctance of the rotor (Eq.22)
    R              = 1/((2/R_tip) + (1/R_rotor))                # reluctance of one coil (Eq.22)
    L_m            = N**2/R                                     # self-inductance of a single coil (Eq.23)
    
    # -----------------------------------------------------------------------------------------
    # 3.2.2 Closed Form Field Solutions
    # ----------------------------------------------------------------------------------------- 
    
    p              = 1                                          # pole count
    Rm             = 1                                          # magnet outer radius
    Rr             = 1                                          # magnet inner radius
    Rs             = 1                                          # stator inner radius
    theta          = 0
    Br             = ((Br*(p/(p + 1))*(1 - (Rr/Rm)**(p + 1)))/(1 - ((Rr**(2*p))/(Rs))))*(((r/Rs)**(p - 1))*((Rm/Rs)**(p + 1)) + ((Rm/r)**(p + 1)))*np.cos(p*theta) # solution for the radial airgap field from iron cored Halbach machines (Eq.24)
    r              = Rs
    Br_Rs_0        = ((Br*(p/(p + 1))*(1 - (Rr/Rm)**(p + 1)))/(1 - ((Rr**(2*p))/(Rs))))*(((r/Rs)**(p - 1))*((Rm/Rs)**(p + 1)) + ((Rm/r)**(p + 1)))*np.cos(p*theta)
    B_sign         = (2/np.pi)*Br_Rs_0                          # average airgap field (Eq.25)
    
    # -----------------------------------------------------------------------------------------
    # 3.2.3 Magnet Remnant Flux Density Br
    # -----------------------------------------------------------------------------------------
    
    Br_20C         = 1                                          # remnant flux density at 20 °C
    T              = 1                                          # temperature in Celsius
    alpha_mag      = 1                                          # material property related to the specific magnet grade
    Br_T           = Br_20C - alpha_mag*(T - 20)*Br_20C         # magnet remnant flux density at temperature (Eq.26)
    
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
    sigma                      = 1                              # material conductivity
    H_e                        = 1                              # peak value of the applied external magnetic field
    n                          = 1                              # number of strands 
    Ns                         = 1                              # number of turns 
    b                          = 1                              # winding width
    Rac                        = (Rdc/2)-(gamma*(ber_gamma*der_bei_gamma - bei_gamma*der_ber_gamma)/(der_ber_gamma**2 + der_bei_gamma**2)) # AC resistance due to skin effect for round conductors (Eq.38)
    P_prox                     = (2*np.pi*gamma/sigma)*((ber_2_gamma*der_ber_gamma + bei_2_gamma*der_ber_gamma)/(ber_gamma**2 + bei_gamma**2))*H_e**2 # proximity loss per unit stack length in a conductor (Eq.40)
    R_ac                       = (Rdc/2)*gamma*(ber_gamma*der_bei_gamma - bei_gamma*der_ber_gamma)/(der_ber_gamma**2 + der_bei_gamma**2) - 2*np.pi*((2*m - 1)**2)*(ber_2_gamma*der_ber_gamma + bei_2_gamma*der_ber_gamma)/(ber_gamma**2 + bei_gamma**2) # AC resistivity of the mth layer of conductors in the slot (Eq.41)
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
    I_d_max_torque             = 0                              # d axis current for max torque      
    I_s                        = 1                              # motor supply current    
    V_q                        = EMF_i + X_d*I_d + R_s*I_q      # q axis voltage of the machine (Eq.45)
    V_d                        = X_q*I_q + R_s*I_d              # d axis voltage of the machine (Eq.46)
    V_ph                       = np.sqrt(V_q**2 + V_d**2)       # peak per phase motor voltage (Eq.44)
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
    
    L_aa                       = 1                              # self inductance
    L_bb                       = 1                              # self inductance
    L_cc                       = 1                              # self inductance
    L_ab                       = 1                              # mutual inductance
    L_ac                       = 1                              # mutual inductance  
    L_ba                       = 1                              # mutual inductance
    L_bc                       = 1                              # mutual inductance     
    X                          = L*2*np.pi*f_elec               # reactance for an inductive load (Eq.53)
    L                          = [[L_aa, L_ab, L_ac],[L_ab, L_bb, L_bc],[L_ac, L_ba, L_cc]] # inductance matrix of a machine (Eq.54)
    L                          = [[L_l+L_m, -L_m/2, -L_m/2],[-L_m/2, L_l+L_m, -L_m/2],[-L_m/2, -L_m/2, L_l+L_m]] # inductance matrix of a machine (Eq.55)
    L_d                        = L_l + (3/2)*L_m                # direct axis inductance (Eq.56)
    N_t                        = (Slots*Layers*N/(2*N_p))       # number of series connected turns per phase (Eq.58)
    L_d                        = ((Slots*Layers)/(2*N_p))*(3/2)*(N**2)*((2*A_tip*mu_0/l_tip + (3/2)*(W_sp*Stack*mu_0)/(l_m + l_g))) # direct axis inductance (Eq.57)
    L_d                        = ((2*N_p)/(Layers*Slots))*(3/2)*(N_t**2)*(2*A_tip*mu_0/l_tip + (3/2)*(W_sp*Stack*mu_0)/(l_m + l_g)) # direct axis inductance (Eq.57)
    
    # -----------------------------------------------------------------------------------------
    # 3.5.3 Resistance
    # -----------------------------------------------------------------------------------------    
    
    R_s                        = rho_copper*(L_Layer*2*N_t/(Fill*A_Layer*Slots*Layers)/(2*N_t*N_p)) # Resistance (Eq.59)    
    A_turn                     = Fill*A_Layer/N                 # Turn Area (Eq.60)    
    A_turn                     = Fill*A_Layer*Slots*Layers/(2*N_t*N_p) # Turn Area (Eq.60) 
    R_s                        = rho_copper*(L_turn*N_t)/A_turn # Resistance (Eq.59)    
    
    # -----------------------------------------------------------------------------------------
    # 3.5.4 Current
    # -----------------------------------------------------------------------------------------    
        
    I_s                        = I_peak_layer/N                 # Current (Eq.61)  
    I_s                        = Slots*Layers*I_peak_layer/(2*N_p*N_t) # Current (Eq.61) 
    I_s                        = (np.pi/2)*I_tot/(2*N_p*N_t)    # Current (Eq.61) 
    
    # -----------------------------------------------------------------------------------------
    # 3.5.5 Turn Count
    # -----------------------------------------------------------------------------------------    
        
    V_ph                       = np.sqrt((omega*D*L*k_w*B_g1*N_t + r_copper*((L_Layer*2*N_t)/(Fill*A_Layer))*I_peak_layer)**2 + ((3/2)*N_t*(2*A_tip*mu_0/l_tip + (3/2)*W_sp*Stack*mu_0/(l_m + l_g))*I_peak_layer)**2) # Peak per phase motor voltage (Eq.62) 
    P                          = (3/2)*(EMF_i)*I_s              # Output power (Eq.63)
    P                          = 0.5*(4/np.pi)*omega*D*L*k_w*B_sign*N_t*N_p*I_s # Output power (Eq.63)
    P                          = omega*D*L*k_w*B_sign*I_tot/2   # Output power (Eq.63) 
    R_loss                     = (N_p*R_a*I_rms**2)/(I**2)      # Resistive losses (Eq.64)
    R_loss                     = (1/(I**2))*(N_p*r_copper*L_Layer*4*N_p*N_t**2)/(Slots*Layers*Fill*A_Layer)*(Slots*Layers*I_peak_layer/(2*N_p*N_t*np.sqrt(2)))**2 # Resistive losses (Eq.64)
    R_loss                     = (1/(I**2))*Slots*Layers*r_copper*(L_Layer/(Fill*A_Layer))*I_rms_layer**2 # Resistive losses (Eq.64)
    
    # -----------------------------------------------------------------------------------------
    # 3.5.6 Power Factor
    # -----------------------------------------------------------------------------------------    
    
    PF                         = P/S                            # Power Factor (Eq.65)
    PF                         = ((3/2)*EMF_i*I_s)/((3/2)*V_ph*I_s) # Power Factor (Eq.65)
    PF                         = EMF_i/V_ph                     # Power Factor (Eq.65)
    
    # -----------------------------------------------------------------------------------------
    # 3.5.7 Modulation Index
    # -----------------------------------------------------------------------------------------    
        
    m_a                        = 2*V_ph/V_bus                   # Modulation index (Eq.66)
    
    # -----------------------------------------------------------------------------------------
    # 4.0 Thermal Considerations
    # -----------------------------------------------------------------------------------------    

    Q                          = Delta_T/R                      # heat through a thermal path (Eq.67)
    
    # -----------------------------------------------------------------------------------------
    # 4.1 Conductive Path Thermal Resistances
    # -----------------------------------------------------------------------------------------    
        
    R                          = l/(k*A)                        # thermal resistance (Eq.68)  
    
    # -----------------------------------------------------------------------------------------
    # 4.2 Fluid Flow Thermal Resistances
    # -----------------------------------------------------------------------------------------    
            
    R                          = 1/(h*A)                        # thermal resistance (Eq.69)
    Nu                         = h*L/(k_f)                      # Nusselt number (Eq.70)
    Nu_laminar                 = 0.453*(Re**0.5)*(Pr**(1/3))    # Laminar Nusselt number (Eq.71)
    Nu_turbulent               = 0.0308*(Re**(4/5))*(Pr**(1/3)) # Turbulent Nusselt number (Eq.71)
    Nu_Re_d_less_than_3000     = 1.051*ln(h_fin/w_channel) + 2.89 # Nusselt number for cooling flow in rectangular ducts and Re_d < 3000 (Eq.72)
    Nu_Re_d_more_than_3000     = ((f/8)*(Re_d - 1000)*Pr)/(1 + 12.7*((f/8)**0.5)*(Pr**(2/3) - 1)) # Nusselt number for cooling flow in rectangular ducts and Re_d >= 3000 (Eq.72)
    f_laminar                  = 64/Re                          # Laminar Moody friction factor (Eq.73)
    f_turbulent                = (0.79*ln(Re) - 1.64)**(-2)     # Turbulent Moody friction factor (Eq.73)
    Delta_P_flow               = ((f*rho*v**2)/(2*D_h))*L_channel # Flow pressure drop (Eq.74)
    Loss_cooling               = Delta_P_flow*V_dot             # Flow pressure drop (Eq.75)
    Nu_Ta_less_than_41         = 2                              # Nusselt number for the airgap convection and Ta < 41 (Eq.76)
    Nu_41_Ta_100               = 0.202*(Ta**(0.63))*(Pr**0.27)  # Nusselt number for the airgap convection and 41 < Ta < 100 (Eq.76)
    Nu_Ta_more_than_100        = 0.386*(Ta**0.5)*(Pr**0.27)     # Nusselt number for the airgap convection and 100 < Ta (Eq.76)
    Nu_laminar_G_0_01          = 7.46*Re**(0.32)                # Nusselt number for laminar flow and G = 0.01 (Eq.77)
    Nu_laminar_G_0_02_0_06     = 0.5*(1 + 5.47*(10**-4)*np.exp(112*G))*(Re**0.5) # Nusselt number for laminar flow and G = 0.02 - 0.06 (Eq.77)
    Nu_laminar_G_more_than_0_06 = 0.35*(Re**0.5)                # Nusselt number for laminar flow and G > 0.06 (Eq.77)
    Nu_turbulent_G_0_01        = 0.044*Re**(0.75)               # Nusselt number for turbulent flow and G = 0.01 (Eq.78)
    Nu_turbulent_G_0_02_0_06   = 0.5*(12.57*np.exp(-33.18*G))*(Re**(0.6 + 25*G**(12/7))) # Nusselt number for turbulent flow and G = 0.02 - 0.06 (Eq.78)
    Nu_turbulent_G_more_than_0_06 = 0.0151*(Re**0.6)            # Nusselt number for turbulent flow and G > 0.06 (Eq.78)
    
    # -----------------------------------------------------------------------------------------
    # 5.0 Mechanical Considerations
    # -----------------------------------------------------------------------------------------    
    
    # -----------------------------------------------------------------------------------------
    # 5.1 Magnet Retention
    # -----------------------------------------------------------------------------------------    
          
    sigma                      = P*((R_2**2 + R_1**2)/(R_1**2 - R_2**2)) # peak hoop stress (Eq. 79)
    sigma                      = P*R_1/t                        # thin-walled hoop stress approximation (Eq. 80)
    P                          = (omega**2)*R_cg*M_mag/(2*np.pi*R_2*L) # pressure(Eq. 81)
    sigma                      = ((omega**2)*R_cg*(M_mag + M_iron)/(t*Stack*P))*k_t # Tensile stress on the iron dovetails holding the magnets (Eq. 82)
    
    # -----------------------------------------------------------------------------------------
    # 5.2 Windage Losses
    # -----------------------------------------------------------------------------------------    
        
    Loss_wind                  = k*C_f*np.pi*rho*(omega**3)*(r**4)*Stack # windage power loss (Eq. 83)
    Re_ag                      = omega*R*ag/nu                  # airgap Reynolds number (Eq. 85)            
    C_f_laminar                = 0.515*((ag/R)**0.3)/(Re_ag**0.5) # Laminar skin friction coefficient (Eq. 84)
    C_f_turbulent              = 0.0325*((ag/R)**0.3)/(Re_ag**0.2) # Turbulent skin friction coefficient (Eq. 84)
    Loss_wind                  = 0.5*C_f*rho*(omega**3)*(r_2**5 - r_1**5) # windage power loss (Eq. 86)
    Re_R                       = omega*(R**2)/nu                # tip speed Reynolds number of the machine (Eq. 88)
    C_f_laminar                = 3.87/(Re_R**0.5)               # Approximated laminar skin friction coefficient (Eq. 87)
    C_f_turbulent              = 0.146/(Re_R**0.2)              # Approximated turbulent skin friction coefficient (Eq. 87)
    
    # -----------------------------------------------------------------------------------------
    # 5.3 Bearing Sizing
    # -----------------------------------------------------------------------------------------    
        
    M_gyro                     = I_gz*omega_motor*omega_aircraft # Gyroscopic moment (Eq. 89)
    L_10                       = ((10**6)/(60*n))*(C/P)**k      # Bearing life with 90 percent reliability (Eq. 90)
    M_f                        = mu*R_m*P                       # Friction moment (Eq. 92)
    
    # -----------------------------------------------------------------------------------------
    # 5.4 Rotor Dynamics
    # -----------------------------------------------------------------------------------------    
            
    N_c                        = (60/(2*np.pi))*np.sqrt((192*EI)/(m*L**3)) # Shaft critical speed for a disk rotor at the center of the shaft (Eq. 93)
    N_c                        = 60*1.57*(np.sqrt(EI/m))/(L**2) # Shaft critical speed for a shaft on its own supported at both ends (Eq. 94)
    
    # -----------------------------------------------------------------------------------------
    # 5.5 Coil Thermo-Mechanical Stress
    # -----------------------------------------------------------------------------------------    
      
    F_shear                    = A_w*E_w*(alpha_w*(T_w - T_0) - alpha_i*(T_i - T_0)) # Force generated at the winding to iron interface (Eq. 96)  
    sigma_shear                = F_shear/A_p                    # shear stress (Eq. 97)  
    
    # -----------------------------------------------------------------------------------------
    # 6.0 Inverter Sizing
    # -----------------------------------------------------------------------------------------    
        
    # -----------------------------------------------------------------------------------------
    # 6.1 MOSFET Losses
    # -----------------------------------------------------------------------------------------    
     
    P_cond                     = 0.5*(l_RMS**2)*R_DS_on         # Conduction loss for a single switch (Eq. 98)  
    P_on                       = f_sw*((((np.sqrt(2))/(np.pi))*I_RMS*V_bus)/(I_DS_test*V_DS_test))*E_on_test # turn-on loss (Eq. 99)  
    P_off                      = f_sw*((((np.sqrt(2))/(np.pi))*I_RMS*V_bus)/(I_DS_test*V_DS_test))*E_off_test # turn-off loss (Eq. 100)
    P_rr                       = 0.25*Q_rr*V_bus*f_sw           # Reverse recovery loss (Eq. 101)
    P_loss_sw                  = P_cond + P_on + P_off + P_rr   # Sum of all losses for a single switching device (Eq. 102)
    
    # -----------------------------------------------------------------------------------------
    # 6.2 Ripple Current
    # -----------------------------------------------------------------------------------------    
    
    I_ripple                   = (V_bus/2)*(m_a/(2*np.sqrt(3)*L))*(1/f_sw) # Ripple current at the output of the inverter (Eq. 103)  
    W_a                        = 2*L*I_0*A_w/(B_max * Fill * A_c) # Minimum winding area within a core (Eq. 104)     
    N                          = I_0*L/(B_max*A_c)              # Number of turns (Eq. 105)    
    N                          = I_0*L/phi                      # Number of turns (Eq. 105)
    l_path                     = 2*(np.pi**0.5)*(W_a**0.5)      # Length of the flux path in the iron (Eq. 107)
    L                          = A_c*mu*N/l_path                # Inductor flux path length (Eq. 106)
    Loss_I2R                   = N*r_copper*(l_w/A_w)*(I_0**2)/2 # Resistive loss per inductor in the winding (Eq. 108)
    
    # -----------------------------------------------------------------------------------------
    # 6.3 DC Link Capacitor
    # -----------------------------------------------------------------------------------------    
    diff                       = max(I_In_rms**2 - I_In_avg**2, 0)  
    I_c_rms                    = np.sqrt(diff)     
    I_c_rms                    = np.sqrt(I_In_rms**2 - I_In_avg**2) # Current in the DC link capacitor (Eq. 109)
    I_In_rms                   = I_ph_rms*np.sqrt(2*(np.sqrt(3)/np.pi)*m_a*(PF**2 + 0.25)) # Root mean squared inverter input current (Eq. 110) 
    I_In_avg                   = (3/4)*I_0*m_a*PF               # Average inverter current (Eq. 111) 
    if I_c_rms > 0:
        C                          = I_c_rms/(dV*f_sw)              # Needed capacitance to limit voltage ripple on the supply side to a desired value (Eq. 112) 
        R_cap_f                    = np.tan(delta)/(2*np.pi*f*C)    # Equivalent series resistance of the capacitor (Eq. 114) 
        Loss_cap                   = R_cap*I_c_rms**2               # Losses in the capacitor due to the current ripple (Eq. 113) 
    else:
        C = np.nan
        R_cap_f = np.nan
        Loss_cap = np.nan           
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