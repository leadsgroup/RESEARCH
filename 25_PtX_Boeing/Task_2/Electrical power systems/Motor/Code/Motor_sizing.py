# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
    B_sign = 1                                                  # average magnitude of the radial flux density produced by the rotor
    D      = 1                                                  # stator inner diameter
    k_w    = 0.95                                               # winding factor
    I_tot  = 1                                                  # total current that passes through the stator in both axial directions
    D      = 1                                                  # stator inner diameter
    L      = 1                                                  # motor stack length
    omega  = 1                                                  # rotor angular velocity
    
    compute_basic_motor_sizing(B_sign, k_w, I_tot, D, L, omega)
    
    return

def compute_basic_motor_sizing(B_sign, k_w, I_tot, D, L, omega):
    
    A_sign         = (k_w*I_tot)/(np.pi*D)                              # stator electrical loading    
    tau            = (np.pi/2)*(B_sign*A_sign)*(D**2)*L                 # torque
    P_version_1    = omega*tau                                          # motor power
    B_g1           = 4*B_sign/np.pi                                     # peak magnetic flux density of the fundamental harmonic produced by the rotor
    K_s1           = A_sign*np.pi/2                                     # peak fundamental value of the linear current density
    P_version_2    = omega*(np.pi/4)*(B_g1*K_s1)*(D**2)*L               # motor power
    torque_density = (1/2)*(B_sign*A_sign)*D                            # torque per unit volume
    power_density  = omega*torque_density                               # power per unit rotor volume 
    
    #P_out + LOSS_speed = 
    
    return

if __name__ == '__main__': 
    main()