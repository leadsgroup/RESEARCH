# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import numpy as np

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
    B_sign = 1  # average magnitude of the radial flux density produced by the rotor
    A_sign = 1  # stator electrical loading
    D      = 1  # stator inner diameter
    L      = 1  # motor stack length
    omega  = 1  # rotor angular velocity
    
    compute_basic_motor_sizing()
    
    return

def compute_basic_motor_sizing(B_sign, A_sign, D, L, omega):
    
    tau = (np.pi/2)*(B_sign*A_sign)*(D**2)*L    # torque
    P   = omega*tau                             # motor power
    
    return