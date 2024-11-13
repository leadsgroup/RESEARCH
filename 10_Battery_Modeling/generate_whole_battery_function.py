import numpy as np
import matplotlib.pyplot as plt
import RCAIDE
from scipy.interpolate import RegularGridInterpolator
import os
import pandas as pd

def main():

    # ------------ Change Values Here-------------------------
    current            = 2 #A 
    discharge_capcity = -5 # mAh
    temperature       = 25 # deg C
    # ------------ Change Values Here-------------------------

    full_path = os.path.join(os.path.dirname(__file__), 'nmc_fit_parameters_data.res')
    fit_parameters_data = RCAIDE.load(full_path)
    interpolator = create_interpolator_from_function(fit_parameters_data)
    a,b,c = interpolator(current,temperature)
    Voltage = compute_voltage(discharge_capcity,a,b,c)
    print(Voltage)
    
    return

def create_interpolator_from_function(params):
    # Define the grid of Amps and Temperatures
    amps_values = [2, 4, 6, 8]  
    temp_values = [0, 10, 20, 30, 40, 50] 
   

    data_a = np.array([params['a'][f'Amps_{amp}'] for amp in amps_values])
    data_b = np.array([params['b'][f'Amps_{amp}'] for amp in amps_values])
    data_c = np.array([params['c'][f'Amps_{amp}'] for amp in amps_values])

    # Create interpolators for each parameter
    interp_a = RegularGridInterpolator((amps_values, temp_values), data_a)
    interp_b = RegularGridInterpolator((amps_values, temp_values), data_b)
    interp_c = RegularGridInterpolator((amps_values, temp_values), data_c)
    
    # Function of functions LOL 
    def interpolator(amps, temp):
        return float(interp_a([amps, temp])), float(interp_b([amps, temp])), float(interp_c([amps, temp]))

    return interpolator

def compute_voltage(x, a, b, c):
    return a * np.sinh(b * x) + c

if __name__ == '__main__': 
    main()    
    #plt.show()