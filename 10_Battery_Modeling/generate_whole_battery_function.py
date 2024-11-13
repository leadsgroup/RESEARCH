import numpy as np
import matplotlib.pyplot as plt
import RCAIDE
from scipy.interpolate import RegularGridInterpolator
import os
import pandas as pd

def main():

    # ------------ Change Values Here-------------------------
    current            = 2 #A 
    state_of_charge   = 0
    temperature       = 25 # deg C
    # ------------ Change Values Here-------------------------

    full_path = os.path.join(os.path.dirname(__file__), 'nmc_fit_parameters_data.res')
    fit_parameters_data = RCAIDE.load(full_path)
    interpolator = create_interpolator_from_function(fit_parameters_data)
    a,b,c,d,e,f,g = interpolator(current,temperature)
    Voltage = compute_voltage(state_of_charge,a,b,c,d,e,f,g)
    print(Voltage)
    
    return

def create_interpolator_from_function(params):
    # Define the grid of Amps and Temperatures
    amps_values = [2, 4, 6, 8]  
    temp_values = [0, 10, 20, 30, 40, 50] 
   

    data_a = np.array([params['a'][f'Amps_{amp}'] for amp in amps_values])
    data_b = np.array([params['b'][f'Amps_{amp}'] for amp in amps_values])
    data_c = np.array([params['c'][f'Amps_{amp}'] for amp in amps_values])
    data_d = np.array([params['d'][f'Amps_{amp}'] for amp in amps_values])
    data_e = np.array([params['e'][f'Amps_{amp}'] for amp in amps_values])
    data_f = np.array([params['f'][f'Amps_{amp}'] for amp in amps_values])
    data_g = np.array([params['g'][f'Amps_{amp}'] for amp in amps_values])
    

    # Create interpolators for each parameter
    interp_a = RegularGridInterpolator((amps_values, temp_values), data_a)
    interp_b = RegularGridInterpolator((amps_values, temp_values), data_b)
    interp_c = RegularGridInterpolator((amps_values, temp_values), data_c)
    interp_d = RegularGridInterpolator((amps_values, temp_values), data_d)
    interp_e = RegularGridInterpolator((amps_values, temp_values), data_e)
    interp_f = RegularGridInterpolator((amps_values, temp_values), data_f)
    interp_g = RegularGridInterpolator((amps_values, temp_values), data_g)
    
    # Function of functions LOL 
    def interpolator(amps, temp):
        return float(interp_a([amps, temp])), float(interp_b([amps, temp])), float(interp_c([amps, temp])),float(interp_d([amps, temp])),float(interp_e([amps, temp])),float(interp_f([amps, temp])),float(interp_g([amps, temp]))

    return interpolator

def compute_voltage(x, a, b, c,d,e,f,g):
    return f*((b-np.sinh(g*x-a))/d+c*x)+e

if __name__ == '__main__': 
    main()    
    #plt.show()