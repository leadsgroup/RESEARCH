import numpy as np
import matplotlib.pyplot as plt
import RCAIDE
from scipy.interpolate import RegularGridInterpolator,RBFInterpolator
import os
import pandas as pd

def main():

    # ------------ Change Values Here-------------------------
    current            = 1.99230576 #A 
    state_of_charge   = 0.95
    temperature       = 25.25405 # deg C
    # ------------ Change Values Here-------------------------

    full_path = os.path.join(os.path.dirname(__file__), 'nmc_fit_parameters_data.res')
    fit_parameters_data = RCAIDE.load(full_path)
    interpolator = create_interpolator_from_function(fit_parameters_data)
    parameter = interpolator(np.array([[current,temperature]]))
    parameter = np.reshape(parameter,(7,1))
    a = parameter[0]
    b = parameter[1]
    c = parameter[2]
    d = parameter[3]
    e = parameter[4]
    f = parameter[5]
    g = parameter[6]

    Voltage = compute_voltage(state_of_charge,a,b,c,d,e,f,g)
    print(Voltage)
    
    return

# def create_interpolator_from_function(params):
#     # Define the grid of Amps and Temperatures
#     amps_values = np.array([2, 4, 6, 8]  )
#     temp_values = np.array([0, 10, 20, 30, 40, 50] )
   

#     data_a = np.array([params['a'][f'Amps_{amp}'] for amp in amps_values])
#     data_b = np.array([params['b'][f'Amps_{amp}'] for amp in amps_values])
#     data_c = np.array([params['c'][f'Amps_{amp}'] for amp in amps_values])
#     data_d = np.array([params['d'][f'Amps_{amp}'] for amp in amps_values])
#     data_e = np.array([params['e'][f'Amps_{amp}'] for amp in amps_values])
#     data_f = np.array([params['f'][f'Amps_{amp}'] for amp in amps_values])
#     data_g = np.array([params['g'][f'Amps_{amp}'] for amp in amps_values])
    

#        # Create interpolators with extrapolation enabled
#     interp_a = RBFInterpolator((amps_values, temp_values), data_a)
#     interp_b = RBFInterpolator((amps_values, temp_values), data_b)
#     interp_c = RBFInterpolator((amps_values, temp_values), data_c)
#     interp_d = RBFInterpolator((amps_values, temp_values), data_d)
#     interp_e = RBFInterpolator((amps_values, temp_values), data_e)
#     interp_f = RBFInterpolator((amps_values, temp_values), data_f)
#     interp_g = RBFInterpolator((amps_values, temp_values), data_g)


    
#     # Function of functions LOL 
#     def interpolator(amps, temp):
#         return float(interp_a([amps, temp])), float(interp_b([amps, temp])), float(interp_c([amps, temp])),float(interp_d([amps, temp])),float(interp_e([amps, temp])),float(interp_f([amps, temp])),float(interp_g([amps, temp]))

#     return interpolator

# def create_interpolator_from_function(params):
#     # Define the grid of Amps and Temperatures
#     amps_values = [2, 4, 6, 8]
#     temp_values = [0, 10, 20, 30, 40, 50]
# # Repeat amps for each temperature, and tile temps for each amp
#     amps_array = np.repeat(amps_values, len(temp_values))
#     temps_array = np.tile(temp_values, len(amps_values))
#     coordinates = np.column_stack((amps_array, temps_array))

    
#     # Stack all parameters into a single 3D array
#     data = np.stack([
#         np.array([params['a'][f'Amps_{amp}'] for amp in amps_values]),
#         np.array([params['b'][f'Amps_{amp}'] for amp in amps_values]),
#         np.array([params['c'][f'Amps_{amp}'] for amp in amps_values]),
#         np.array([params['d'][f'Amps_{amp}'] for amp in amps_values]),
#         np.array([params['e'][f'Amps_{amp}'] for amp in amps_values]),
#         np.array([params['f'][f'Amps_{amp}'] for amp in amps_values]),
#         np.array([params['g'][f'Amps_{amp}'] for amp in amps_values])
#     ], axis=-1)  # Shape: (amps_values, temp_values, 7)
#     flattened_data = data.reshape(-1, 7)  # Shape: (24, 7)
#     # Create a single interpolator that interpolates across all parameters
#     interpolator = RegularGridInterpolator(coordinates, flattened_data)

#     return interpolator

def create_interpolator_from_function(params):
    # Define the grid of Amps and Temperatures as separate arrays
    amps_values = [2, 4, 6, 8]
    temp_values = [0, 10, 20, 30, 40, 50]
    
    # Stack all parameters into a single 3D array (amps_values, temp_values, 7)
    data = np.stack([
        np.array([params['a'][f'Amps_{amp}'] for amp in amps_values]),
        np.array([params['b'][f'Amps_{amp}'] for amp in amps_values]),
        np.array([params['c'][f'Amps_{amp}'] for amp in amps_values]),
        np.array([params['d'][f'Amps_{amp}'] for amp in amps_values]),
        np.array([params['e'][f'Amps_{amp}'] for amp in amps_values]),
        np.array([params['f'][f'Amps_{amp}'] for amp in amps_values]),
        np.array([params['g'][f'Amps_{amp}'] for amp in amps_values])
    ], axis=-1)  # Shape: (4, 6, 7)

    # Create the interpolator across amps and temp, interpolating each of the 7 parameters
    interpolator = RegularGridInterpolator((amps_values, temp_values), data, bounds_error=False, fill_value=None)

    return interpolator

def compute_voltage(x, a, b, c,d,e,f,g):
    return f*((b-np.sinh(g*x-a))/d+c*x)+e

if __name__ == '__main__': 
    main()    
    #plt.show()