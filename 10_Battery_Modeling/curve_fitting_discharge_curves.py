import numpy as np
import matplotlib.pyplot as plt
import RCAIDE
from scipy.optimize import curve_fit
import os
import pandas as pd

def main():
    full_path = os.path.join(os.path.dirname(__file__), 'nmc_structured_raw_data.res')

    raw_data = RCAIDE.load(full_path)
    plot_voltage_data(raw_data)
    fit_parameters = fit_and_plot_sinh_curves(raw_data)
    RCAIDE.save(fit_parameters, os.path.join(os.path.dirname(__file__), 'nmc_fit_parameters_data.res'))
    
    return


def fit_and_plot_sinh_curves(raw_data):
    """
    Function to loop through the current levels and temperatures in the raw_data dictionary,
    fit sinh curves to each set of data, and plot the results.
    
    Parameters:
        raw_data (dict): Nested dictionary with current, temperature, and curve data.
    
    Returns:
        dict: Dictionary containing mean errors for each current and temperature combination.
    """
    
    # Initialize a dictionary to store mean errors for each current and temperature combination
    mean_errors = {}
    fit_params = {'a': {}, 'b': {}, 'c': {}, 'd': {}, 'e': {}, 'f': {}, 'g': {}}

    # Iterate through each current level
    for amp_key, temp_dict in raw_data['Current'].items():
        # Set up a subplot for each temperature level
        n_temps = len(temp_dict)
        fig, axs = plt.subplots(n_temps, 1, figsize=(10, 6 * n_temps), constrained_layout=True)
        fig.suptitle(f'Sinh Curve Fitting for Current Level: {amp_key}', fontsize=16)
        fit_params['a'][amp_key] = []
        fit_params['b'][amp_key] = []
        fit_params['c'][amp_key] = []
        fit_params['d'][amp_key] = []
        fit_params['e'][amp_key] = []
        fit_params['f'][amp_key] = []
        fit_params['g'][amp_key] = []

        # Iterate through each temperature level for the given current
        for i, (temp_key, data) in enumerate(temp_dict.items()):
            curve_data = data['curve']
            discharge_capacity = np.array([point[0] for point in curve_data])/3300
            voltage = np.array([point[1] for point in curve_data])
            

            initial_guess = [7.2533586, 1.88260754, -4.5058183, 1.37190552,  3.86687984,  0.2072,  9.35559597]
            params, covariance = curve_fit(sinh_model, discharge_capacity, voltage, p0=initial_guess,maxfev=1000000,xtol=1e-16, gtol=1e-16)
            a, b, c,d,e,f,g = params
            fit_params['a'][amp_key].append(a)
            fit_params['b'][amp_key].append(b)
            fit_params['c'][amp_key].append(c)
            fit_params['d'][amp_key].append(d)
            fit_params['e'][amp_key].append(e)
            fit_params['f'][amp_key].append(f)
            fit_params['g'][amp_key].append(g)

            
            # Generate model predictions
            x_fit = np.linspace(min(discharge_capacity), max(discharge_capacity), 200)
            y_fit = sinh_model(x_fit, a, b, c,d,e,f,g)
            
            # Calculate the mean error
            y_pred = sinh_model(discharge_capacity, a, b, c,d,e,f,g)
            mean_error = np.mean(np.abs(voltage - y_pred))
            mean_errors[f"{amp_key}_{temp_key}"] = mean_error
            
            # Plotting
            if n_temps == 1:
                ax = axs  # single subplot case
            else:
                ax = axs[i]
            ax.scatter(discharge_capacity, voltage, label='Data')
            ax.plot(x_fit, y_fit,label=rf'Fit: $y = {f:.4f} \left( \frac{{{b:.4f} - \sinh({g:.4f} \cdot x - {a:.4f})}}{{{d:.4f}}} + {c:.4f} \cdot x \right) + {e:.4f}$', color='red')
            ax.set_xlabel('State of Charge')
            ax.set_ylabel('Voltage (V)')
            ax.set_title(f'Temperature: {temp_key}, Mean Error: {mean_error:.3f}')
            ax.legend()

    return fit_params   
def sinh_model(x, a, b, c,d,e,f,g):
    return f*((b-np.sinh(g*x-a))/d+c*x)+e



def plot_voltage_data(raw_data):
    # Plotting each current and temperature curve
    for amp_key, temp_dict in raw_data['Current'].items():
        plt.figure(figsize=(10, 6))
        for temp_key, data in temp_dict.items():
            curve_data = data['curve']
            discharge_capacity = [point[0] for point in curve_data]
            voltage = [point[1] for point in curve_data]
            plt.plot(discharge_capacity, voltage, label=temp_key)
        
        plt.xlabel('Discharge Capacity')
        plt.ylabel('Voltage')
        plt.title(f'Voltage vs Discharge Capacity for {amp_key}')
        plt.legend(title='Temperature')



if __name__ == '__main__': 
    main()    
    plt.show()