from SALib.sample import saltelli
from SALib.analyze import sobol
import pandas as pd
import numpy as np
import os

def main():
    filename = 'consolidated_exit_conditions.xlsx'
    data = load_excel_data(filename)
    sobol_sensitivity_analysis(data)
    return

def load_excel_data(file_name, sheet_name="Sheet1"):
    file_path = os.path.join(os.path.dirname(__file__), file_name)
    return pd.read_excel(file_path, sheet_name=sheet_name)


def sobol_sensitivity_analysis(data):
    # Normalize the input variables using min-max normalization
    HAS_power = (data['HAS_power'] - data['HAS_power'].min()) / (data['HAS_power'].max() - data['HAS_power'].min())
    HEX_power = (data['HEX_power'] - data['HEX_power'].min()) / (data['HEX_power'].max() - data['HEX_power'].min())
    RES_dimensions = (data['RES_dimensions'] - data['RES_dimensions'].min()) / (data['RES_dimensions'].max() - data['RES_dimensions'].min())

    # Replace original columns with normalized values
    data['HAS_power'] = HAS_power
    data['HEX_power'] = HEX_power
    data['RES_dimensions'] = RES_dimensions

    # Prepare the Sobol problem dictionary
    problem = {
        'num_vars': 3,
        'names': ['HAS_power', 'HEX_power', 'RES_dimensions'],
        'bounds': [[0, 1], [0, 1], [0, 1]]  # Normalized bounds
    }

    # Generate Sobol samples
    param_values = saltelli.sample(problem, 2048)

    # Interpolate the Cycle Day values based on normalized input parameter values
    response_variable = np.array([
        interpolate_cycle_day(data, sample) for sample in param_values
    ])
    
    # Perform Sobol sensitivity analysis
    Si = sobol.analyze(problem, response_variable, calc_second_order=True,print_to_console=True)
    
    # Display results
    print("\nSobol Indices:")
    print("S1 (First Order):", Si['S1'])
    print("ST (Total Order):", Si['ST'])
    print("S2 (Second Order):", Si['S2'])
    
    Si.plot()
    return 

def interpolate_cycle_day(data, params):
    """
    Interpolate the Cycle Day value based on the sampled input parameters.
    """
    has_power, hex_power, res_dimensions = params

    # Find the nearest data point in the normalized dataframe
    nearest_row = data.loc[
        ((data['HAS_power'] - has_power).abs() +
         (data['HEX_power'] - hex_power).abs() +
         (data['RES_dimensions'] - res_dimensions).abs()).idxmin()
    ]
    
    
    return nearest_row['Cycle Day']

if __name__ == "__main__":
    main()