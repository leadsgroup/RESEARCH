# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from   RCAIDE.Framework.Core                                                             import Units   
from   RCAIDE.Library.Plots                                                              import *     

# python imports 
import numpy             as np  
import matplotlib.pyplot as plt  
import os   
import pickle
import pandas as pd

# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------

def main():
    
    excel_file_name  = "C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Postprocess/0_25_0____0____0_35_0____0____0_45_0____0____4__0__0__4_1_0___0___4_2_0___0___Baseline_stick_fixed_cruise_Opt_Results.xlsx"
    pickle_file_name = "0_25_0____0____0_35_0____0____0_45_0____0____4__0__0__4_1_0___0___4_2_0___0___Baseline_Opt_Vehicle"
      
    excel_data       = read_results(excel_file_name)
    pickle_data      = load_results(pickle_file_name)  

    debug = 0
    
    return


def load_results(filename):  
    # Define the directory where the file is located
    current_dir = 'C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Postprocess'
    
    # Combine directory and filename to get full path
    load_file = os.path.join(current_dir, filename + '.pkl')
    
    # Open and load the pickle file
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    
    return results

def read_results(file_name):
    # Load the Excel file
    excel_data = pd.ExcelFile(file_name)
    
    # Check sheet names
    print("Available sheets:", excel_data.sheet_names)
    
    # Load data from the 'Stick_Fixed' sheet
    data = excel_data.parse("Stick_Fixed")
    
    return data

if __name__ == '__main__': 
    main()
    plt.show()