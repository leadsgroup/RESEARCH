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
from   scipy.interpolate import interp1d
from   sklearn.linear_model import LinearRegression

def main():
    
    # ------------------------------------------------------------------
    #   Load Results
    # ------------------------------------------------------------------    
    
    excel_file_name_1                   = "C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations/025_00_00_40_00_00_Baseline_stick_fixed_cruise_Opt_Results.xlsx"
    Baseline_Opt_Vehicle_1_pkl          = "025_00_00_40_00_00_Baseline_Opt_Vehicle"
    Optimized_Vehicle_1_pkl             = "025_00_00_40_00_00_Optimized_Vehicle"
    Output_Stick_Fixed_1_pkl            = "025_00_00_40_00_00_Output_Stick_Fixed"
    Planform_Optimization_Problem_1_pkl = "025_00_00_40_00_00_Planform_Optimization_Problem"
    excel_data_1                        = read_results(excel_file_name_1                  )
    Optimized_Vehicle_1                 = load_results(Optimized_Vehicle_1_pkl            )
    Output_Stick_Fixed_1                = load_results(Output_Stick_Fixed_1_pkl           )
                               
    excel_file_name_2                   = "C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations/025_00_00_40_00_00_Baseline_stick_fixed_cruise_Opt_Results.xlsx"
    Baseline_Opt_Vehicle_2_pkl          = "025_00_00_40_00_00_Baseline_Opt_Vehicle"
    Optimized_Vehicle_2_pkl             = "025_00_00_40_00_00_Optimized_Vehicle"
    Output_Stick_Fixed_2_pkl            = "025_00_00_40_00_00_Output_Stick_Fixed"
    Planform_Optimization_Problem_2_pkl = "025_00_00_40_00_00_Planform_Optimization_Problem"
    excel_data_2                        = read_results(excel_file_name_2                  )
    Optimized_Vehicle_2                 = load_results(Optimized_Vehicle_2_pkl            )
    Output_Stick_Fixed_2                = load_results(Output_Stick_Fixed_2_pkl           )    
    
    # ------------------------------------------------------------------
    #   Assign Variables 
    # ------------------------------------------------------------------      
    
    mw_root_twist                       = np.array([excel_data_1.mw_root_twist             [0], excel_data_2.mw_root_twist             [0]])
    mw_tip_twist                        = np.array([excel_data_1.mw_tip_twist              [0], excel_data_2.mw_tip_twist              [0]])
    vt_span                             = np.array([excel_data_1.vt_span                   [0], excel_data_2.vt_span                   [0]])
    vt_AR                               = np.array([excel_data_1.vt_AR                     [0], excel_data_2.vt_AR                     [0]])
    ht_span                             = np.array([excel_data_1.ht_span                   [0], excel_data_2.ht_span                   [0]])
    ht_AR                               = np.array([excel_data_1.ht_AR                     [0], excel_data_2.ht_AR                     [0]])
    AoA                                 = np.array([excel_data_1.AoA                       [0], excel_data_2.AoA                       [0]])
    CD                                  = np.array([excel_data_1.CD                        [0], excel_data_2.CD                        [0]])
    CM_residual                         = np.array([excel_data_1.CM_residual               [0], excel_data_2.CM_residual               [0]])
    spiral_criteria                     = np.array([excel_data_1.spiral_criteria           [0], excel_data_2.spiral_criteria           [0]])
    static_margin                       = np.array([excel_data_1.static_margin             [0], excel_data_2.static_margin             [0]])
    CM_alpha                            = np.array([excel_data_1.CM_alpha                  [0], excel_data_2.CM_alpha                  [0]])
    phugoid_damping_ratio               = np.array([excel_data_1.phugoid_damping_ratio     [0], excel_data_2.phugoid_damping_ratio     [0]])
    short_period_damping_ratio          = np.array([excel_data_1.short_period_damping_ratio[0], excel_data_2.short_period_damping_ratio[0]])
    dutch_roll_frequency                = np.array([excel_data_1.dutch_roll_frequency      [0], excel_data_2.dutch_roll_frequency      [0]])
    dutch_roll_damping_ratio            = np.array([excel_data_1.dutch_roll_damping_ratio  [0], excel_data_2.dutch_roll_damping_ratio  [0]])
    spiral_doubling_time                = np.array([excel_data_1.spiral_doubling_time      [0], excel_data_2.spiral_doubling_time      [0]])
    run_time                            = np.array([excel_data_1.run_time                  [0], excel_data_2.run_time                  [0]])
                                        
    x_CG                                = np.array([Optimized_Vehicle_1.mass_properties.center_of_gravity        [0][0], Optimized_Vehicle_2.mass_properties.center_of_gravity        [0][0]])
    I_xx                                = np.array([Optimized_Vehicle_1.mass_properties.moments_of_inertia.tensor[0, 0], Optimized_Vehicle_2.mass_properties.moments_of_inertia.tensor[0, 0]])
    I_yy                                = np.array([Optimized_Vehicle_1.mass_properties.moments_of_inertia.tensor[1, 1], Optimized_Vehicle_2.mass_properties.moments_of_inertia.tensor[1, 1]])
    I_zz                                = np.array([Optimized_Vehicle_1.mass_properties.moments_of_inertia.tensor[2, 2], Optimized_Vehicle_2.mass_properties.moments_of_inertia.tensor[2, 2]])
    

    C_m_alpha                           = np.array([-0.01, -0.02, -0.04,  -0.05,  -0.02,  -0.05,  -0.04, -0.05])
    x_CG                                = np.array([1,       1.1,   1.2,    1.3,    1.4,    1.5,    1.6,   1.7])
    
    # ------------------------------------------------------------------
    #   Interpolation 
    # ------------------------------------------------------------------ 

    #interpolating_function = interp1d(x_CG, C_m_alpha, kind='linear', fill_value="extrapolate")
    #x_dense          = np.linspace(x_CG.min(), x_CG.max(), 100)
    #C_m_alpha_new    = interpolating_function(x_dense)
    #d_C_m_alpha_dx   = np.gradient(C_m_alpha_dense, x_dense)
    
    # ------------------------------------------------------------------
    #   Regression 
    # ------------------------------------------------------------------ 
    
    x_CG_reshaped    = x_CG.reshape(-1, 1)
    linear_regressor = LinearRegression()
    linear_regressor.fit(x_CG_reshaped, C_m_alpha)
    slope            = linear_regressor.coef_[0]
    intercept        = linear_regressor.intercept_
    
    def linear_regression_function(x):
        return slope * x + intercept
    
    x_dense          = np.linspace(x_CG.min(), x_CG.max(), 100)
    C_m_alpha_new    = linear_regression_function(x_dense)   
    d_C_m_alpha_dx   = np.full_like(x_dense, slope)
    
    # ------------------------------------------------------------------
    #   Plots
    # ------------------------------------------------------------------ 
    
    fig, ax1 = plt.subplots(figsize=(8, 6))
    ax1.plot(x_dense, C_m_alpha_new, '-', label=r"$C_{m_{\alpha}}$ (Interpolation)", color='blue')
    ax1.scatter(x_CG, C_m_alpha, color='blue', label="Original Points", zorder=5)
    ax1.set_xlabel(r"$x_{CG}$", fontsize=14)
    ax1.set_ylabel(r"$C_{m_{\alpha}}$", fontsize=14, color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid()
    ax2 = ax1.twinx()
    ax2.plot(x_dense, d_C_m_alpha_dx, '--', label=r"$\frac{dC_{m_{\alpha}}}{dx_{CG}}$ (Slope)", color='red')
    ax2.set_ylabel(r"$\frac{dC_{m_{\alpha}}}{dx_{CG}}$", fontsize=14, color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax1.legend(loc='upper left', fontsize=10)
    ax2.legend(loc='upper right', fontsize=10)
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)    
    plt.show()

    return


def load_results(filename):  
    # Define the directory where the file is located
    current_dir = 'C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations'
    
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