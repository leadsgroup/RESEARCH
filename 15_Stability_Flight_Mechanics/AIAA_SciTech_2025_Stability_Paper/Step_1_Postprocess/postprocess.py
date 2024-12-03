# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from   RCAIDE.Framework.Core import Units   
from   RCAIDE.Library.Plots  import *     

# python imports 
import numpy                 as np  
import matplotlib.pyplot     as plt  
import os   
import pickle
import pandas                as pd
from   sklearn.linear_model  import LinearRegression

def main():
    
    # ------------------------------------------------------------------
    #   Load Results
    # ------------------------------------------------------------------  
    
    excel_files = [
        "C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations/025_00_00_40_00_00_Baseline_stick_fixed_cruise_Opt_Results.xlsx",
        "C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations/025_00_00_40_00_00_Baseline_stick_fixed_cruise_Opt_Results.xlsx"
    ]
    pkl_files = {
        "Optimized_Vehicle": [
            "025_00_00_40_00_00_Optimized_Vehicle",
            "025_00_00_40_00_00_Optimized_Vehicle"
        ],
        "Output_Stick_Fixed": [
            "025_00_00_40_00_00_Output_Stick_Fixed",
            "025_00_00_40_00_00_Output_Stick_Fixed"
        ]
    }
    
    # Initialize lists to store data
    mw_root_twist, mw_tip_twist, vt_span, vt_AR, ht_span, ht_AR, AoA, CD, CM_residual, spiral_criteria, static_margin = [], [], [], [], [], [], [], [], [], [], []
    C_m_alpha, phugoid_damping_ratio, short_period_damping_ratio, dutch_roll_frequency, dutch_roll_damping_ratio, spiral_doubling_time, run_time = [], [], [], [], [], [], []
    x_CG, I_xx, I_yy, I_zz = [], [], [], []
    LonModes, LatModes, phugoidDamping, phugoidTimeDoubleHalf, shortPeriodDamping, shortPeriodTimeDoubleHalf = [], [], [], [], [], []
    
    for i in range(len(excel_files)):
        excel_data         = read_results(excel_files[i])
        Optimized_Vehicle  = load_results(pkl_files["Optimized_Vehicle"][i])
        Output_Stick_Fixed = load_results(pkl_files["Output_Stick_Fixed"][i])
        results            = run_mission(Optimized_Vehicle)
        
        mw_root_twist.append(excel_data.mw_root_twist[0])
        mw_tip_twist.append(excel_data.mw_tip_twist[0])
        vt_span.append(excel_data.vt_span[0])
        vt_AR.append(excel_data.vt_AR[0])
        ht_span.append(excel_data.ht_span[0])
        ht_AR.append(excel_data.ht_AR[0])
        AoA.append(excel_data.AoA[0])
        CD.append(excel_data.CD[0])
        CM_residual.append(excel_data.CM_residual[0])
        spiral_criteria.append(excel_data.spiral_criteria[0])
        static_margin.append(excel_data.static_margin[0])
        C_m_alpha.append(results.segments[0].conditions.static_stability.derivatives.CM_alpha[0][0])
        phugoid_damping_ratio.append(excel_data.phugoid_damping_ratio[0])
        short_period_damping_ratio.append(excel_data.short_period_damping_ratio[0])
        dutch_roll_frequency.append(excel_data.dutch_roll_frequency[0])
        dutch_roll_damping_ratio.append(excel_data.dutch_roll_damping_ratio[0])
        spiral_doubling_time.append(excel_data.spiral_doubling_time[0])
        run_time.append(excel_data.run_time[0])
    
        x_CG.append(Optimized_Vehicle.mass_properties.center_of_gravity[0][0])
        I_xx.append(Optimized_Vehicle.mass_properties.moments_of_inertia.tensor[0, 0])
        I_yy.append(Optimized_Vehicle.mass_properties.moments_of_inertia.tensor[1, 1])
        I_zz.append(Optimized_Vehicle.mass_properties.moments_of_inertia.tensor[2, 2])
    
        LonModes.append(results.segments[0].conditions.dynamic_stability.LongModes.LongModes[0])
        LatModes.append(results.segments[0].conditions.dynamic_stability.LatModes.LatModes[0])
        phugoidDamping.append(results.segments[0].conditions.dynamic_stability.LongModes.phugoidDamping[0])
        phugoidTimeDoubleHalf.append(results.segments[0].conditions.dynamic_stability.LongModes.phugoidTimeDoubleHalf[0])
        shortPeriodDamping.append(results.segments[0].conditions.dynamic_stability.LongModes.shortPeriodDamping[0])
        shortPeriodTimeDoubleHalf.append(results.segments[0].conditions.dynamic_stability.LongModes.shortPeriodTimeDoubleHalf[0])
    
    mw_root_twist              = np.array(mw_root_twist)
    mw_tip_twist               = np.array(mw_tip_twist)
    vt_span                    = np.array(vt_span)
    vt_AR                      = np.array(vt_AR)
    ht_span                    = np.array(ht_span)
    ht_AR                      = np.array(ht_AR)
    AoA                        = np.array(AoA)
    CD                         = np.array(CD)
    CM_residual                = np.array(CM_residual)
    spiral_criteria            = np.array(spiral_criteria)
    static_margin              = np.array(static_margin)
    C_m_alpha                  = np.array(C_m_alpha)
    phugoid_damping_ratio      = np.array(phugoid_damping_ratio)
    short_period_damping_ratio = np.array(short_period_damping_ratio)
    dutch_roll_frequency       = np.array(dutch_roll_frequency)
    dutch_roll_damping_ratio   = np.array(dutch_roll_damping_ratio)
    spiral_doubling_time       = np.array(spiral_doubling_time)
    run_time                   = np.array(run_time)
    x_CG                       = np.array(x_CG)
    I_xx                       = np.array(I_xx)
    I_yy                       = np.array(I_yy)
    I_zz                       = np.array(I_zz)
    LonModes                   = np.array(LonModes)
    LatModes                   = np.array(LatModes)
    phugoidDamping             = np.array(phugoidDamping)
    phugoidTimeDoubleHalf      = np.array(phugoidTimeDoubleHalf)
    shortPeriodDamping         = np.array(shortPeriodDamping)
    shortPeriodTimeDoubleHalf  = np.array(shortPeriodTimeDoubleHalf)
    
    # ------------------------------------------------------------------
    #   Regression CG
    # ------------------------------------------------------------------ 
    
    x_CG_reshaped      = x_CG.reshape(-1, 1)
    linear_regressor1  = LinearRegression()
    linear_regressor1.fit(x_CG_reshaped, C_m_alpha)
    slope1             = linear_regressor1.coef_[0]
    intercept1         = linear_regressor1.intercept_
    
    def linear_regression_function1(x):
        return slope1 * x + intercept1
    
    x_dense1          = np.linspace(x_CG.min(), x_CG.max(), 100)
    C_m_alpha_new1    = linear_regression_function1(x_dense1)   
    d_C_m_alpha_dx1   = np.full_like(x_dense1, slope1)
    
    # ------------------------------------------------------------------
    #   CG Plot
    # ------------------------------------------------------------------ 
    
    plot_parameters = plot_style()
    fig, ax1 = plt.subplots(figsize=(8, 6))
    ax1.plot(x_dense1, C_m_alpha_new1, plot_parameters.line_style[0], 
             label=r"$C_{m_{\alpha}}$ (Interpolation)", color='blue', linewidth=plot_parameters.line_width)
    ax1.scatter(x_CG, C_m_alpha, color='blue', label="Original Points", 
                zorder=5, s=plot_parameters.marker_size**2, marker=plot_parameters.markers[0])
    ax1.set_xlabel(r"$x_{CG}$", fontsize=plot_parameters.axis_font_size)
    ax1.set_ylabel(r"$C_{m_{\alpha}}$", fontsize=plot_parameters.axis_font_size, color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid()
    ax2 = ax1.twinx()
    ax2.plot(x_dense1, d_C_m_alpha_dx1, plot_parameters.line_style[1], 
             label=r"$\frac{dC_{m_{\alpha}}}{dx_{CG}}$ (Slope)", color='red', linewidth=plot_parameters.line_width)
    ax2.set_ylabel(r"$\frac{dC_{m_{\alpha}}}{dx_{CG}}$", fontsize=plot_parameters.axis_font_size, color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)   
        
    # ------------------------------------------------------------------
    #   Regression MOI
    # ------------------------------------------------------------------ 
    
    I_yy_reshaped    = I_yy.reshape(-1, 1)
    linear_regressor2 = LinearRegression()
    linear_regressor2.fit(I_yy_reshaped, C_m_alpha)
    slope2            = linear_regressor2.coef_[0]
    intercept2        = linear_regressor2.intercept_
    
    def linear_regression_function2(x):
        return slope2 * x + intercept2
    
    x_dense2          = np.linspace(I_yy.min(), I_yy.max(), 100)
    C_m_alpha_new2    = linear_regression_function2(x_dense2)   
    d_C_m_alpha_dx2   = np.full_like(x_dense2, slope2)    
    
    # ------------------------------------------------------------------
    #   MOI Plot
    # ------------------------------------------------------------------ 
    
    plot_parameters = plot_style()
    
    fig, ax1 = plt.subplots(figsize=(8, 6))
    ax1.plot(x_dense2, C_m_alpha_new2, plot_parameters.line_style[0], 
             label=r"$C_{m_{\alpha}}$ (Interpolation)", color='blue', linewidth=plot_parameters.line_width)
    ax1.scatter(I_yy, C_m_alpha, color='blue', label="Original Points", 
                zorder=5, s=plot_parameters.marker_size**2, marker=plot_parameters.markers[0])
    ax1.set_xlabel(r"$I_{yy}$", fontsize=plot_parameters.axis_font_size)
    ax1.set_ylabel(r"$C_{m_{\alpha}}$", fontsize=plot_parameters.axis_font_size, color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid()
    ax2 = ax1.twinx()
    ax2.plot(x_dense2, d_C_m_alpha_dx2, plot_parameters.line_style[1], 
             label=r"$\frac{dC_{m_{\alpha}}}{dI_{yy}}$ (Slope)", color='red', linewidth=plot_parameters.line_width)
    ax2.set_ylabel(r"$\frac{dC_{m_{\alpha}}}{dI_{yy}}$", fontsize=plot_parameters.axis_font_size, color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    fig.tight_layout()
    fig.subplots_adjust(right=0.85)
    
    # ------------------------------------------------------------------
    #   Root locus Plot
    # ------------------------------------------------------------------
    
    plot_parameters = plot_style()
    real_parts = LatModes.real
    imaginary_parts = LatModes.imag
    

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ['blue', 'green']
    labels = [f"Static Margin {sm:.4f}" for sm in static_margin]
    for i in range(len(static_margin)):
        ax.plot(real_parts[i], imaginary_parts[i], plot_parameters.line_style[0], 
                label=labels[i], color=colors[i], linewidth=plot_parameters.line_width, 
                marker=plot_parameters.markers[i])
        for j in range(1, len(real_parts[i])):
            ax.annotate('', xy=(real_parts[i][j], imaginary_parts[i][j]),
                        xytext=(real_parts[i][j-1], imaginary_parts[i][j-1]),
                        arrowprops=dict(facecolor=colors[i], arrowstyle='->', lw=0.5))
    
    ax.set_xlabel("n (s$^{-1}$)", fontsize=plot_parameters.axis_font_size)
    ax.set_ylabel("$\omega$ (rad/s)", fontsize=plot_parameters.axis_font_size)
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Horizontal line at Imaginary=0
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')  # Vertical line at Real=0
    ax.grid()
    ax.legend(loc="best", fontsize=plot_parameters.legend_fontsize)
    plt.tight_layout()
    plt.show()

    return


def load_results(filename):  
    # Define the directory where the file is located
    #current_dir = '/Users/aidanmolloy/Documents/LEADS/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations'
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

def run_mission(vehicle):
    configs  = configs_setup(vehicle)

    # vehicle analyses
    analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(analyses)
    missions = missions_setup(mission) 

    results = missions.base_mission.evaluate()
    return(results)

def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle() 
    
    # ------------------------------------------------------------------
    #  Weights
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics         = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics)


    # ------------------------------------------------------------------
    #  Stability Analysis
    stability         = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method() 
    stability.vehicle = vehicle
    analyses.append(stability)    

    # ------------------------------------------------------------------
    #  Energy
    energy          = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Framework.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses    

def configs_setup(vehicle):
    '''
    The configration set up below the scheduling of the nacelle angle and vehicle speed.
    Since one prop_rotor operates at varying flight conditions, one must perscribe  the 
    pitch command of the prop_rotor which us used in the variable pitch model in the analyses
    Note: low pitch at take off & low speeds, high pitch at cruise
    '''
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config                                                       = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag                                                   = 'base'     
    configs.append(base_config) 
 
    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------
    config                                            = RCAIDE.Library.Components.Configs.Config(vehicle)
    config.tag                                        = 'forward_flight'   
    vector_angle                                      = 0.0 * Units.degrees   
    config.networks.electric.busses.lift_rotor_bus.active  = False   
    bus = config.networks.electric.busses.prop_rotor_bus  
    for propulsor in  bus.propulsors:
        propulsor.rotor.orientation_euler_angles =  [0, vector_angle, 0]
        propulsor.rotor.pitch_command   = propulsor.rotor.cruise.design_pitch_command  
    configs.append(config)     
    
    return configs
 

# ----------------------------------------------------------------------
#   Define the Mission
# ---------------------------------------------------------------------- 
def mission_setup(analyses): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission        = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag    = 'baseline_mission' 
    
    # unpack Segments module
    Segments       = RCAIDE.Framework.Mission.Segments

    # base segment           
    base_segment   = Segments.Segment()
    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Cruise 
    #------------------------------------------------------------------------------------------------------------------------------------  
    segment                                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                                      = "Cruise"  
    segment.analyses.extend( analyses.forward_flight )                             
    segment.altitude                                                 = 1000.0 * Units.ft  
    segment.air_speed                                                = 110.  * Units['mph']  
    segment.distance                                                 = 40 *Units.nmi 
    segment.true_course                                              = 0 * Units.degree
                                                                     
    # define flight dynamics to model                                
    segment.flight_dynamics.force_x                                  = True  
    segment.flight_dynamics.force_z                                  = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['prop_rotor_propulsor_1','prop_rotor_propulsor_2','prop_rotor_propulsor_3',
                                                                         'prop_rotor_propulsor_4','prop_rotor_propulsor_5','prop_rotor_propulsor_6'] ] 
    segment.assigned_control_variables.body_angle.active             = True                
         
    mission.append_segment(segment)   
   
    return mission 


def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions 

def plot_style():

    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18,
                  'figure.dpi': 100}
    plt.rcParams.update(parameters)
    
    class Data:
        pass
    
    plot_parameters = Data()
    plot_parameters.line_width = 1
    plot_parameters.line_style = ['-', '-']
    plot_parameters.marker_size = 4
    plot_parameters.legend_fontsize = '12'
    plot_parameters.legend_title_font_size = 14
    plot_parameters.axis_font_size = 16
    plot_parameters.title_font_size = 16
    plot_parameters.markers = ['o', 'x', 'o', 'v', 'P', 'p', '^', 'D', '*']
    plot_parameters.color = 'black'
    
    return plot_parameters

if __name__ == '__main__': 
    main()
    plt.show()