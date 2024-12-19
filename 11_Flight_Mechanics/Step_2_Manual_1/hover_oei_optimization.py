
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from   RCAIDE.Framework.Core import Units,Data  
from   RCAIDE.Library.Methods.Propulsors.Converters.Rotor.compute_rotor_performance import compute_rotor_performance 
from   RCAIDE.Library.Plots  import *  
import numpy                 as np
from   scipy.optimize        import fsolve
from   scipy                 import interpolate
from   scipy.optimize        import minimize 
from   scipy.interpolate     import interp1d
import os   
import pickle  
import pylab as plt 
from   copy import deepcopy
import sys 
import pandas                as pd

sys.path.append(os.path.join(os.path.split(os.path.split(os.path.split(sys.path[0])[0])[0])[0], 'Aircraft'))

from   Stopped_Rotor.Stopped_Rotor                                        import vehicle_setup as SR_vehicle_setup   
from   Tiltrotor.Tiltrotor                                                import vehicle_setup as TR_vehicle_setup   
from   Tiltwing.Tiltwing                                                  import vehicle_setup as TW_vehicle_setup   
from   Hexacopter.Hexacopter                                              import vehicle_setup as HC_vehicle_setup   
from   Tilt_Stopped_Rotor.Tilt_Stopped_Rotor_Conv_Tail                    import vehicle_setup as TSR_vehicle_setup 

def main():
    
    # -------------------------------------------------------------
    #  Select aircraft model
    # -------------------------------------------------------------     
    
    aircraft_model  = 'TSR'
    
    if aircraft_model == 'SR':
        vehicle =  SR_vehicle_setup(redesign_rotors=False)
    if aircraft_model == 'TR':
        vehicle =  TR_vehicle_setup(redesign_rotors=False)
    if aircraft_model == 'TW':
        vehicle =  TW_vehicle_setup(redesign_rotors=False)
    if aircraft_model == 'HC':
        vehicle =  HC_vehicle_setup(redesign_rotors=False)
    if aircraft_model == 'TSR':
        vehicle =  TSR_vehicle_setup(redesign_rotors=False)
        
    # -------------------------------------------------------------
    #  Deepcopy vehicle and find excel file location
    # -------------------------------------------------------------         
        
    case_vehicle  = deepcopy(vehicle)
    
    excel_file = "C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations/025_00_00_40_00_00_Baseline_stick_fixed_cruise_Opt_Results.xlsx"
    
    # >>>>>>>>>>>>>>>>>>>>> ADD OTHER FILES
    
    
    # delete control surfaces if they have been defined 
    #for wing in case_vehicle.wings:
        #for control_surface in wing.control_surfaces:
            #del wing.control_surfaces[control_surface.tag]
            
    # -------------------------------------------------------------
    #  Sweeps of battery locations (UPDATE WITH LAT VALUES LATER)
    # -------------------------------------------------------------            
    
    # >>>>>>>>>>>>>>>>>>>>> ADD LATERAL CASES
    
    # prop rotor battery module (first module)
    #                 CG: X,    Y,    Z,  wing_flag 
    CG_bat_1 = np.array([[0.25, 0.,   0., True], # Change wing flag back to false
                         [0.35, 0.,   0., False],
                         [0.45, 0.,   0., False],
                         [1.9375, 0.25, 1.0, True],
                         [1.9375, 0.35, 1.0, True],
                         [1.9375, 0.45, 1.0, True]])
    
    # lift rotor battery modules
    #                 CG: X,    Y,    Z 
    CG_bat_2 = np.array([[4.0,  0.,   0.],
                         [4.1,  0.,   0.],
                         [4.2,  0.,   0.]])    
    
    RPM_1        = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    RPM_2        = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    RPM_3        = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    RPM_4        = np.zeros([len(CG_bat_1),len(CG_bat_2)])  
    
    # -------------------------------------------------------------
    #  Hover OEI
    # -------------------------------------------------------------
    
    #bus                            = RCAIDE.Library.Components.Energy.Distributors.Electrical_Bus() 
    #electric_rotor                 = RCAIDE.Library.Components.Propulsors.Electric_Rotor() 
    #rotor                          = F8745_D4_Propeller() 
    #electric_rotor.rotor           = rotor  
    #bus.propulsors.append(electric_rotor)     
    #segment.state.conditions.energy[bus.tag] = Conditions()
    #segment.state.conditions.noise[bus.tag]  = Conditions()
    #electric_rotor.append_operating_conditions(segment,bus) 
    #for tag, item in  electric_rotor.items(): 
        #if issubclass(type(item), RCAIDE.Library.Components.Component):
            #item.append_operating_conditions(segment,bus,electric_rotor)  
    ## Run BEMT
    #segment.state.conditions.expand_rows(ctrl_pts)
    #rotor_conditions             =  segment.state.conditions.energy[bus.tag][electric_rotor.tag][rotor.tag]     
    #rotor_conditions.omega[:,0]  = test_omega    
    

    #state.conditions.energy[disributor.tag][propulsor.tag][rotor.tag].omega
    
    # Find relation between RPM and Thrust and Power for the two types of rotor  
    
    omega                            = np.array([500, 1000, 1500, 2000, 2500])
    
    thrust_prop                      = Data()
    thrust_lift                      = Data()
    power_prop                       = Data()
    power_lift                       = Data()
    
    for i in range(len(omega)):
        thrust_prop[i], power_prop[i] = compute_rotor_performance(propulsor_prop,state,disributor)
        thrust_lift[i], power_lift[i] = compute_rotor_performance(propulsor_lift,state,disributor)
    
    # Create the functions describing the thrust as function of RPM for both rotor types
    
    thrust_fz_prop            = interpolate.interp1d(omega, thrust_prop, kind='linear', fill_value="extrapolate", bounds_error=False)    
    power_fz_prop             = interpolate.interp1d(omega, power_prop, kind='linear', fill_value="extrapolate", bounds_error=False)
    thrust_fz_lift            = interpolate.interp1d(omega, thrust_lift, kind='linear', fill_value="extrapolate", bounds_error=False)
    power_fz_lift             = interpolate.interp1d(omega, power_lift, kind='linear', fill_value="extrapolate", bounds_error=False)   
 
    for i in range(len(CG_bat_1)):
        for j in range(len(CG_bat_2)):
            
            # -------------------------------------------------------------
            #  Assign new battery locations
            # -------------------------------------------------------------            
            
            # >>>>>>>>>>>>>>>>>>>>> ADD LATERAL CASES
            
            # prop rotor battery modules      
            case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_1[i]])
            case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_1[i,0] + case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.length, 
                                                                                                                 CG_bat_1[i,1], 
                                                                                                                 CG_bat_1[i,2]])
            # lift rotor battery modules 
            case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_2[j]])
            case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_2[j,0], 
                                                                                                                 CG_bat_2[j,1], 
                                                                                                                 CG_bat_2[j,2] + case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.height])
            # -------------------------------------------------------------
            #  Load specs from Stick Fixed Optimization 
            # -------------------------------------------------------------               
            
            excel_data         = read_results(excel_file)   # Change to for loop when add all results from stick fixed         
            
            # -------------------------------------------------------------
            #  Assign specs from Stick Fixed Optimization 
            # -------------------------------------------------------------            
            
            case_vehicle.wings.main_wing.spans.projected       = excel_data.mw_span[0]*10
            case_vehicle.wings.main_wing.aspect_ratio          = excel_data.mw_AR[0]*10
            case_vehicle.wings.main_wing.twists.root           = excel_data.mw_root_twist[0]
            case_vehicle.wings.main_wing.twists.tip            = excel_data.mw_tip_twist[0]
            case_vehicle.wings.vertical_tail.spans.projected   = excel_data.vt_span[0]
            case_vehicle.wings.vertical_tail.aspect_ratio      = excel_data.vt_AR[0]
            case_vehicle.wings.horizontal_tail.aspect_ratio    = excel_data.ht_AR[0]*10
            case_vehicle.wings.horizontal_tail.spans.projected = excel_data.ht_span[0]*10
            case_vehicle.wings.horizontal_tail.twists.tip      = excel_data.ht_tip_twist[0]
            # ADD control surfaces
            
            static_variable = Data()
            
            # Import data from vehicle
            
            # General 

            static_variable.W                                          = case_vehicle.mass_properties.max_takeoff*9.81
            static_variable.rho                                        = 1.225
            static_variable.S                                          = case_vehicle.wings.main_wing.areas.reference           
            static_variable.wing_span                                  = case_vehicle.wings.main_wing.spans.projected 
            static_variable.vertical_tail_span                         = case_vehicle.wings.vertical_tail.spans.projected 
            static_variable.C_Y_beta                                   = derivatives.CY_beta
            static_variable.C_L_beta                                   = derivatives.CL_beta
            static_variable.C_N_beta                                   = derivatives.CN_beta 
            static_variable.V_stall                                    = 48.624
            static_variable.V                                          = 1.13 *static_variable.V_stall
            static_variable.C_w                                        = static_variable.W/(0.5*static_variable.rho*(static_variable.V**2)*static_variable.S) 
            static_variable.aileron_chord_fraction                     = case_vehicle.wings.main_wing.control_surfaces.aileron.chord_fraction
            static_variable.rudder_chord_fraction                      = case_vehicle.wings.vertical_tail.control_surfaces.rudder.chord_fraction
            static_variable.wing_tip_chord                             = case_vehicle.wings.main_wing.chords.tip  
            static_variable.wing_root_chord                            = case_vehicle.wings.main_wing.chords.root
            static_variable.vertical_tail_tip_chord                    = case_vehicle.wings.vertical_tail.chords.tip 
            static_variable.vertical_tail_root_chord                   = case_vehicle.wings.vertical_tail.chords.root

            # HOVER OEI
                                                                         
            static_variable.thrust_fz_prop                             = thrust_fz_prop                                                                               
            static_variable.thrust_fz_lift                             = thrust_fz_lift
            static_variable.power_fz_prop                              = power_fz_prop                                                                               
            static_variable.power_fz_lift                              = power_fz_lift
            static_variable.arm1                                       = 0 # FIX
            static_variable.arm2                                       = 0 # FIX            
            static_variable.lat1                                       = 0 # FIX
            static_variable.lat2                                       = 0 # FIX
            static_variable.lat3                                       = 0 # FIX            
            
            RPM_1[i,j], RPM_2[i,j], RPM_3[i,j], RPM_4[i,j] = optimization(static_variable)

            debug = 0
    
    return RPM_1, RPM_2, RPM_3, RPM_4

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
    
def optimization(static_variable):
    
    # define optimizer bounds 
    RPM_1_lower_bound        = 100                # [RPM]  
    RPM_1_upper_bound        = 5000               # [RPM]
    RPM_2_lower_bound        = 100                # [RPM]  
    RPM_2_upper_bound        = 5000               # [RPM]
    RPM_3_lower_bound        = 100                # [RPM]  
    RPM_3_upper_bound        = 5000               # [RPM]
    RPM_4_lower_bound        = 100                # [RPM]  
    RPM_4_upper_bound        = 5000               # [RPM]    
    
   
    x0         = [2500, 2500, 2500, 2500]
    
    args       = [static_variable]
    
    hard_cons  = [{'type':'eq', 'fun': hard_constraint_hover_oei_Fz,'args': args},
                  {'type':'eq', 'fun': hard_constraint_hover_oei_Mx,'args': args},
                  {'type':'eq', 'fun': hard_constraint_hover_oei_Mz,'args': args}]
    
    bnds       = ((RPM_1_lower_bound,        RPM_1_upper_bound),                 
                  (RPM_2_lower_bound,        RPM_2_upper_bound),                 
                  (RPM_3_lower_bound,        RPM_3_upper_bound),                 
                  (RPM_4_lower_bound,        RPM_4_upper_bound)) 
    
    # try hard constraints to find optimum motor parameters
    sol = minimize(objective, x0, args= (static_variable) , method='trust-constr', bounds=bnds, tol=1e-6, constraints=hard_cons) 
    
    if sol.success == False:
        print('\n Optimum control surfaces design failed.')
        print(sol)
    
    print(sol)
    
    RPM_1_hover_oei     = sol.x[0]
    RPM_2_hover_oei     = sol.x[1]   
    RPM_3_hover_oei     = sol.x[2] 
    RPM_4_hover_oei     = sol.x[3]   
    
    return RPM_1_hover_oei, RPM_2_hover_oei, RPM_3_hover_oei, RPM_4_hover_oei

# objective function
def objective(x,static_variable):
    
    Power = 3*static_variable.power_fz_prop(x[7]) + 3*static_variable.power_fz_prop(x[8]) + 3*static_variable.power_fz_lift(x[9]) + 2*static_variable.power_fz_lift(x[10])
    
    print(Power)
        
    return Power
    
def hard_constraint_hover_oei_Fz(x,static_variable):     

    W               = static_variable.W   
    
    T1              = static_variable.thrust_fz_prop(x[0])
    T2              = static_variable.thrust_fz_prop(x[0])
    T3              = static_variable.thrust_fz_prop(x[0])
    T4              = static_variable.thrust_fz_prop(x[1])
    T5              = static_variable.thrust_fz_prop(x[1])
    T6              = static_variable.thrust_fz_prop(x[1])
    
    T7              = static_variable.thrust_fz_lift(x[2])
    T8              = static_variable.thrust_fz_lift(x[2])
    T9              = static_variable.thrust_fz_lift(x[2])
    T10             = static_variable.thrust_fz_lift(x[3])
    T11             = static_variable.thrust_fz_lift(x[3])
    
    T               = T1 + T2 + T3 + T4 + T5 + T6 + T7 + T8 + T9 + T10 + T11
    
    res             = T - W    
    
    return res

def hard_constraint_hover_oei_Mx(x,static_variable):     

    arm1            = static_variable.arm1
    arm2            = static_variable.arm2

    T1              = static_variable.thrust_fz_prop(x[0])
    T2              = static_variable.thrust_fz_prop(x[0])
    T3              = static_variable.thrust_fz_prop(x[0])
    T4              = static_variable.thrust_fz_prop(x[1])
    T5              = static_variable.thrust_fz_prop(x[1])
    T6              = static_variable.thrust_fz_prop(x[1])
    
    T7              = static_variable.thrust_fz_lift(x[2])
    T8              = static_variable.thrust_fz_lift(x[2])
    T9              = static_variable.thrust_fz_lift(x[2])
    T10             = static_variable.thrust_fz_lift(x[3])
    T11             = static_variable.thrust_fz_lift(x[3])    

    res             = arm1*(T1 + T2 + T3 + T4 + T5 + T6) - arm2*(T7 + T8 + T9 + T10 + T11)    
    
    return res

def hard_constraint_hover_oei_Mz(x,static_variable):     

    lat1            = static_variable.lat1
    lat2            = static_variable.lat2
    lat3            = static_variable.lat3

    T1              = static_variable.thrust_fz_prop(x[0])
    T2              = static_variable.thrust_fz_prop(x[0])
    T3              = static_variable.thrust_fz_prop(x[0])
    T4              = static_variable.thrust_fz_prop(x[1])
    T5              = static_variable.thrust_fz_prop(x[1])
    T6              = static_variable.thrust_fz_prop(x[1])
    
    T7              = static_variable.thrust_fz_lift(x[2])
    T8              = static_variable.thrust_fz_lift(x[2])
    T9              = static_variable.thrust_fz_lift(x[2])
    T10             = static_variable.thrust_fz_lift(x[3])
    T11             = static_variable.thrust_fz_lift(x[3])    

    res             = lat1*(T4 + T9 - T10 - T3) + lat2*(T8 + T5 - T11 - T2) + lat3*(T7 + T6 - T1)    
    
    return res  

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