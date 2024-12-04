
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from   RCAIDE.Framework.Core import Units,Data  
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
    cruise_velocity = 150  * Units['mph']
    cruise_altitude = 1000*Units.feet
    
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
    
    # delete control surfaces if they have been defined 
    for wing in case_vehicle.wings:
        for control_surface in wing.control_surfaces:
            del wing.control_surfaces[control_surface.tag]
            
    # -------------------------------------------------------------
    #  Sweeps of battery locations (UPDATE WITH LAT VALUES LATER)
    # -------------------------------------------------------------            
    
    # prop rotor battery module (first module)
    #                 CG: X,    Y,    Z 
    CG_bat_1 = np.array([[0.25, 0.,   0.],
                         [0.35, 0.,   0.],
                         [0.45, 0.,   0.],
                         [0.25, 0.25, 0.],
                         [0.25, 0.35, 0.],
                         [0.25, 0.45, 0.]])
    
    # lift rotor battery modules
    #                 CG: X,    Y,    Z 
    CG_bat_2 = np.array([[4.0,  0.,   0.],
                         [4.1,  0.,   0.],
                         [4.2,  0.,   0.]])   
 
    for i in range(len(CG_bat_1)):
        for j in range(len(CG_bat_2)):
            
            # -------------------------------------------------------------
            #  Assign new battery locations
            # -------------------------------------------------------------            
            
            # prop rotor battery modules      
            case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_1[i]])
            case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_1[i,0] + case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.length, 
                                                                                                                 CG_bat_1[i,1], 
                                                                                                                 CG_bat_1[i,2]])
            # lift rotor battery modules 
            case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_2[j]])
            case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_2[i,0], 
                                                                                                                 CG_bat_2[i,1], 
                                                                                                                 CG_bat_2[i,2] + case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.height])
            
            # -------------------------------------------------------------
            #  Load specs from Stick Fixed Optimization 
            # -------------------------------------------------------------               
            
            excel_data         = read_results(excel_file)            
            
            # -------------------------------------------------------------
            #  Assign specs from Stick Fixed Optimization 
            # -------------------------------------------------------------            
            
            #case_vehicle.wings.main_wing.spans.projected = excel_data.mw_span
            #case_vehicle.wings.main_wing            = excel_data.mw_AR
            #case_vehicle.wings.main_wing            = excel_data.mw_root_twist
            #case_vehicle.wings.main_wing            = excel_data.mw_tip_twist
            #case_vehicle.wings.vertical_tail         = excel_data.vt_span
            #case_vehicle.wings.vertical_tail         = excel_data.vt_AR
            #case_vehicle.wings.horizontal_tail       = excel_data.ht_AR
            #case_vehicle.wings.horizontal_tail       = excel_data.ht_span
            #case_vehicle.wings.horizontal_tail       = excel_data.ht_tip_twist
            
            # -------------------------------------------------------------
            #  Aileron & Rudder
            # -------------------------------------------------------------
            
            # Create arrays of aileron and rudder sizes
            
            aileron_size                              = np.array([0.1, 0.5, 0.9])
            rudder_size                               = np.array([0.1, 0.5, 0.9])
            
            # Append aileron and rudder
            
            case_vehicle                              = setup_rudder_aileron(case_vehicle)
            
            # Compute the control derivatives related to the control surfaces
            
            derivatives, results                      = compute_rudder_aileron_derivatives(aileron_size, rudder_size, case_vehicle, seg_num=0)
            
            C_Y_delta_a                               = derivatives.C_Y_delta_a
            C_L_delta_a                               = derivatives.C_L_delta_a
            C_N_delta_a                               = derivatives.C_N_delta_a
            C_Y_delta_r                               = derivatives.C_Y_delta_r
            C_L_delta_r                               = derivatives.C_L_delta_r
            C_N_delta_r                               = derivatives.C_N_delta_r          
            
            # Create the functions describing the control derivatives as functions of control surface size
            
            static_variable = Data()
            
            static_variable.C_Y_delta_a_fz            = interpolate.interp1d(aileron_size, C_Y_delta_a)
            static_variable.C_L_delta_a_fz            = interpolate.interp1d(aileron_size, C_L_delta_a)
            static_variable.C_N_delta_a_fz            = interpolate.interp1d(aileron_size, C_N_delta_a)   
            static_variable.C_Y_delta_r_fz            = interpolate.interp1d(rudder_size,  C_Y_delta_r)
            static_variable.C_L_delta_r_fz            = interpolate.interp1d(rudder_size,  C_L_delta_r)
            static_variable.C_N_delta_r_fz            = interpolate.interp1d(rudder_size,  C_N_delta_r)
            
            # Import data from vehicle
            
            # General 
            
            static_variable.W                                          = 2700*9.81
            static_variable.rho                                        = 1.225
            static_variable.S                                          = 15.629           
            static_variable.wing_span                                  = 11.82855 
            static_variable.vertical_tail_span                         = 2
            static_variable.aileron_chord                              = 0.1
            static_variable.rudder_chord                               = 0.1  
            static_variable.C_Y_beta                                   = 1 
            static_variable.C_L_beta                                   = 1
            static_variable.C_N_beta                                   = 1 
            static_variable.V                                          = 110. * Units['mph']
            static_variable.C_w                                        = static_variable.W/(0.5*static_variable.rho*(static_variable.V**2)*static_variable.S)    
            
            # CROSSWIND 
         
            static_variable.C_Y_0_cw                                   = 0                       
            static_variable.C_L_0_cw                                   = 0    
            static_variable.C_N_0_cw                                   = 0  
            static_variable.v_cw                                       = 35  * Units['ft/s']                          # [ft/s] 
            static_variable.beta_cw                                    = np.arcsin(static_variable.v_cw/static_variable.V) # [rad]  
            static_variable.phi_cw                                     = 0                                                                    #[-]
            
            # OEI
            
            static_variable.T_oei                                      = 525              # [N]
            static_variable.arm_oei                                    = 9.2              # [m]
            static_variable.beta_oei                                   = 0 * np.pi/180  
            static_variable.C_Y_0_oei                                  = 0                        
            static_variable.C_L_0_oei                                  = 0                                                        
            static_variable.C_N_0_oei                                  = -(static_variable.T_oei*static_variable.arm_oei)/(0.5*static_variable.rho*(static_variable.V**2)*static_variable.S*static_variable.wing_span)
            static_variable.phi_oei                                    = 0
            
            delta_a_lower_bound_cw, delta_r_lower_bound_cw, delta_a_lower_bound_oei, delta_r_lower_bound_oei, aileron_span_lower_bound, rudder_span_lower_bound = optimization(static_variable)
    
    return delta_a_size, delta_r_size 

def interpolation(x, y, kind='linear'):
    # ------------------------------------------------------------------
    #   Interpolation
    # ------------------------------------------------------------------
    
    interpolating_function = interp1d(x, y, kind=kind, fill_value="extrapolate")
    
    x_new = np.linspace(x.min(), x.max(), 100)
    y_new = interpolating_function(x_new)
    
    if kind == 'linear':
        dy_dx = np.gradient(y_new, x_new)
    else:
        dy_dx = np.gradient(y_new, x_new)
    
    return interpolating_function, x_new, y_new, dy_dx

def setup_rudder_aileron(vehicle):
    
    mw_wing                       = vehicle.wings.main_wing 
    aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.9 
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.2
    mw_wing.append_control_surface(aileron)     
    
    vt_wing                      = vehicle.wings.vertical_tail
    rudder                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
    rudder.tag                   = 'rudder'
    rudder.span_fraction_start   = 0.2
    rudder.span_fraction_end     = 0.8
    rudder.deflection            = 0.0  * Units.deg
    rudder.chord_fraction        = 0.2
    vt_wing.append_control_surface(rudder) 
    
    return vehicle

def compute_rudder_aileron_derivatives(aileron_size, rudder_size, vehicle,  seg_num=0):           
    # seg_num is the number of the mission segment correspdoning to the segment for which control surface derivatives are being found
    
    CN_delta_a    =  np.zeros(np.size(aileron_size))
    CL_delta_a    =  np.zeros(np.size(aileron_size))
    CY_delta_a    =  np.zeros(np.size(aileron_size))
    
    CN_delta_r     = np.zeros(np.size(rudder_size))
    CL_delta_r     = np.zeros(np.size(rudder_size))
    CY_delta_r     = np.zeros(np.size(rudder_size))
    
    for i in range(len(aileron_size)):
        vehicle.wings.main_wing.control_surfaces.aileron.chord_fraction = aileron_size[i]
        results =  evalaute_aircraft(vehicle)
            
        # store properties
        CN_delta_a[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CN_delta_a[0, 0]
        CL_delta_a[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CL_delta_a[0, 0]
        CY_delta_a[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CY_delta_a[0, 0]
        
    for i in range(len(rudder_size)):
        vehicle.wings.vertical_tail.rudder.chord_fraction = rudder_size[i]
        results =  evalaute_aircraft(vehicle)
            
        # store properties
        CN_delta_r[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CN_delta_r[0, 0]
        CL_delta_r[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CL_delta_r[0, 0]
        CY_delta_r[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CY_delta_r[0, 0]
        
    # Pack data for return 
    derivatives = {'delta_r':{'CN_delta_r': CN_delta_r,'CL_delta_r': CL_delta_r, 'CY_delta_r': CY_delta_r}, 'delta_a':{'CN_delta_a': CN_delta_a,'CL_delta_a': CL_delta_a, 'CY_delta_a': CY_delta_a}}
                
    return derivatives, results
                 
                
def evalaute_aircraft(vehicle): 
        
    # Set up configs
    configs  = configs_setup(vehicle)

    # vehicle analyses
    analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(analyses)
    missions = missions_setup(mission) 
     
    results = missions.base_mission.evaluate()
    
          
    return results
 
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
    analyses        = RCAIDE.Framework.Analyses.Vehicle() 

    # ------------------------------------------------------------------
    #  Weights
    # ------------------------------------------------------------------
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    # ------------------------------------------------------------------
    aerodynamics                                      = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle                              = vehicle  
    aerodynamics.settings.use_surrogate               = False 
    aerodynamics.settings.trim_aircraft               = False 
    analyses.append(aerodynamics)
    
    # ------------------------------------------------------------------
    #  Stability Analysis
    stability         = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method() 
    stability.vehicle = vehicle
    analyses.append(stability)
    
    # ------------------------------------------------------------------
    #  Energy
    # ------------------------------------------------------------------
    energy     = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle  
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    # ------------------------------------------------------------------
    planet     = RCAIDE.Framework.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    # ------------------------------------------------------------------
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses

  
def configs_setup(vehicle): 
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
 
def base_mission_setup(analyses,max_speed_multiplier,cruise_velocity,cruise_altitude,angle_of_attack,sideslip_angle,bank_angle, roll_rate,pitch_rate,yaw_rate):   
    '''
    This sets up the nominal cruise of the aircraft
    '''
     
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'base_mission'
  
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    #   Cruise Segment: constant Speed, constant altitude 
    segment                           = Segments.Untrimmed.Untrimmed()
    segment.analyses.extend( analyses )   
    segment.tag                       = "cruise"
    segment.angle_of_attack           = angle_of_attack
    segment.sideslip_angle            = sideslip_angle
    segment.bank_angle                = bank_angle
    segment.altitude                  = cruise_altitude
    segment.air_speed                 = cruise_velocity * max_speed_multiplier
    segment.roll_rate                 = roll_rate 
    segment.pitch_rate                = pitch_rate  
    segment.yaw_rate                  = yaw_rate   

    segment.flight_dynamics.force_x   = True    
    segment.flight_dynamics.force_z   = True    
    segment.flight_dynamics.force_y   = True     
    segment.flight_dynamics.moment_y  = True 
    segment.flight_dynamics.moment_x  = True
    segment.flight_dynamics.moment_z  = True
    
    mission.append_segment(segment)     
    
    return mission

def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions  

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
    aileron_span_lower_bound = 0.75               # [%]    # just changing the upper and lower bound of the lower bound of the control surface (the inner limit closer to root of wing)
    aileron_span_upper_bound = 0.90               # [%]  
    rudder_span_lower_bound  = 0.05               # [%]  
    rudder_span_upper_bound  = 0.90               # [%]
    delta_a_lower_bound_cw   = -30 * np.pi/180    # [rad]
    delta_a_upper_bound_cw   = 30  * np.pi/180    # [rad]
    delta_r_lower_bound_cw   = -30 * np.pi/180    # [rad]
    delta_r_upper_bound_cw   = 30  * np.pi/180    # [rad]                     
    delta_a_lower_bound_oei  = -30 * np.pi/180    # [rad]
    delta_a_upper_bound_oei  = 30  * np.pi/180    # [rad]
    delta_r_lower_bound_oei  = -30 * np.pi/180    # [rad]
    delta_r_upper_bound_oei  = 30  * np.pi/180    # [rad]   
    
    x0         = [0.3, 0.1, 0.3, 0.1, 0.3, 0.3]
    
    args       = [static_variable]
    
    hard_cons  = [{'type':'eq', 'fun': hard_constraint_Y_cw,'args':  args},
                  {'type':'eq', 'fun': hard_constraint_L_cw,'args':  args},
                  {'type':'eq', 'fun': hard_constraint_N_cw,'args':  args},
                  {'type':'eq', 'fun': hard_constraint_Y_oei,'args': args},
                  {'type':'eq', 'fun': hard_constraint_L_oei,'args': args},
                  {'type':'eq', 'fun': hard_constraint_N_oei,'args': args}]
    
    bnds       = ((delta_a_lower_bound_cw,   delta_a_upper_bound_cw  ), 
                  (delta_r_lower_bound_cw ,  delta_r_upper_bound_cw  ),
                  (delta_a_lower_bound_oei,  delta_a_upper_bound_oei ),
                  (delta_r_lower_bound_oei,  delta_r_upper_bound_oei ),
                  (aileron_span_lower_bound, aileron_span_upper_bound),                 
                  (rudder_span_lower_bound,  rudder_span_upper_bound )) 
    
    # try hard constraints to find optimum motor parameters
    sol = minimize(objective, x0, args= (static_variable) , method='SLSQP', bounds=bnds, tol=1e-6, constraints=hard_cons) 
    
    if sol.success == False:
        print('\n Optimum control surfaces design failed.')
    
    delta_a_cw   = sol.x[0]
    delta_r_cw   = sol.x[1]    
    delta_a_oei  = sol.x[2] 
    delta_r_oei  = sol.x[3] 
    aileron_span = sol.x[4] 
    rudder_span  = sol.x[5] 
    
    return delta_a_cw, delta_r_cw, delta_a_oei, delta_r_oei, aileron_span, rudder_span 

# objective function
def objective(design_variables,static_variable):
    aileron_span    = static_variable.wing_span*(0.95 - design_variables[4])
    rudder_span     = static_variable.vertical_tail_span*(0.95 - design_variables[5]) 
    
    return (aileron_span*static_variable.aileron_chord)*(rudder_span*static_variable.rudder_chord)

# hard constraint
def hard_constraint_Y_cw(design_variables,static_variable):     
    
    aileron_span    = static_variable.wing_span*(0.95 - design_variables[4])
    rudder_span     = static_variable.vertical_tail_span*(0.95 - design_variables[5])
    
    C_w             = static_variable.C_w
    phi_cw          = static_variable.phi_cw
    C_Y_0_cw        = static_variable.C_Y_0_cw
    C_Y_beta        = static_variable.C_Y_beta  
    beta_cw         = static_variable.beta_cw 
    C_Y_delta_a_fz  = static_variable.C_Y_delta_a_fz
    C_Y_delta_r_fz  = static_variable.C_Y_delta_r_fz
    
    res = C_w*np.sin(phi_cw) + C_Y_0_cw + C_Y_beta*beta_cw + C_Y_delta_a_fz(aileron_span)*design_variables[0] + C_Y_delta_r_fz(rudder_span)*design_variables[1]
    
    return res

def hard_constraint_L_cw(design_variables,static_variable): 
    
    aileron_span    = static_variable.wing_span*(0.95 - design_variables[4])
    rudder_span     = static_variable.vertical_tail_span*(0.95 - design_variables[5])
    
    C_L_0_cw        = static_variable.C_L_0_cw
    C_L_beta        = static_variable.C_L_beta 
    beta_cw         = static_variable.beta_cw 
    C_L_delta_a_fz  = static_variable.C_L_delta_a_fz
    C_L_delta_r_fz  = static_variable.C_L_delta_r_fz
    
    res = C_L_0_cw + C_L_beta*beta_cw + C_L_delta_a_fz(aileron_span)*design_variables[0] +  C_L_delta_r_fz(rudder_span)*design_variables[1]
    
    return res

def hard_constraint_N_cw(design_variables,static_variable): 
    
    aileron_span    = static_variable.wing_span*(0.95 - design_variables[4])
    rudder_span     = static_variable.vertical_tail_span*(0.95 - design_variables[5])
    
    C_N_0_cw        = static_variable.C_N_0_cw
    C_N_beta        = static_variable.C_N_beta
    beta_cw         = static_variable.beta_cw 
    C_N_delta_a_fz  = static_variable.C_N_delta_a_fz
    C_N_delta_r_fz  = static_variable.C_N_delta_r_fz
    
    res = C_N_0_cw + C_N_beta*beta_cw + C_N_delta_a_fz(aileron_span)*design_variables[0] +  C_N_delta_r_fz(rudder_span)*design_variables[1]
    
    return res

def hard_constraint_Y_oei(design_variables,static_variable): 
    
    aileron_span     = static_variable.wing_span*(0.95 - design_variables[4])
    rudder_span      = static_variable.vertical_tail_span*(0.95 - design_variables[5])
    
    C_w              = static_variable.C_w
    phi_oei          = static_variable.phi_oei
    C_Y_0_oei        = static_variable.C_Y_0_oei
    C_Y_beta         = static_variable.C_Y_beta 
    beta_oei         = static_variable.beta_oei     
    C_Y_delta_a_fz   = static_variable.C_Y_delta_a_fz
    C_Y_delta_r_fz   = static_variable.C_Y_delta_r_fz
    
    res = C_w*np.sin(phi_oei) + C_Y_0_oei + C_Y_beta*beta_oei + C_Y_delta_a_fz(aileron_span)*design_variables[2] + C_Y_delta_r_fz(rudder_span)*design_variables[3]
    
    return res

def hard_constraint_L_oei(design_variables,static_variable): 
    
    aileron_span     = static_variable.wing_span*(0.95 - design_variables[4])
    rudder_span      = static_variable.vertical_tail_span*(0.95 - design_variables[5])
    
    C_L_0_oei        = static_variable.C_L_0_oei
    C_L_beta         = static_variable.C_L_beta 
    beta_oei         = static_variable.beta_oei 
    C_L_delta_a_fz   = static_variable.C_L_delta_a_fz
    C_L_delta_r_fz   = static_variable.C_L_delta_r_fz  
    
    res = C_L_0_oei + C_L_beta*beta_oei + C_L_delta_a_fz(aileron_span)*design_variables[2] +  C_L_delta_r_fz(rudder_span)*design_variables[3]
    
    return res

def hard_constraint_N_oei(design_variables,static_variable): 
    
    aileron_span     = static_variable.wing_span*(0.95 - design_variables[4])
    rudder_span      = static_variable.vertical_tail_span*(0.95 - design_variables[5])
    
    C_N_0_oei        = static_variable.C_N_0_oei
    C_N_beta         = static_variable.C_N_beta
    beta_oei         = static_variable.beta_oei 
    C_N_delta_a_fz   = static_variable.C_N_delta_a_fz
    C_N_delta_r_fz   = static_variable.C_N_delta_r_fz 
    
    res = C_N_0_oei + C_N_beta*beta_oei + C_N_delta_a_fz(aileron_span)*design_variables[2] +  C_N_delta_r_fz(rudder_span)*design_variables[3]
    
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