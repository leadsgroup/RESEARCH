
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from   RCAIDE.Framework.Core import Units,Data  
from   RCAIDE.Library.Plots  import *  
import numpy                 as np
import  scipy as sp
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
    filenames     =  ["025_00_00_40_00_00_Baseline_stick_fixed_cruise_Opt_Results.xlsx"]
    
    
    # delete control surfaces if they have been defined 
    for wing in case_vehicle.wings:
        for control_surface in wing.control_surfaces:
            del wing.control_surfaces[control_surface.tag]
            
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
    
    delta_a_cw   = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    delta_r_cw   = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    delta_a_oei  = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    delta_r_oei  = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    delta_a_roll = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    aileron_span = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    rudder_span  = np.zeros([len(CG_bat_1),len(CG_bat_2)])
    
    # Create arrays of aileron and rudder  
    aileron_span_fraction_start               = np.array([0.6, 0.7, 0.8])  # inboard aileron span percentage
    rudder_span_fraction_start                = np.array([0.2, 0.3, 0.4])  # inboard rudder span percentage
    
    # Append aileron and rudder 
    case_vehicle  = append_rudder_and_aileron(case_vehicle)    
   
    idx =  0
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
        
         
            separator          = os.path.sep
            relative_path      = sys.path[0] + separator        
            excel_file         = relative_path + filenames[idx]         
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
            
            # -------------------------------------------------------------
            #  Aileron & Rudder
            # -------------------------------------------------------------
            
            # Compute the control derivatives related to the control surface dimensions
            
            derivatives, results                      = compute_rudder_aileron_derivatives(aileron_span_fraction_start, rudder_span_fraction_start, case_vehicle, seg_num=0)
            
            CY_delta_a                                = derivatives.CY_delta_a
            CL_delta_a                                = derivatives.CL_delta_a
            CN_delta_a                                = derivatives.CN_delta_a
            CY_delta_r                                = derivatives.CY_delta_r
            CL_delta_r                                = derivatives.CL_delta_r
            CN_delta_r                                = derivatives.CN_delta_r          
            
            # Create the functions describing the control derivatives as functions of control surface size
            
            static_variable = Data()
            
            static_variable.C_Y_delta_a_fz            = interpolate.interp1d(aileron_span_fraction_start, CY_delta_a, kind='linear', fill_value="extrapolate", bounds_error=False)
            static_variable.C_L_delta_a_fz            = interpolate.interp1d(aileron_span_fraction_start, CL_delta_a, kind='linear', fill_value="extrapolate", bounds_error=False)
            static_variable.C_N_delta_a_fz            = interpolate.interp1d(aileron_span_fraction_start, CN_delta_a, kind='linear', fill_value="extrapolate", bounds_error=False)   
            static_variable.C_Y_delta_r_fz            = interpolate.interp1d(rudder_span_fraction_start,  CY_delta_r, kind='linear', fill_value="extrapolate", bounds_error=False)
            static_variable.C_L_delta_r_fz            = interpolate.interp1d(rudder_span_fraction_start,  CL_delta_r, kind='linear', fill_value="extrapolate", bounds_error=False)
            static_variable.C_N_delta_r_fz            = interpolate.interp1d(rudder_span_fraction_start,  CN_delta_r, kind='linear', fill_value="extrapolate", bounds_error=False)
            
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
                        
            # CROSSWIND 
         
            static_variable.C_Y_0_cw                                   = 0                       
            static_variable.C_L_0_cw                                   = 0    
            static_variable.C_N_0_cw                                   = 0  
            static_variable.v_cw                                       = 35  * Units['ft/s']                          # [ft/s] 
            static_variable.beta_cw                                    = np.arcsin(static_variable.v_cw/static_variable.V) # [rad]  
            static_variable.phi_cw                                     = 5 * Units.degrees                                        #[-]
            
            # OEI
            
            static_variable.T_oei                                      = 525              # [N]
            static_variable.arm_oei                                    = case_vehicle.networks.electric.busses.prop_rotor_bus.propulsors.prop_rotor_propulsor_6.rotor.origin[0][1]              # [m]
            static_variable.beta_oei                                   = 0 * Units.degrees 
            static_variable.C_Y_0_oei                                  = 0                        
            static_variable.C_L_0_oei                                  = 0                                                        
            static_variable.C_N_0_oei                                  = -(static_variable.T_oei*static_variable.arm_oei)/(0.5*static_variable.rho*(static_variable.V**2)*static_variable.S*static_variable.wing_span)
            static_variable.phi_oei                                    = 5 * Units.degrees 
            
            # AILERON RATE OF ROLL
            
            static_variable.CS_angle                                   = 60* Units.degrees # [rad]
            static_variable.CS_time                                    = 5             # [sec]
            static_variable.CD_R                                       = 0.9          
            static_variable.Ixx                                        = case_vehicle.mass_properties.moments_of_inertia.tensor[0,0]
            static_variable.S_w                                        = case_vehicle.wings.main_wing.areas.reference
            static_variable.S_h                                        = case_vehicle.wings.horizontal_tail.areas.reference
            static_variable.S_v                                        = case_vehicle.wings.vertical_tail.areas.reference
            
            delta_a_cw[i,j], delta_r_cw[i,j], delta_a_oei[i,j], delta_r_oei[i,j],delta_a_roll[i,j], aileron_span[i,j], rudder_span[i,j] = optimization(static_variable)
            
            idx += 1
    
    return delta_a_cw, delta_r_cw, delta_a_oei, delta_r_oei, aileron_span, rudder_span

def append_rudder_and_aileron(vehicle):
    
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

def compute_rudder_aileron_derivatives(aileron_span_fraction_start, rudder_span_fraction_start, vehicle,  seg_num=0):           
    # seg_num is the number of the mission segment correspdoning to the segment for which control surface derivatives are being found
    
    CN_delta_a    =  np.zeros(np.size(aileron_span_fraction_start))
    CL_delta_a    =  np.zeros(np.size(aileron_span_fraction_start))
    CY_delta_a    =  np.zeros(np.size(aileron_span_fraction_start))
    
    CN_delta_r     = np.zeros(np.size(rudder_span_fraction_start))
    CL_delta_r     = np.zeros(np.size(rudder_span_fraction_start))
    CY_delta_r     = np.zeros(np.size(rudder_span_fraction_start))
    
    for i in range(len(aileron_span_fraction_start)):
        vehicle.wings.main_wing.control_surfaces['aileron'].span_fraction_start = aileron_span_fraction_start[i]
        results =  evalaute_aircraft(vehicle)
            
        # store properties
        CN_delta_a[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CN_delta_a[0, 0]
        CL_delta_a[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CL_delta_a[0, 0]
        CY_delta_a[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CY_delta_a[0, 0]
        
    for i in range(len(rudder_span_fraction_start)):
        vehicle.wings.vertical_tail.control_surfaces['rudder'].span_fraction_start = rudder_span_fraction_start[i]
        results =  evalaute_aircraft(vehicle)
            
        # store properties
        CN_delta_r[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CN_delta_r[0, 0]
        CL_delta_r[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CL_delta_r[0, 0]
        CY_delta_r[i]     =  results.segments[seg_num].conditions.static_stability.derivatives.CY_delta_r[0, 0]
        
    CN_beta =  results.segments[seg_num].conditions.static_stability.derivatives.CN_beta[0, 0]
    CL_beta =  results.segments[seg_num].conditions.static_stability.derivatives.CL_beta[0, 0]
    CY_beta =  results.segments[seg_num].conditions.static_stability.derivatives.CY_beta[0, 0]    
        
    # Pack data for return 
    derivatives = Data(  
         CN_delta_r =  CN_delta_r,
         CL_delta_r = CL_delta_r,
         CY_delta_r =  CY_delta_r , 
         CN_delta_a =  CN_delta_a,
         CL_delta_a =  CL_delta_a, 
         CY_delta_a =  CY_delta_a,
         CN_beta    =  CN_beta,
         CL_beta    =  CL_beta,
         CY_beta    =  CY_beta) 
                
    return derivatives, results
                 
                
def evalaute_aircraft(vehicle): 
    
    # Set up vehicle configs
    configs  = configs_setup(vehicle)

    # create analyses
    analyses = analyses_setup(configs)

    # mission analyses
    mission  = base_mission_setup(analyses) 

    # create mission instances (for multiple types of missions)
    missions = missions_setup(mission) 

    # mission analysis 
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
    #  Energy
    # ------------------------------------------------------------------
    energy     = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle  
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    # ------------------------------------------------------------------
    planet     = RCAIDE.Framework.Analyses.Planets.Earth()
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
 
def base_mission_setup(analyses):   
    '''
    This sets up the nominal cruise of the aircraft
    '''
     
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'base_mission'
  
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    #   Cruise Segment: constant Speed, constant altitude 
    segment                           = Segments.Untrimmed.Untrimmed()
    segment.analyses.extend( analyses.forward_flight )   
    segment.tag                       = "cruise" 
    segment.altitude                  = 1000 * Units.feet 
    segment.air_speed                 = 150  * Units['mph']

    segment.flight_dynamics.force_x   = True    
    segment.flight_dynamics.force_z   = True    
    segment.flight_dynamics.force_y   = True     
    segment.flight_dynamics.moment_y  = True 
    segment.flight_dynamics.moment_x  = True
    segment.flight_dynamics.moment_z  = True
    
    mission.append_segment(segment)     
    
    return mission

def missions_setup(mission): 
 
    missions     = RCAIDE.Framework.Mission.Missions()
    
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
    aileron_span_lower_bound = 0.5                    # [%]    # just changing the upper and lower bound of the lower bound of the control surface (the inner limit closer to root of wing)
    aileron_span_upper_bound = 0.90                   # [%]  
    rudder_span_lower_bound  = 0.05                   # [%]  
    rudder_span_upper_bound  = 0.90                   # [%]        
    delta_a_lower_bound_cw   = -30 * Units.degrees    # [rad]
    delta_a_upper_bound_cw   = 30  * Units.degrees    # [rad]
    delta_r_lower_bound_cw   = -30 * Units.degrees    # [rad]
    delta_r_upper_bound_cw   = 30  * Units.degrees    # [rad]                
    delta_a_lower_bound_oei  = -30 * Units.degrees    # [rad]
    delta_a_upper_bound_oei  = 30  * Units.degrees    # [rad]
    delta_r_lower_bound_oei  = -30 * Units.degrees    # [rad]
    delta_r_upper_bound_oei  = 30  * Units.degrees    # [rad]                    
    delta_a_lower_bound_roll = -30 * Units.degrees    # [rad]
    delta_a_upper_bound_roll = 30  * Units.degrees    # [rad]  
    
   
    x0          = [  0.9, 0.9,  0.5, 0.5] #, 0.5, 0.5, 0.5]
    
    args       = [static_variable]
    
    hard_cons  = [{'type':'eq', 'fun': con_Y_cw,'args':  args},
                  {'type':'eq', 'fun': con_L_cw,'args':  args},
                  {'type':'eq', 'fun': con_N_cw,'args':  args},
                  #{'type':'eq', 'fun': con_Y_oei,'args': args},
                  #{'type':'eq', 'fun': con_L_oei,'args': args},
                  #{'type':'eq', 'fun': con_N_oei,'args': args},
                  #{'type':'ineq', 'fun': con_roll_rate,'args': args}
                  ]
    
    bnds       = ( 
                  (aileron_span_lower_bound, aileron_span_upper_bound),                 
                  (rudder_span_lower_bound,  rudder_span_upper_bound ),
                  (delta_a_lower_bound_cw,   delta_a_upper_bound_cw  ), 
                  (delta_r_lower_bound_cw ,  delta_r_upper_bound_cw  ),
                  #(delta_a_lower_bound_oei,  delta_a_upper_bound_oei ),
                  #(delta_r_lower_bound_oei,  delta_r_upper_bound_oei ),
                  #(delta_a_lower_bound_roll,  delta_a_upper_bound_roll),
                  ) 
    
    # try hard constraints to find optimum motor parameters
    hard_cons  = ({'type':'eq', 'fun': con_Y_cw,'args':  args},
                  {'type':'eq', 'fun': con_L_cw,'args':  args},
                  {'type':'eq', 'fun': con_N_cw,'args':  args},
                  #{'type':'eq', 'fun': con_Y_oei,'args': args},
                  #{'type':'eq', 'fun': con_L_oei,'args': args},
                  #{'type':'eq', 'fun': con_N_oei,'args': args},
                  #{'type':'ineq', 'fun': con_roll_rate,'args': args}
                  )
        
    sol = sp.optimize.fmin_slsqp(objective,x0, args= static_variable ,f_eqcons=hard_cons,bounds=bnds, epsilon = 1e-8, acc  = 1e-6)
    
    #sol = minimize(objective, x0, args= (static_variable) , method='SLSQP', bounds=bnds, tol=1e-6, constraints=hard_cons)  
    print(sol)
        
    aileron_lower_bound = sol.x[0] 
    rudder_lower_bound  = sol.x[1] 
    delta_a_cw          = sol.x[2]
    delta_r_cw          = sol.x[3]   
    #delta_a_oei         = sol.x[4] 
    #delta_r_oei         = sol.x[5] 
    #delta_a_roll        = sol.x[6]  
        
    return aileron_lower_bound, rudder_lower_bound ,  delta_a_cw, delta_r_cw # , delta_a_oei, delta_r_oei , delta_a_roll 

# objective function
def objective(x,static_variable):
    
    aileron_span                     = static_variable.wing_span*(0.95 - x[0])
    delta_chord                      = static_variable.wing_root_chord - static_variable.wing_tip_chord 
    aileron_start_chord_wing_segment = static_variable.wing_root_chord - delta_chord * x[0]
    aileron_end_chord_wing_segment   = static_variable.wing_root_chord - delta_chord * 0.95
    aileron_start_chord              = static_variable.aileron_chord_fraction* aileron_start_chord_wing_segment 
    aileron_end_chord                = static_variable.aileron_chord_fraction* aileron_end_chord_wing_segment 
    Area_aileron                     = aileron_span  * ( aileron_start_chord +  aileron_end_chord) / 2
    

    rudder_span                     = static_variable.vertical_tail_span*(0.95 - x[1])
    delta_chord                     = static_variable.vertical_tail_root_chord - static_variable.vertical_tail_tip_chord 
    rudder_start_chord_wing_segment = static_variable.vertical_tail_root_chord - delta_chord * x[1]
    rudder_end_chord_wing_segment   = static_variable.vertical_tail_root_chord - delta_chord * 0.95
    rudder_start_chord              = static_variable.rudder_chord_fraction* rudder_start_chord_wing_segment 
    rudder_end_chord                = static_variable.rudder_chord_fraction* rudder_end_chord_wing_segment 
    Area_rudder                     = rudder_span  * ( rudder_start_chord +  rudder_end_chord) / 2
      
    Total_Area = Area_aileron + Area_rudder
    print('Objective: ', Total_Area)
    return Total_Area

# hard constraint
def con_Y_cw(x,static_variable):     
    
    C_w             = static_variable.C_w
    phi_cw          = static_variable.phi_cw
    C_Y_0_cw        = static_variable.C_Y_0_cw
    C_Y_beta        = static_variable.C_Y_beta  
    beta_cw         = static_variable.beta_cw 
    C_Y_delta_a_fz  = static_variable.C_Y_delta_a_fz
    C_Y_delta_r_fz  = static_variable.C_Y_delta_r_fz
    
    res = C_w*np.sin(phi_cw) + C_Y_0_cw + C_Y_beta*beta_cw + C_Y_delta_a_fz(x[0])*x[2] + C_Y_delta_r_fz(x[1])*x[3]
    
    return res

def con_L_cw(x,static_variable): 

    C_L_0_cw        = static_variable.C_L_0_cw
    C_L_beta        = static_variable.C_L_beta 
    beta_cw         = static_variable.beta_cw 
    C_L_delta_a_fz  = static_variable.C_L_delta_a_fz
    C_L_delta_r_fz  = static_variable.C_L_delta_r_fz
    
    res = C_L_0_cw + C_L_beta*beta_cw + C_L_delta_a_fz(x[0])*x[2] +  C_L_delta_r_fz(x[1])*x[3]
    
    return res

def con_N_cw(x,static_variable):  
    
    C_N_0_cw        = static_variable.C_N_0_cw
    C_N_beta        = static_variable.C_N_beta
    beta_cw         = static_variable.beta_cw 
    C_N_delta_a_fz  = static_variable.C_N_delta_a_fz
    C_N_delta_r_fz  = static_variable.C_N_delta_r_fz
    
    res = C_N_0_cw + C_N_beta*beta_cw + C_N_delta_a_fz(x[0])*x[2] +  C_N_delta_r_fz(x[1])*x[3]
    
    return res

def con_Y_oei(x,static_variable):   
    
    C_w              = static_variable.C_w
    phi_oei          = static_variable.phi_oei
    C_Y_0_oei        = static_variable.C_Y_0_oei
    C_Y_beta         = static_variable.C_Y_beta 
    beta_oei         = static_variable.beta_oei     
    C_Y_delta_a_fz   = static_variable.C_Y_delta_a_fz
    C_Y_delta_r_fz   = static_variable.C_Y_delta_r_fz
    
    res = C_w*np.sin(phi_oei) + C_Y_0_oei + C_Y_beta*beta_oei + C_Y_delta_a_fz(x[5])*x[2] + C_Y_delta_r_fz(x[6])*x[3]
    
    return res

def con_L_oei(x,static_variable): 
    
    C_L_0_oei        = static_variable.C_L_0_oei
    C_L_beta         = static_variable.C_L_beta 
    beta_oei         = static_variable.beta_oei 
    C_L_delta_a_fz   = static_variable.C_L_delta_a_fz
    C_L_delta_r_fz   = static_variable.C_L_delta_r_fz  
    
    res = C_L_0_oei + C_L_beta*beta_oei + C_L_delta_a_fz(x[5])*x[2] +  C_L_delta_r_fz(x[6])*x[3]
    
    return res

def con_N_oei(x,static_variable): 
    
    C_N_0_oei        = static_variable.C_N_0_oei
    C_N_beta         = static_variable.C_N_beta
    beta_oei         = static_variable.beta_oei 
    C_N_delta_a_fz   = static_variable.C_N_delta_a_fz
    C_N_delta_r_fz   = static_variable.C_N_delta_r_fz 
    
    res              = C_N_0_oei + C_N_beta*beta_oei + C_N_delta_a_fz( x[5])*x[2] +  C_N_delta_r_fz(x[6])*x[3]
    
    return res

def con_roll_rate(x,static_variable):     
    
    rho             = static_variable.rho
    V_app           = 1.3*static_variable.V_stall
    b               = static_variable.wing_span
    CS_angle        = static_variable.CS_angle
    CS_time         = static_variable.CS_time 
    CD_R            = static_variable.CD_R  
    C_L_delta_a_fz  = static_variable.C_L_delta_a_fz
    Ixx             = static_variable.Ixx
    S_w             = static_variable.S_w
    S_h             = static_variable.S_h
    S_v             = static_variable.S_v
    
    y_D             = 0.4*b/2
    C_L_Rolling     = C_L_delta_a_fz(x[5])*x[4]
    L_A             = 0.5*rho*(V_app**2)*S_w*b*C_L_Rolling
    P_ss            = np.sqrt(2*np.abs(L_A)/(rho*(S_w + S_h + S_v)*CD_R*(y_D)**3))
    Phi_1           = (Ixx/(rho*((y_D)**3)*(S_w + S_h + S_v)*CD_R))*np.log(P_ss**2)
    P_dot           = (P_ss**2)/(2*np.abs(Phi_1)) 
    
    res             = CS_time - np.sqrt(2 * CS_angle / P_dot)
    
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