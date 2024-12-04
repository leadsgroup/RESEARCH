
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

def main():
    
    # -------------------------------------------------------------
    #  Aileron & Rudder
    # -------------------------------------------------------------
    
    # Create arrays of aileron and rudder sizes
    
    aileron_size                              = np.array([0.5, 0.6])
    rudder_size                               = np.array([0.5, 0.6])
    
    # Load optimized vehicle
    
    #Optimized_Vehicle_1_pkl                   = "025_00_00_40_00_00_Optimized_Vehicle"
    #vehicle                                   = load_results(Optimized_Vehicle_1_pkl)
    
    # Append aileron and rudder
    
    #vehicle                                   = setup_rudder_aileron(vehicle)
    
    # Compute the control derivatives related to the control surfaces
    
    #derivatives, results                      = compute_rudder_aileron_derivatives(aileron_size, rudder_size, vehicle, seg_num=0)
    
    C_Y_delta_a_cw                            = np.array([0.5, 0.6])
    C_L_delta_a_cw                            = np.array([0.5, 0.6])
    C_N_delta_a_cw                            = np.array([0.5, 0.6])     
    C_Y_delta_a_oei                           = np.array([0.5, 0.6])
    C_L_delta_a_oei                           = np.array([0.5, 0.6])
    C_N_delta_a_oei                           = np.array([0.5, 0.6])
    C_Y_delta_r_cw                            = np.array([0.5, 0.6])
    C_L_delta_r_cw                            = np.array([0.5, 0.6])
    C_N_delta_r_cw                            = np.array([0.5, 0.6])     
    C_Y_delta_r_oei                           = np.array([0.5, 0.6])
    C_L_delta_r_oei                           = np.array([0.5, 0.6])
    C_N_delta_r_oei                           = np.array([0.5, 0.6])    
    
    # Create the functions describing the control derivatives as functions of control surface size
    
    static_variable = Data()
    
    static_variable.C_Y_delta_a_cw            = interpolate.interp1d(aileron_size, C_Y_delta_a_cw)
    static_variable.C_L_delta_a_cw            = interpolate.interp1d(aileron_size, C_L_delta_a_cw)
    static_variable.C_N_delta_a_cw            = interpolate.interp1d(aileron_size, C_N_delta_a_cw)
    static_variable.C_Y_delta_a_oei           = interpolate.interp1d(aileron_size, C_Y_delta_a_oei)
    static_variable.C_L_delta_a_oei           = interpolate.interp1d(aileron_size, C_L_delta_a_oei)
    static_variable.C_N_delta_a_oei           = interpolate.interp1d(aileron_size, C_N_delta_a_oei)    
    static_variable.C_Y_delta_r_cw            = interpolate.interp1d(rudder_size,  C_Y_delta_r_cw)
    static_variable.C_L_delta_r_cw            = interpolate.interp1d(rudder_size,  C_L_delta_r_cw)
    static_variable.C_N_delta_r_cw            = interpolate.interp1d(rudder_size,  C_N_delta_r_cw)
    static_variable.C_Y_delta_r_oei           = interpolate.interp1d(rudder_size,  C_Y_delta_r_oei)
    static_variable.C_L_delta_r_oei           = interpolate.interp1d(rudder_size,  C_L_delta_r_oei)
    static_variable.C_N_delta_r_oei           = interpolate.interp1d(rudder_size,  C_N_delta_r_oei)
    
    # Import data from vehicle
    
    # CROSSWIND 
    
    static_variable.W_cw                                          = 2700*9.81
    static_variable.rho_cw                                        = 1.225
    static_variable.S_cw                                          = 15.629           
    static_variable.wing_span_cw                                  = 11.82855 
    static_variable.vertical_tail_span_cw                         = 2
    static_variable.aileron_chord_cw                              = 0.1
    static_variable.rudder_chord_cw                               = 0.1
    static_variable.T_cw                                          = 525                           # [N]
    static_variable.arm_cw                                        = 9.2                           # [m]
    static_variable.delta_a_cw                                    = 30 * np.pi/180                # [rad/s] 
    static_variable.delta_r_cw                                    = 30 * np.pi/180                # [rad/s]
    static_variable.beta_cw                                       = 5 * np.pi/180                 # [rad/s]
    static_variable.C_l_0_cw                                      = 0
    static_variable.C_l_beta_cw                                   = 1
    static_variable.C_n_beta_cw                                   = 1                                                          
    static_variable.C_n_0_cw                                      = 0
    static_variable.V_cw                                          = 110. * Units['mph']
    static_variable.v_cw                                          = 35                            # [ft/s] 
    static_variable.beta_cw                                       = np.arcsin(static_variable.v_cw/static_variable.V_cw) # [rad]
    static_variable.C_w_cw                                        = static_variable.W_cw/(0.5*static_variable.rho_cw*(static_variable.V_cw**2)*static_variable.S_cw)                                                                      #[-]
    
    # OEI
    
    static_variable.W_oei                                         = 2700*9.81
    static_variable.rho_oei                                       = 1.225
    static_variable.V_oei                                         = 110.  * Units['mph']
    static_variable.S_oei                                         = 15.629           
    static_variable.b_oei                                         = 11.82855
    static_variable.aileron_chord_oei                             = 0.1
    static_variable.rudder_chord_oei                              = 0.1
    static_variable.T_oei                                         = 525              # [N]
    static_variable.arm_oei                                       = 9.2              # [m]
    static_variable.delta_a_oei                                   = 30 * np.pi/180   # [rad/s] 
    static_variable.delta_r_oei                                   = 30 * np.pi/180   # [rad/s]
    static_variable.beta_oei                                      = 0 * np.pi/180    # [rad/s]
    static_variable.C_l_0_oei                                     = 0
    static_variable.C_l_beta_oei                                  = 1
    static_variable.C_n_beta_oei                                  = 1                                                          
    static_variable.C_n_0_oei                                     = -(static_variable.T_oei*static_variable.arm_oei)/(0.5*static_variable.rho_oei*(static_variable.V_oei**2)*static_variable.S_oei*static_variable.b_oei)
    static_variable.C_w_oei                                       = static_variable.W_oei/(0.5*static_variable.rho_oei*(static_variable.V_oei**2)*static_variable.S_oei)
    
    delta_a_lower_bound_cw, delta_r_lower_bound_cw, delta_a_lower_bound_oei, delta_r_lower_bound_oei, aileron_span_lower_bound, rudder_span_lower_bound = optimization(static_variable)
    
    debug = 0
    
    # Setup the system of equations to solve
    
    #def system(vars):
        #x, y = vars  
        #eq1 = C_l_beta*beta + fa1(x)*delta_a + fr1(y)*delta_r
        #eq2 = C_n_0 + C_n_beta*beta + fa2(x)*delta_a + fr2(y)*delta_r
        #return [eq1, eq2]
    
    ## Solve the system

    #initial_guess              = [1.0, 2.0]
    #solution                   = fsolve(system, initial_guess)
    #delta_a_size, delta_r_size = solution    
    
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
    aileron_span_lower_bound = 0.75               # [%]    # just changing the upper and lower bound of the lower bound of the control surface (the limit closer to root of wing)
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
    
    args       = static_variable
    
    hard_cons  = [{'type':'eq', 'fun': hard_constraint_Y_cw,'args': args},
                  {'type':'eq', 'fun': hard_constraint_L_cw,'args': args},
                  {'type':'eq', 'fun': hard_constraint_N_cw,'args': args},
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
    sol = minimize(objective, [0.5, 0.1], args= (static_variable) , method='SLSQP', bounds=bnds, tol=1e-6, constraints=hard_cons) 
    
    if sol.success == False:
        print('\n Optimum motor design failed.')
    
    delta_a_lower_bound_cw   = sol.x[0]
    delta_r_lower_bound_cw   = sol.x[1]    
    delta_a_lower_bound_oei  = sol.x[2] 
    delta_r_lower_bound_oei  = sol.x[3] 
    aileron_span_lower_bound = sol.x[4] 
    rudder_span_lower_bound  = sol.x[5] 
    
    return delta_a_lower_bound_cw, delta_r_lower_bound_cw, delta_a_lower_bound_oei, delta_r_lower_bound_oei, aileron_span_lower_bound, rudder_span_lower_bound 

# objective function
def objective(design_variables,static_variable):
    aileron_span_cw    = static_variable.wing_span_cw*(0.95 - design_variables[4])
    rudder_span_cw     = static_variable.vertical_tail_span_cw*(0.95 - design_variables[5])    
    return (aileron_span_cw*static_variable.aileron_chord_cw)*(rudder_span_cw*static_variable.rudder_chord_cw)

# hard constraint
def hard_constraint_Y_cw(design_variables,static_variable): 
    
    aileron_span    = static_variable.wing_span_cw*(0.95 - design_variables[4])
    rudder_span     = static_variable.vertical_tail_span_cw*(0.95 - design_variables[5])
    
    C_w_cw          = static_variable.C_w_cw
    phi_cw          = static_variable.phi_cw
    C_Y_0_cw        = static_variable.C_Y_0_cw
    C_Y_beta_cw     = static_variable.C_Y_beta_cw
    beta_cw         = static_variable.beta_cw 
    C_Y_delta_a_cw  = static_variable.C_Y_beta_cw
    C_Y_delta_r_cw  = static_variable.C_Y_beta_cw
    
    res = C_w_cw*np.sin(phi_cw) + C_Y_0_cw + C_Y_beta_cw*beta_cw + C_Y_delta_a_cw(aileron_span)*design_variables[0] + C_Y_delta_r_cw(rudder_span)*design_variables[1]
    
    return res

def hard_constraint_L_cw(design_variables,static_variable): 
    
    aileron_span    = static_variable.wing_span_cw*(0.95 - design_variables[4])
    rudder_span     = static_variable.vertical_tail_span_cw*(0.95 - design_variables[5])
    
    C_L_0_cw        = static_variable.C_L_0_cw
    C_L_beta_cw     = static_variable.C_L_beta_cw
    beta_cw         = static_variable.beta_cw 
    C_L_delta_a_cw  = static_variable.C_L_beta_cw
    C_L_delta_r_cw  = static_variable.C_L_beta_cw
    
    res = C_L_0_cw + C_L_beta_cw*beta_cw + C_L_delta_a_cw(aileron_span)*design_variables[0] +  C_L_delta_r_cw(rudder_span)*design_variables[1]
    
    return res

def hard_constraint_N_cw(design_variables,static_variable): 
    
    aileron_span    = static_variable.wing_span_cw*(0.95 - design_variables[4])
    rudder_span     = static_variable.vertical_tail_span_cw*(0.95 - design_variables[5])
    
    C_N_0_cw        = static_variable.C_N_0_cw
    C_N_beta_cw     = static_variable.C_N_beta_cw
    beta_cw         = static_variable.beta_cw 
    C_N_delta_a_cw  = static_variable.C_N_beta_cw
    C_N_delta_r_cw  = static_variable.C_N_beta_cw
    
    res = C_N_0_cw + C_N_beta_cw*beta_cw + C_N_delta_a_cw(aileron_span)*design_variables[0] +  C_N_delta_r_cw(rudder_span)*design_variables[1]
    
    return res

def hard_constraint_Y_oei(design_variables,static_variable): 
    
    aileron_span     = static_variable.wing_span_oei*(0.95 - design_variables[4])
    rudder_span      = static_variable.vertical_tail_span_oei*(0.95 - design_variables[5])
    
    C_w_oei          = static_variable.C_w_oei
    phi_oei          = static_variable.phi_oei
    C_Y_0_oei        = static_variable.C_Y_0_oei
    C_Y_beta_oei     = static_variable.C_Y_beta_oei
    beta_oei         = static_variable.beta_oei 
    C_Y_delta_a_oei  = static_variable.C_Y_beta_oei
    C_Y_delta_r_oei  = static_variable.C_Y_beta_oei
    
    res = C_w_oei*np.sin(phi_oei) + C_Y_0_oei + C_Y_beta_oei*beta_oei + C_Y_delta_a_oei(aileron_span)*design_variables[2] + C_Y_delta_r_oei(rudder_span)*design_variables[3]
    
    return res

def hard_constraint_L_oei(design_variables,static_variable): 
    
    aileron_span     = static_variable.wing_span_oei*(0.95 - design_variables[4])
    rudder_span      = static_variable.vertical_tail_span_oei*(0.95 - design_variables[5])
    
    C_L_0_oei        = static_variable.C_L_0_oei
    C_L_beta_oei     = static_variable.C_L_beta_oei
    beta_oei         = static_variable.beta_oei 
    C_L_delta_a_oei  = static_variable.C_L_beta_oei
    C_L_delta_r_oei  = static_variable.C_L_beta_oei
    
    res = C_L_0_oei + C_L_beta_oei*beta_oei + C_L_delta_a_oei(aileron_span)*design_variables[2] +  C_L_delta_r_oei(rudder_span)*design_variables[3]
    
    return res

def hard_constraint_N_oei(design_variables,static_variable): 
    
    aileron_span     = static_variable.wing_span_oei*(0.95 - design_variables[4])
    rudder_span      = static_variable.vertical_tail_span_oei*(0.95 - design_variables[5])
    
    C_N_0_oei        = static_variable.C_N_0_oei
    C_N_beta_oei     = static_variable.C_N_beta_oei
    beta_oei         = static_variable.beta_oei 
    C_N_delta_a_oei  = static_variable.C_N_beta_oei
    C_N_delta_r_oei  = static_variable.C_N_beta_oei
    
    res = C_N_0_oei + C_N_beta_oei*beta_oei + C_N_delta_a_oei(aileron_span)*design_variables[2] +  C_N_delta_r_oei(rudder_span)*design_variables[3]
    
    return res
 
if __name__ == '__main__': 
    main()     