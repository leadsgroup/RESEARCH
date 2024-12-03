# --------------------------------------
# Imports
# --------------------------------------
import RCAIDE
from   RCAIDE.Framework.Core import Units
import os 
import pickle  
import numpy                 as np

def main():
    
    # ----------------------------------------------------------------------------
    # -1  load optimized_vehicle
    # ----------------------------------------------------------------------------    
    
    Optimized_Vehicle_1_pkl               = "025_00_00_40_00_00_Optimized_Vehicle"
    Optimized_Vehicle_1                   = load_results(Optimized_Vehicle_1_pkl)
    
    # ----------------------------------------------------------------------------
    # -2  create aircraft with any elevator 
    # ----------------------------------------------------------------------------     
    
    elevator                              = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                          = 'elevator'
    elevator.span_fraction_start          = 0.1
    elevator.span_fraction_end            = 0.9
    elevator.deflection                   = 0.0  * Units.deg
    elevator.chord_fraction               = 0.35
    Optimized_Vehicle_1.wings.horizontal_tail.control_surfaces.append(elevator)     
    
    # ----------------------------------------------------------------------------
    # -3  run 1st mission: Pull up 
    # ---------------------------------------------------------------------------- 
    
    results = run_mission(Optimized_Vehicle_1)
    
    C_m_alpha    = 1
    C_L_trim     = 1
    C_L_alpha    = 1
    C_m_0        = 1
    C_L_alpha    = 1
    C_m_delta_e  = 1
    C_L_delta_e  = 1
    C_m_alpha    = 1
    
    det          = C_L_alpha * C_m_delta_e - C_L_delta_e * C_m_alpha

    delta_e_trim = -(C_m_alpha*C_L_trim + C_L_alpha*C_m_0)/det
    
    delta_e_trim_MAX = - (C_m_0 * C_L_alpha)/det
    
    stability_criterion = ((4*C_L_alpha*W)/(det*(rho*S*v**3)))*(static_margin)

    C_L_trim        = (2*W*np.cos(gamma))/(rho*S*V**2)                                                      #[-]                                                                                           
    C_D             = C_D_0 + k*C_L_trim**2                                                                 #[-]    
    
    if stability_criterion <= 0:
        print("Error: Unstable")
    else:
        print("Stable") 

    
    #W               = 8000                                                                                  #[lbs]       
    #S               = 200                                                                                   #[ft**2]
    #MAC             = 6                                                                                     #[ft]
    #z_p             = -2                                                                                    #[ft]
    #S_t             = 40                                                                                    #[ft**2]
    #l_sign_t        = 6                                                                                     #[ft]
    #h_n_wb          = 0.25                                                                                  #[-]
    #a_t             = 4.5                                                                                   #[rad**-1]
    #a_e             = 2                                                                                     #[rad**-1]
    #i_t             = 0.047                                                                                 #[rad]    
    #C_L_alpha_wb    = 5.0                                                                                   #[rad**-1]
    #C_m_ac_wb       = -0.1                                                                                  #[-]
    #C_D_0           = 0.03                                                                                  #[-]
    #k               = 0.1                                                                                   #[-]
    #epsilon_0       = 0.01                                                                                  #[rad]
    #depsilon_dalpha = 0.3                                                                                   #[-]
    #V               = 421.952                                                                               #[ft/s]
    #rho             = 0.001756                                                                              #[slugs/ft**3]
    #gamma           = 0*np.pi/180                                                                           #[deg -> rad]
    #h               = 0.46                                                                                  #[-]
    #C_m_p_alpha     = 0                                                                                     #[rad**-1]  
    #V_H             = l_sign_t*S_t/(MAC*S)                                                                  #[-]
    #C_L_alpha       = C_L_alpha_wb + a_t*(S_t/S)*(1 - depsilon_dalpha)                                      #[rad**-1]
    #h_n_p           = h_n_wb + (a_t/C_L_alpha)*V_H*(1 - depsilon_dalpha) - (1/C_L_alpha)*C_m_p_alpha        #[-]
    #C_L_delta_e     = a_e*S_t/S                                                                             #[rad**-1]
    #C_L_trim        = (2*W*np.cos(gamma))/(rho*S*V**2)                                                      #[-]                                                                                           
    #C_D             = C_D_0 + k*C_L_trim**2                                                                 #[-]
    D               = 0.5*rho*(V**2)*S*C_D                                                                  #[lbs]
    T               = D + W*gamma                                                                           #[lbs]
    C_T             = T/(0.5*rho*(V**2)*S)                                                                  #[-]
    C_m_0_p         = C_T*z_p/MAC                                                                           #[-]
    C_m_0           = C_m_ac_wb + C_m_0_p + V_H*a_t*(epsilon_0 + i_t)*(1 + (1 -                             
                      depsilon_dalpha) * ((a_t*S_t)/(C_L_alpha*S)))                                         #[-]
    delta_e_trim    = (- C_m_0 - C_L_trim*(h - h_n_p))/(C_L_delta_e*(h_n_p - h_n_wb) - a_e*V_H)*180/np.pi   #[deg]
        
    print("The elevator deflection required to trim the aircraft in cruise configuration is:", "%0.3f" % delta_e_trim, "deg.") 
    
    # ---------------------------------------------------------------------------------------------------------------
    
    # b) Assume the aircraft is symmetric (but not in the trim condition). Compute the bank angle,
    #    aileron deflection, and rudder deflection necessary to trim the aircraft with 
    
    W                    = 8000                                                                              #[lbs]       
    b                    = 30                                                                                #[ft]        
    S                    = 200                                                                               #[ft**2]
    S_f                  = 40                                                                                #[ft**2]
    a_f                  = 4                                                                                 #[rad**-1]
    z_f                  = -3                                                                                #[ft]
    l_f                  = 18                                                                                #[ft]
    a_r                  = 1                                                                                 #[rad**-1]
    v_f_v_squared        = 0.85                                                                              #[-]
    dsigma_dbeta         = -0.1                                                                              #[-]
    C_Y_beta_wb          = -0.5                                                                              #[-]
    C_l_beta_wb          = -0.15                                                                             #[-]
    C_n_beta_wb          = -0.3                                                                              #[-]
    C_Y_delta_a          = 0                                                                                 #[-]
    C_l_delta_a          = -0.2                                                                              #[-]
    C_n_delta_a          = 0                                                                                 #[-]    
    C_n_beta             = C_n_beta_wb + ((S_f*l_f)/(S*b))*v_f_v_squared*a_f*(1 - dsigma_dbeta)     
    V                    = 200                                                                               #[ft/s]  
    v                    = 35                                                                                #[ft/s] 
    
    # crosswind. Assume standard sea-level conditions 
    
    rho                  = 0.002377                                                                          #[slug/ft**3] 
      
    beta                 = np.arcsin(v/V)                                                                    #[rad] 
    C_w                  = W/(0.5*rho*(V**2)*S)                                                              #[-] 
    C_l_beta             = a_f*(1 - dsigma_dbeta)*((S_f*z_f)/(S*b))*v_f_v_squared + C_l_beta_wb              #[-] 
    C_l_delta_r          = -a_r*((S_f*z_f)/(S*b))*v_f_v_squared                                              #[-] 
    C_n_delta_r          = -((S_f*l_f)/(S*b))*v_f_v_squared*a_r                                              #[-] 
    C_Y_beta             = -a_f*(1 - dsigma_dbeta)*(S_f/S)*v_f_v_squared + C_Y_beta_wb                       #[-] 
    C_Y_delta_r          = a_r*(S_f/S)*v_f_v_squared                                                         #[-] 
                         
    A                    = np.array ([[C_w, C_Y_delta_a, C_Y_delta_r],[0, C_l_delta_a, C_l_delta_r],[0, C_n_delta_a, C_n_delta_r]])
    B                    = np.array([-C_Y_beta*beta,-C_l_beta*beta,-C_n_beta*beta])
    C                    = np.linalg.solve(A,B)
                         
    phi                  = np.arcsin(C[0])*180/np.pi                                                         #[deg] 
    delta_a              = C[1]*180/np.pi                                                                    #[deg] 
    delta_r              = C[2]*180/np.pi                                                                    #[deg] 
    
    print("The bank angle, aileron deflection, and rudder deflection necessary to trim the aircraft are respectively: \nphi =",
          "%0.3f" % phi, "deg, delta_a =", "%0.3f" % delta_a, "deg, delta_r =", "%0.3f" % delta_r, "deg." )      
        
    '''
    RUN LOOP TO GET DERIVATIVE FUNCTIONS OF CONTROL SURFACES     
    '''
    vehicle = append_control_surfaces(vehicle)
     
    
    '''
    RUN LOOP TO GET DERIVATIVE FUNCTIONS OF CONTROL SURFACES     
    '''
    derivatives = compute_control_surface_derivatives(vehicle)        
    
    
    
    
    
    
    debug = 0
    
    return



def load_results(filename):  

    #current_dir = '/Users/aidanmolloy/Documents/LEADS/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations'
    current_dir = 'C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations'
    load_file = os.path.join(current_dir, filename + '.pkl')
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    
    return results

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



if __name__ == '__main__': 
    main()
    plt.show()