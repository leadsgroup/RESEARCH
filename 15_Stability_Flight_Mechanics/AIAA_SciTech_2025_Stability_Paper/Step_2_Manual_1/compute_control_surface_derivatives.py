
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from   RCAIDE.Framework.Core import Units,Data  
from   RCAIDE.Library.Plots  import *  
import numpy                 as np
from   scipy.optimize        import fsolve
from   sklearn.linear_model  import LinearRegression
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
    
    derivatives_aileron1                      = np.array([0.5, 0.6])
    derivatives_aileron2                      = np.array([0.5, 0.6])
    derivatives_rudder1                       = np.array([0.5, 0.6])     
    derivatives_rudder2                       = np.array([0.5, 0.6]) 
    
    # Create the functions describing the control derivatives as functions of control surface size
    
    fa1                                       = regression(aileron_size, derivatives_aileron1)
    fa2                                       = regression(aileron_size, derivatives_aileron2)
    fr1                                       = regression(rudder_size, derivatives_rudder1)    
    fr2                                       = regression(rudder_size, derivatives_rudder2)
    
    # Import data from vehicle
    
    W                                         = 2700*9.81
    rho                                       = 1.225
    v                                         = 110.  * Units['mph']
    S                                         = 15.629           
    b                                         = 11.82855
    T                                         = 525              # [N]
    arm                                       = 9.2              # [m]
    delta_a                                   = 30 * np.pi/180   # [rad/s] 
    delta_r                                   = 30 * np.pi/180   # [rad/s]
    beta                                      = 5 * np.pi/180    # [rad/s]
    C_l_0                                     = 0
    C_l_beta                                  = 1
    C_n_beta                                  = 1                                                          
    C_n_0                                     = -(T*arm)/(0.5*rho*(v**2)*S*b)
    
    #delta_r                                   = (-(C_l_0/C_l_delta_r) - (C_l_beta/C_l_delta_r)*beta + (C_l_delta_a/C_l_delta_r)*(C_n_0/C_n_delta_a) + (C_l_delta_a/C_l_delta_r)*(C_n_beta/C_n_delta_a)*beta)/(1 - ((C_l_delta_a)/(C_l_delta_r))*((C_n_delta_r)/(C_n_delta_a)))
    #delta_a                                   = (-C_n_0 - C_n_beta*beta - C_n_delta_r*delta_r)/C_n_delta_a
    
    # Setup the system of equations to solve
    
    def system(vars):
        x, y = vars  
        eq1 = C_l_beta*beta + fa1(x)*delta_a + fr1(y)*delta_r
        eq2 = C_n_0 + C_n_beta*beta + fa2(x)*delta_a + fr2(y)*delta_r
        return [eq1, eq2]
    
    # Solve the system

    initial_guess              = [1.0, 2.0]
    solution                   = fsolve(system, initial_guess)
    delta_a_size, delta_r_size = solution    
    
    return delta_a_size, delta_r_size 

def regression(x,y):
    
    # ------------------------------------------------------------------
    #   Regression 
    # ------------------------------------------------------------------ 
    
    x_reshaped      = x.reshape(-1, 1)
    linear_regressor  = LinearRegression()
    linear_regressor.fit(x_reshaped, y)
    slope             = linear_regressor.coef_[0]
    intercept         = linear_regressor.intercept_
    
    def linear_regression_function1(x):
        return slope * x + intercept
    
    x_new            = np.linspace(x.min(), x.max(), 100)
    y_new            = linear_regression_function1(x_new)   
    dy_dx            = np.full_like(x_new, slope)    
    
    return linear_regression_function1

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
 
def mission_setup(analyses): 
    
    '''
    This sets up the nominal cruise of the aircraft
    '''
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
    segment.true_course                                              = 30 * Units.degree
                                                                     
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

from scipy.optimize import minimize 

# ----------------------------------------------------------------------------------------------------------------------
#  design motor 
# ----------------------------------------------------------------------------------------------------------------------     
def design_motor(motor):
    ''' Sizes a DC motor to obtain the best combination of speed constant and resistance values
    by sizing the motor for a design RPM value. Note that this design RPM value can be compute
    from design tip mach. The following properties are computed.  
      motor.speed_constant  (float): speed-constant [untiless] 
      motor.resistance      (float): resistance     [ohms]        
    
    Assumptions:
        None 
    
    Source:
        None
    
    Args:
        motor.no_load_current  (float): no-load current  [A]
        motor.nominal_voltage  (float): nominal voltage  [V]
        motor.angular_velocity (float): angular velocity [radians/s]    
        motor.efficiency       (float): efficiency       [unitless]
        motor.design_torque    (float): design torque    [Nm]
       
    Returns:
       None 
    
    '''     
    # design properties of the motor 
    io    = motor.no_load_current
    v     = motor.nominal_voltage
    omeg  = motor.angular_velocity     
    etam  = motor.efficiency 
    Q     = motor.design_torque 
    
    # define optimizer bounds 
    KV_lower_bound  = 0.01
    Res_lower_bound = 0.001
    KV_upper_bound  = 100
    Res_upper_bound = 10
    
    args       = (v , omeg,  etam , Q , io ) 
    hard_cons  = [{'type':'eq', 'fun': hard_constraint_1,'args': args},{'type':'eq', 'fun': hard_constraint_2,'args': args}] 
    slack_cons = [{'type':'eq', 'fun': slack_constraint_1,'args': args},{'type':'eq', 'fun': slack_constraint_2,'args': args}]  
    bnds       = ((KV_lower_bound, KV_upper_bound), (Res_lower_bound , Res_upper_bound)) 
    
    # try hard constraints to find optimum motor parameters
    sol = minimize(objective, [0.5, 0.1], args=(v , omeg,  etam , Q , io) , method='SLSQP', bounds=bnds, tol=1e-6, constraints=hard_cons) 
    
    if sol.success == False:
        # use slack constraints if optimizer fails and motor parameters cannot be found 
        print('\n Optimum motor design failed. Using slack constraints')
        sol = minimize(objective, [0.5, 0.1], args=(v , omeg,  etam , Q , io) , method='SLSQP', bounds=bnds, tol=1e-6, constraints=slack_cons) 
        if sol.success == False:
            assert('\n Slack contraints failed')  
    
    motor.speed_constant   = sol.x[0]
    motor.resistance       = sol.x[1]    
    
    return   
  
# objective function
def objective(x, v , omeg,  etam , Q , io ): 
    return (v - omeg/x[0])/x[1]   

# hard efficiency constraint
def hard_constraint_1(x, v , omeg,  etam , Q , io ): 
    return etam - (1- (io*x[1])/(v - omeg/x[0]))*(omeg/(v*x[0]))   

# hard torque equality constraint
def hard_constraint_2(x, v , omeg,  etam , Q , io ): 
    return ((v - omeg/x[0])/x[1] - io)/x[0] - Q  

# slack efficiency constraint 
def slack_constraint_1(x, v , omeg,  etam , Q , io ): 
    return abs(etam - (1- (io*x[1])/(v - omeg/x[0]))*(omeg/(v*x[0]))) - 0.2

# slack torque equality constraint 
def slack_constraint_2(x, v , omeg,  etam , Q , io ): 
    return  abs(((v - omeg/x[0])/x[1] - io)/x[0] - Q) - 200 

 
if __name__ == '__main__': 
    main()     