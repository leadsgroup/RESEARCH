# Procedure.py 
# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------     
import RCAIDE
from RCAIDE.Framework.Core import Units, Data
import numpy as np
from RCAIDE.Framework.Analyses.Process import Process    
from RCAIDE.Framework.Mission.Common   import Results, State  
from RCAIDE.Framework.Analyses.Process                     import Process    
from RCAIDE.Framework.Mission.Common                       import Results, State 
from RCAIDE.Library.Methods.Stability.Common               import compute_dynamic_flight_modes   
from RCAIDE.Library.Methods.Weights.Moment_of_Inertia      import compute_aircraft_moment_of_inertia
from RCAIDE.Library.Methods.Weights.Center_of_Gravity      import compute_vehicle_center_of_gravity  
from RCAIDE.Library.Methods.Aerodynamics.Athena_Vortex_Lattice.run_AVL_analysis  import run_AVL_analysis 

# Routines  
import Missions 
from copy import  deepcopy

# ----------------------------------------------------------------------        
#   Setup
# ----------------------------------------------------------------------    

def stick_fixed_stability_and_drag_procedure(): 
    procedure                 = Process()
    procedure.modify_vehicle  = modify_stick_fixed_vehicle   
    procedure.post_process    = longitudinal_static_stability_and_drag_post_process   
        
    return procedure 

def elevator_sizing_setup(): 
    procedure = Process()  
    procedure.post_process  = elevator_sizing_post_process   
    return procedure   

def aileron_rudder_sizing_setup(): 
    procedure = Process()  
    procedure.post_process  = aileron_rudder_sizing_post_process   
    return procedure   

def flap_sizing_setup(): 
    procedure = Process()    
    procedure.post_process  = flap_sizing_post_process 
    return procedure  

# ----------------------------------------------------------------------      
#   Modify Vehicle 
# ----------------------------------------------------------------------  

def modify_stick_fixed_vehicle(nexus): 
    '''
    This function takes the updated design variables and modifies the aircraft 
    '''
    # Pull out the vehicles
    vehicle = nexus.vehicle_configurations.stick_fixed_cruise   
        
    # Update Wing    
    for wing in vehicle.wings:
        if type(wing) ==  RCAIDE.Library.Components.Wings.Main_Wing: 
            vehicle.reference_area        = wing.areas.reference            
        Sref                     = wing.areas.reference      # fixed 
        span                     = wing.spans.projected      # optimization input
        taper                    = wing.taper                # optimization input
        croot                    = 2*Sref/((taper+1)*span)   # set by Sref and current design point
        ctip                     = taper * croot             # set by Sref and current design point 
        wing.chords.root         = croot
        wing.chords.tip          = ctip 
        
        # Wing Segments
        if 'Segments' in wing:
            for seg in wing.Segments:
                seg.twist = (wing.twists.tip-wing.twists.root)*seg.percent_span_location  + wing.twists.root
                
        # update remaning wing properties  
        wing_segmented_planform(wing)  
        wing.areas.wetted             = wing.areas.reference  * 2 
        wing.areas.exposed            = wing.areas.reference  * 2  
                
      
    # Update MOI 
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weight_analysis.vehicle                       = vehicle 
    results                                       = weight_analysis.evaluate()

    compute_vehicle_center_of_gravity(weight_analysis.vehicle)     
    CG_location      = vehicle.mass_properties.center_of_gravity 
    compute_aircraft_moment_of_inertia(weight_analysis.vehicle, CG_location)
    
    # moment of intertia matrix cannot be negative so set to 0
    vehicle.mass_properties.moments_of_inertia.tensor[0, 2] = 0
    vehicle.mass_properties.moments_of_inertia.tensor[2, 0] = 0

    # Update Mission  
    nexus.missions = Missions.stick_fixed_stability_setup(nexus.analyses,weight_analysis.vehicle,nexus.cruise_velocity, nexus.cruise_altitude)      
    
    # diff the new data
    vehicle.store_diff() 
    
    return nexus   

 
def longitudinal_static_stability_and_drag_post_process(nexus): 
    '''
    This function analyses and post processes the aircraft at cruise conditions. 
    The objective of is to minimize the drag  of a trimmed aircraft 
    '''
    summary                                                 = nexus.summary 
    vehicle                                                 = nexus.vehicle_configurations.stick_fixed_cruise 
    g                                                       = 9.81   
    L                                                       = g*vehicle.mass_properties.max_takeoff
    S                                                       = vehicle.reference_area
    atmosphere                                              = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data                                               = atmosphere.compute_values(altitude = nexus.missions['stick_fixed_cruise'].segments['cruise'].altitude )       
                                     

    # ------------------------------------------------------------------------------------------------------------------------  
    # Stick Fixed 
    # ------------------------------------------------------------------------------------------------------------------------       
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions                                   
    stick_fixed_conditions                                          = Results()
    stick_fixed_conditions.expand_rows(2)
    stick_fixed_conditions.freestream.density[:,0]                  = atmo_data.density[0,0]
    stick_fixed_conditions.freestream.gravity[:,0]                  = g          
    stick_fixed_conditions.freestream.speed_of_sound[:,0]           = atmo_data.speed_of_sound[0,0] 
    stick_fixed_conditions.freestream.velocity[:,0]                 = nexus.missions['stick_fixed_cruise'].segments['cruise'].air_speed 
    stick_fixed_conditions.frames.inertial.velocity_vector[:,0]     = nexus.missions['stick_fixed_cruise'].segments['cruise'].air_speed 
    stick_fixed_conditions.freestream.mach_number                   = stick_fixed_conditions.freestream.velocity/stick_fixed_conditions.freestream.speed_of_sound
    stick_fixed_conditions.freestream.dynamic_pressure              = 0.5 * stick_fixed_conditions.freestream.density *  (stick_fixed_conditions.freestream.velocity ** 2)
    stick_fixed_conditions.aerodynamics.angles.beta[:,0]            = 0.0 
    stick_fixed_conditions.aerodynamics.angles.alpha[:,0]           = 0.0 
    stick_fixed_conditions.aerodynamics.coefficients.lift.total     = np.array([[L, L*1.01]]).T /(S*(0.5*stick_fixed_conditions.freestream.density*(stick_fixed_conditions.freestream.velocity**2))) 
    stick_fixed_conditions.static_stability.coefficients.roll[:,0]  =  0.0
    stick_fixed_conditions.static_stability.coefficients.pitch[:,0] =  0.0

    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses                                          
    stability_stick_fixed                                   = RCAIDE.Framework.Analyses.Stability.Athena_Vortex_Lattice() 
    stability_stick_fixed.settings.filenames.avl_bin_name   = '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' # eg. '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' 
    stability_stick_fixed.settings.trim_aircraft            = True
    stability_stick_fixed.vehicle                           = nexus.vehicle_configurations.stick_fixed_cruise

    # ------------------------------------------------------------------------------------------------------------------------  
    # Run AVL      
    run_AVL_analysis(stability_stick_fixed,stick_fixed_conditions) 
    stick_fixed_state = State()
    stick_fixed_state.conditions = stick_fixed_conditions     
    compute_dynamic_flight_modes(stick_fixed_state,stability_stick_fixed.settings,vehicle)

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------      
    summary.CD              = stick_fixed_conditions.aerodynamics.coefficients.drag.induced.total[0,0]  # OBJECTIVE FUNCTION MULTIPLEID BY 10 
    summary.CM_residual     = abs(stick_fixed_conditions.static_stability.coefficients.pitch[0,0])
    summary.spiral_criteria = stick_fixed_conditions.static_stability.spiral_criteria[0,0]
    NP                      = stick_fixed_conditions.static_stability.neutral_point[0,0]
    cg                      = vehicle.mass_properties.center_of_gravity[0][0]
    MAC                     = vehicle.wings.main_wing.chords.mean_aerodynamic
    summary.static_margin   = (NP - cg)/MAC
    summary.CM_alpha        = stick_fixed_conditions.static_stability.derivatives.CM_alpha[0,0]  
 
    if np.count_nonzero(vehicle.mass_properties.moments_of_inertia.tensor) > 0:  
        summary.phugoid_damping_ratio       = stick_fixed_conditions.dynamic_stability.LongModes.phugoidDamping[0,0] 
        summary.short_period_damping_ratio  = stick_fixed_conditions.dynamic_stability.LongModes.shortPeriodDamping[0,0] 
        summary.dutch_roll_frequency        = stick_fixed_conditions.dynamic_stability.LatModes.dutchRollFreqHz[0,0]* (2 * np.pi)  # converting to omega
        summary.dutch_roll_damping_ratio    = stick_fixed_conditions.dynamic_stability.LatModes.dutchRollDamping[0,0]
        summary.spiral_doubling_time        = stick_fixed_conditions.dynamic_stability.LatModes.spiralTimeDoubleHalf[0,0] 
        print("Drag Coefficient           : " + str(summary.CD))
        print("Moment Coefficient         : " + str(summary.CM_residual))
        print("Static Margin              : " + str(summary.static_margin))
        print("CM alpla                   : " + str(summary.CM_alpha))   
        print("Phugoid Damping Ratio      : " + str(summary.phugoid_damping_ratio))
        print("Short Period Damping Ratio : " + str(summary.short_period_damping_ratio))
        print("Dutch Roll Frequency       : " + str(summary.dutch_roll_frequency))
        print("Dutch Roll Damping Ratio   : " + str(summary.dutch_roll_damping_ratio))
        print("Spiral Doubling Time       : " + str(summary.spiral_doubling_time)) 
        print("Spiral Criteria            : " + str(summary.spiral_criteria))
        print("\n\n") 

    else: 
        summary.phugoid_damping_ratio    = 0 
        summary.dutch_roll_frequency     = 0
        summary.dutch_roll_damping_ratio = 0
        summary.spiral_doubling_time     = 0
        summary.spiral_criteria          = 0 
        print("Drag Coefficient         : " + str(summary.CD))
        print("Moment Coefficient       : " + str(summary.CM_residual))
        print("Static Margin            : " + str(summary.static_margin))
        print("CM alpla                 : " + str(summary.CM_alpha))    
        print("Spiral Criteria          : " + str(summary.spiral_criteria))
        print("\n\n")
        
    return nexus  
    
def elevator_sizing_post_process(nexus): 
    '''
    This function analyses and post processes the aircraft at the flight conditions required to size
    the elevator. These conditions are:
    1) Stick pull maneuver with a load factor of 3.0
    2) Stick push maneuver with a load factor of -1
    ''' 
    summary                                            = nexus.summary  
    g                                                  = 9.81 
    vehicle                                            = nexus.vehicle_configurations.elevator_sizing 
    m                                                  = vehicle.mass_properties.max_takeoff
    S                                                  = vehicle.reference_area 
    V_trim                                             = vehicle.trim_airspeed   
    V_max                                              = nexus.missions['elevator_sizing'].segments['cruise'].air_speed 
    atmosphere                                         = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data                                          = atmosphere.compute_values(altitude = nexus.missions['elevator_sizing'].segments['cruise'].altitude )
     
    q            = 0.5*(V_trim**2)*atmo_data.density[0,0] 
    CL_pull_man  = vehicle.maxiumum_load_factor*m*g/(S*q)  
    CL_push_man  = vehicle.minimum_load_factor*m*g/(S*q)
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # Pull Up Maneuver 
    # ------------------------------------------------------------------------------------------------------------------------       
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions            
    pull_up_conditions                                               = Results()
    pull_up_conditions.freestream.density                            = np.array([[atmo_data.density[0,0]]])
    pull_up_conditions.freestream.gravity                            = np.array([[g]])           
    pull_up_conditions.freestream.speed_of_sound                     = np.array([[atmo_data.speed_of_sound[0,0]]]) 
    pull_up_conditions.freestream.dynamic_pressure                   = 0.5 * pull_up_conditions.freestream.density *  (pull_up_conditions.freestream.velocity ** 2)    
    pull_up_conditions.aerodynamics.angles.beta                      = np.array([[0.0]])
    pull_up_conditions.aerodynamics.angles.alpha                     = np.array([[0.0]])
    pull_up_conditions.static_stability.coefficients.roll            = np.array([[0.0]])
    pull_up_conditions.static_stability.coefficients.pitch           = np.array([[0.0]])  
    pull_up_conditions.aerodynamics.coefficients.lift.total          = CL_pull_man
    pull_up_conditions.freestream.velocity                           = np.array([[V_max]])
    pull_up_conditions.frames.inertial.velocity_vector[:,0]          = V_max
    pull_up_conditions.freestream.mach_number                        = pull_up_conditions.freestream.velocity/pull_up_conditions.freestream.speed_of_sound
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses          
    stability_pull_maneuver                                      = RCAIDE.Framework.Analyses.Stability.Athena_Vortex_Lattice() 
    stability_pull_maneuver.settings.filenames.avl_bin_name      = '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' # eg. '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35'
    stability_pull_maneuver.settings.trim_aircraft               = True    
    stability_pull_maneuver.settings.number_of_spanwise_vortices = 40
    stability_pull_maneuver.vehicle                              = vehicle

    # ------------------------------------------------------------------------------------------------------------------------  
    # Run AVL   
    run_AVL_analysis(stability_pull_maneuver,pull_up_conditions)

    # ------------------------------------------------------------------------------------------------------------------------  
    # Push Over Maneuver 
    # ------------------------------------------------------------------------------------------------------------------------          
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions    
    push_over_conditions                                               = Results()
    push_over_conditions.freestream.density                            = np.array([[atmo_data.density[0,0]]]) 
    push_over_conditions.freestream.gravity                            = np.array([[g]])           
    push_over_conditions.freestream.speed_of_sound                     = np.array([[atmo_data.speed_of_sound[0,0]]]) 
    push_over_conditions.freestream.dynamic_pressure                   = 0.5 * push_over_conditions.freestream.density *  (push_over_conditions.freestream.velocity ** 2) 
    push_over_conditions.aerodynamics.angles.beta                      = np.array([[0.0]])
    push_over_conditions.aerodynamics.angles.alpha                     = np.array([[0.0]])
    push_over_conditions.static_stability.coefficients.roll            = np.array([[0.0]])
    push_over_conditions.static_stability.coefficients.pitch           = np.array([[0.0 ]])  
    push_over_conditions.aerodynamics.coefficients.lift.total          = CL_push_man
    push_over_conditions.freestream.velocity                           = V_trim
    push_over_conditions.frames.inertial.velocity_vector               = V_trim
    push_over_conditions.freestream.mach_number                        = push_over_conditions.freestream.velocity/push_over_conditions.freestream.speed_of_sound

    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses                                    
    stability_push_maneuver                                      = RCAIDE.Framework.Analyses.Stability.Athena_Vortex_Lattice() 
    stability_push_maneuver.settings.filenames.avl_bin_name      =  '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' # eg. '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35'
    stability_push_maneuver.settings.trim_aircraft               = True 
    stability_push_maneuver.settings.number_of_spanwise_vortices = 40
    stability_push_maneuver.vehicle                              = vehicle
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run AVL      
    run_AVL_analysis(stability_push_maneuver,push_over_conditions)   

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------ 
    AoA_pull                                                     = pull_up_conditions.aerodynamics.angles.alpha[0,0]
    elevator_pull_deflection                                     = pull_up_conditions.static_stability.control_surfaces_cases['case_0001_0001'].control_surfaces.elevator.deflection 
    AoA_push                                                     = push_over_conditions.aerodynamics.angles.alpha[0,0]
    elevator_push_deflection                                     = push_over_conditions.static_stability.control_surfaces_cases['case_0001_0001'].control_surfaces.elevator.deflection
    
    summary.elevator_pull_deflection  = abs(elevator_pull_deflection)
    summary.elevator_push_deflection  = abs(elevator_push_deflection) 
    
    # compute control surface area 
    control_surfaces = ['elevator'] 
    total_control_surface_area = compute_control_surface_areas(control_surfaces,vehicle)  
    summary.elevator_surface_area =  total_control_surface_area
    

    print("Elevator Area      : " + str(summary.elevator_surface_area))
    print("Aircraft CL Pull   : " + str(CL_pull_man[0, 0]))
    print("Aircraft AoA Pull  : " + str(AoA_pull))
    print("Elevator Pull Defl.: " + str(elevator_pull_deflection)) 
    print("Aircraft CL Push   : " + str(CL_push_man[0, 0]))
    print("Aircraft AoA Push  : " + str(AoA_push))
    print("Elevator Push Defl.: " + str(elevator_push_deflection)) 
    print("\n\n")     
         
    return nexus    

 

def aileron_rudder_sizing_post_process(nexus):  
    '''
    This function analyses and post processes the aircraft at the flight conditions required to size
    the aileron and rudder. These conditions are:
    1) A controlled roll at a  rate of 0.07
    2) Trimmed flight in a 20 knot crosswind
    ''' 
    summary                                                      = nexus.summary  
    g                                                            = 9.81 
    vehicle                                                      = nexus.vehicle_configurations.aileron_rudder_sizing 
    CL_trim                                                      = vehicle.trim_cl 
    V_crosswind                                                  = vehicle.crosswind_velocity
                                          
    atmosphere                                                   = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data                                                    = atmosphere.compute_values(altitude = nexus.missions['aileron_sizing'].segments['cruise'].altitude )
    
 
    # ------------------------------------------------------------------------------------------------------------------------  
    # Roll Maneuver 
    # ------------------------------------------------------------------------------------------------------------------------       
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions      
    roll_conditions                                               = Results()
    roll_conditions.freestream.density                            = np.array([[atmo_data.density[0,0]]])  
    roll_conditions.freestream.gravity                            = np.array([[g ]])           
    roll_conditions.freestream.speed_of_sound                     = np.array([[atmo_data.speed_of_sound[0,0]]])  
    roll_conditions.freestream.dynamic_pressure                   = 0.5 * roll_conditions.freestream.density *  (roll_conditions.freestream.velocity ** 2)    
    roll_conditions.aerodynamics.angles.beta                      = np.array([[0.0]]) 
    roll_conditions.aerodynamics.angles.alpha                     = np.array([[0.0]]) 
    roll_conditions.aerodynamics.coefficients.lift.total          = CL_trim 
    roll_conditions.freestream.velocity                           = nexus.missions['aileron_sizing'].segments['cruise'].air_speed 
    roll_conditions.frames.inertial.velocity_vector[:, 0]         = nexus.missions['aileron_sizing'].segments['cruise'].air_speed 
    roll_conditions.freestream.mach_number                        = roll_conditions.freestream.velocity/roll_conditions.freestream.speed_of_sound
    roll_conditions.static_stability.coefficients.roll            = np.array([[ 0.07]]) 
    roll_conditions.static_stability.coefficients.pitch           = np.array([[0.0]]) 
    roll_conditions.aerodynamics.angles.beta                      = np.array([[0.0]])  

    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses       
    stability_roll_maneuver                                      = RCAIDE.Framework.Analyses.Stability.Athena_Vortex_Lattice() 
    stability_roll_maneuver.settings.filenames.avl_bin_name      = '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' # eg. '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35'
    stability_roll_maneuver.settings.number_of_spanwise_vortices = 40 
    stability_roll_maneuver.settings.trim_aircraft               = True
    stability_roll_maneuver.vehicle                              = vehicle

    # ------------------------------------------------------------------------------------------------------------------------  
    # Run AVL   
    run_AVL_analysis(stability_roll_maneuver,roll_conditions)   
    aileron_roll_deflection                                      = roll_conditions.static_stability.control_surfaces_cases['case_0001_0001'].control_surfaces.aileron.deflection 

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------      
    summary.aileron_roll_deflection =  abs(aileron_roll_deflection) 
    if vehicle.rudder_flag: 
        rudder_roll_deflection  = roll_conditions.static_stability.control_surfaces_cases['case_0001_0001'].control_surfaces.rudder.deflection
        summary.rudder_roll_deflection = abs(rudder_roll_deflection) 
    else:
        rudder_roll_deflection = 0
        summary.rudder_roll_deflection = 0       
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # Crosswind Maneuver 
    # ------------------------------------------------------------------------------------------------------------------------          
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions            
    crosswind_conditions                                                = Results()
    crosswind_conditions.freestream.density                             = np.array([[atmo_data.density[0,0]]])   
    crosswind_conditions.freestream.gravity                             = np.array([[g ]])         
    crosswind_conditions.freestream.speed_of_sound                      = np.array([[atmo_data.speed_of_sound[0,0]]])    
    crosswind_conditions.freestream.velocity                            = np.array([[nexus.missions['aileron_sizing'].segments['cruise'].air_speed]])
    crosswind_conditions.freestream.dynamic_pressure                    = 0.5 * crosswind_conditions.freestream.density *  (crosswind_conditions.freestream.velocity ** 2) 
    crosswind_conditions.freestream.mach_number                         = crosswind_conditions.freestream.velocity/crosswind_conditions.freestream.speed_of_sound   
    crosswind_conditions.aerodynamics.angles.beta                       = np.array([[0.0]])  
    crosswind_conditions.aerodynamics.angles.alpha                      = np.array([[0.0]])     
    crosswind_conditions.static_stability.coefficients.roll             = np.array([[0.0]])  
    crosswind_conditions.static_stability.coefficients.pitch            = np.array([[0.0]]) 
    crosswind_conditions.aerodynamics.coefficients.lift.total           = CL_trim 
    crosswind_conditions.aerodynamics.angles.beta                       = np.array([[np.tan(V_crosswind/nexus.missions['aileron_sizing'].segments['cruise'].air_speed) ]])# beta
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses      
    stability_cross_wind_maneuver = RCAIDE.Framework.Analyses.Stability.Athena_Vortex_Lattice() 
    stability_cross_wind_maneuver.settings.filenames.avl_bin_name ='/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' # eg. '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35'
    stability_cross_wind_maneuver.settings.trim_aircraft          =  True
    stability_cross_wind_maneuver.vehicle                         = vehicle
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run AVL       
    run_AVL_analysis(stability_cross_wind_maneuver,crosswind_conditions)
    aileron_cross_wind_deflection                                 = crosswind_conditions.static_stability.control_surfaces_cases['case_0001_0001'].control_surfaces.aileron.deflection 

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------      
    # criteria 
    summary.aileron_crosswind_deflection = abs(aileron_cross_wind_deflection)  

    if vehicle.rudder_flag: 
        rudder_cross_wind_deflection  = crosswind_conditions.static_stability.control_surfaces_cases['case_0001_0001'].control_surfaces.rudder.deflection
        summary.rudder_crosswind_deflection =  abs(rudder_cross_wind_deflection)  
    else:
        rudder_cross_wind_deflection = 0
        summary.rudder_crosswind_deflection = 0  
        
    # compute control surface area 
    control_surfaces = ['aileron','rudder'] 
    total_control_surface_area = compute_control_surface_areas(control_surfaces,vehicle)   
    summary.aileron_rudder_surface_area =  total_control_surface_area

    print("Total Rudder Aileron Surface Area : " + str(summary.aileron_rudder_surface_area)) 
    print("Aileron Roll Defl                 : " + str(aileron_roll_deflection)) 
    print("Rudder Roll Defl                  : " + str(rudder_roll_deflection))  
    print("Aileron Crosswind Defl            : " + str(aileron_cross_wind_deflection)) 
    print("Rudder  Crosswind Defl            : " + str(rudder_cross_wind_deflection )) 
    print("\n\n")     
  
    return nexus     


def flap_sizing_post_process(nexus): 
    '''
    This function analyses and post processes the aircraft at the flight conditions required to size
    the flap. These conditions are:
    1) A comparison of clean and deployed flap at 12 deg. angle of attack
    ''' 
    summary                                                  = nexus.summary  
    g                                                        = 9.81 
    vehicle                                                  = nexus.vehicle_configurations.flap_sizing  
    V_max                                                    = nexus.missions['flap_sizing'].segments['cruise'].air_speed 
    atmosphere                                               = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data                                                = atmosphere.compute_values(altitude = nexus.missions['flap_sizing'].segments['cruise'].altitude )
    

    # ------------------------------------------------------------------------------------------------------------------------  
    # Flaps Not Deployed 
    # ------------------------------------------------------------------------------------------------------------------------          
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions            
    no_flap_conditions                                           = Results()
    no_flap_conditions.freestream.density                        = np.array([[atmo_data.density[0,0] ]]) 
    no_flap_conditions.freestream.gravity                        = np.array([[g ]])           
    no_flap_conditions.freestream.speed_of_sound                 = np.array([[atmo_data.speed_of_sound[0,0]  ]]) 
    no_flap_conditions.freestream.velocity                       = np.array([[V_max]])  
    no_flap_conditions.freestream.mach_number                    = no_flap_conditions.freestream.velocity/no_flap_conditions.freestream.speed_of_sound
    no_flap_conditions.aerodynamics.angles.beta                  = np.array([[0.0]]) 
    no_flap_conditions.aerodynamics.coefficients.lift.total      = None
    no_flap_conditions.aerodynamics.angles.alpha                 = np.array([[12.0]])*Units.degrees
    no_flap_conditions.static_stability.coefficients.roll        = np.array([[0.0]]) 
    no_flap_conditions.static_stability.coefficients.pitch       = np.array([[0.0]]) 

    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions                  
    stability_no_flap = RCAIDE.Framework.Analyses.Stability.Athena_Vortex_Lattice() 
    stability_no_flap.settings.filenames.avl_bin_name        =  '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' # eg. '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' 
    stability_no_flap.settings.number_of_spanwise_vortices   = 40
    stability_no_flap.settings.trim_aircraft                 = False
    vehicle.wings.main_wing.control_surfaces.flap.deflection = 0.0
    stability_no_flap.vehicle                                = vehicle
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run AVL       
    run_AVL_analysis(stability_no_flap,no_flap_conditions)
    CL_12_deg_no_flap                                        = no_flap_conditions.aerodynamics.coefficients.lift.total[0,0]  
      

    # ------------------------------------------------------------------------------------------------------------------------  
    # Flaps  Deployed 
    # ------------------------------------------------------------------------------------------------------------------------   
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions            
    w_flap_conditions = deepcopy(no_flap_conditions)
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions                     
    stability_flap = RCAIDE.Framework.Analyses.Stability.Athena_Vortex_Lattice() 
    stability_flap.settings.filenames.avl_bin_name           = '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' # eg. '/Users/matthewclarke/Documents/LEADS/CODES/AVL/avl3.35' 
    stability_flap.settings.number_of_spanwise_vortices      = 40
    stability_flap.settings.trim_aircraft                    = False
    vehicle.wings.main_wing.control_surfaces.flap.deflection = 40*Units.degrees 
    stability_flap.vehicle                                   = vehicle

    # ------------------------------------------------------------------------------------------------------------------------  
    # Run AVL           
    run_AVL_analysis(stability_flap,w_flap_conditions)
    CL_12_deg_flap  = w_flap_conditions.aerodynamics.coefficients.lift.total[0,0]     


    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------      
    # critera     
    flap_criteria  = (CL_12_deg_flap-CL_12_deg_no_flap) - 0.95*(CL_12_deg_flap-CL_12_deg_no_flap) 
    
    # compute control surface area 
    control_surfaces = ['flap'] 
    total_control_surface_area = compute_control_surface_areas(control_surfaces,vehicle)  
    summary.flap_surface_area  = total_control_surface_area
    summary.flap_criteria      = flap_criteria

    print("Flap Area     : " + str(summary.flap_surface_area))
    print("Flap Criteria : " + str(flap_criteria))  # https://aviation.stackexchange.com/questions/48715/how-is-the-area-of-flaps-determined
    print("\n\n")     
         
    return nexus    


def  compute_control_surface_areas(control_surfaces,vehicle): 
    '''
    This function computes the control suface area used in the objectives of the 
    control surface sizing scripts 
    '''
    total_control_surface_area = 0 
    for cs_idx in range(len(control_surfaces)):
        for wing in vehicle.wings:
            if getattr(wing,'control_surfaces',False):  
                for CS in wing.control_surfaces:
                    if CS.tag == control_surfaces[cs_idx]:
                        if wing.Segments: 
                            num_segs = len(wing.Segments)
                            segment_names = list(wing.Segments.keys())   
                            for seg_id in range(num_segs-1): 
                                current_seg =  segment_names[seg_id]
                                next_seg    =  segment_names[seg_id+1] 
                                if (CS.span_fraction_start >= wing.Segments[current_seg].percent_span_location) \
                                   and (CS.span_fraction_end <= wing.Segments[next_seg].percent_span_location): 
                                    root_chord             = wing.Segments[current_seg].root_chord_percent*wing.chords.root
                                    tip_chord              = wing.Segments[next_seg].root_chord_percent*wing.chords.root
                                    span                   = (wing.Segments[next_seg].percent_span_location-wing.Segments[current_seg].percent_span_location)*wing.spans.projected
                                    rel_start_percent_span = CS.span_fraction_start - wing.Segments[current_seg].percent_span_location
                                    rel_end_percent_span   = CS.span_fraction_end - wing.Segments[current_seg].percent_span_location
                                    chord_fraction         = CS.chord_fraction 
                                    area = conpute_control_surface_area(root_chord,tip_chord,span,rel_start_percent_span,rel_end_percent_span,chord_fraction)
                                    total_control_surface_area += area
                        else: 
                            root_chord             = wing.chords.root
                            tip_chord              = wing.chords.tip
                            span                   = wing.spans.projected
                            rel_start_percent_span = CS.span_fraction_start  
                            rel_end_percent_span   = CS.span_fraction_end 
                            chord_fraction         = CS.chord_fraction 
                            area = conpute_control_surface_area(root_chord,tip_chord,span,rel_start_percent_span,rel_end_percent_span,chord_fraction)                            
                            total_control_surface_area += area 
         
    return total_control_surface_area

def conpute_control_surface_area(root_chord,tip_chord,span,rel_start_percent_span,rel_end_percent_span,chord_fraction):  
    '''
    This is a simple function that computes the area of a single control surface
    '''
    cs_start_chord =  (root_chord +   ((tip_chord-root_chord)/span)*(rel_start_percent_span*span))*chord_fraction
    cs_end_chord   =  (root_chord +   ((tip_chord-root_chord)/span)*(rel_end_percent_span*span))*chord_fraction
    cs_span        = (rel_end_percent_span-rel_start_percent_span)*span
    cs_area        = 0.5*(cs_start_chord+cs_end_chord)*cs_span
    return cs_area
    

def linear_discretize(x_pivs,chi,pivot_points):
    
    chi_pivs       = np.zeros(len(x_pivs))
    chi_pivs[0]    = chi[0]
    chi_pivs[-1]   = chi[-1]
    locations      = np.array(pivot_points)*(chi[-1]-chi[0]) + chi[0]
    
    # vectorize 
    chi_2d = np.repeat(np.atleast_2d(chi).T,len(pivot_points),axis = 1)
    pp_2d  = np.repeat(np.atleast_2d(locations),len(chi),axis = 0)
    idxs   = (np.abs(chi_2d - pp_2d)).argmin(axis = 0) 
    
    chi_pivs[1:-1] = chi[idxs]
    
    x_bar  = np.interp(chi,chi_pivs, x_pivs) 
    
    return x_bar 


def spline_discretize(x_pivs,chi,pivot_points):
    chi_pivs       = np.zeros(len(x_pivs))
    chi_pivs[0]    = chi[0]
    chi_pivs[-1]   = chi[-1]
    locations      = np.array(pivot_points)*(chi[-1]-chi[0]) + chi[0]
    
    # vectorize 
    chi_2d = np.repeat(np.atleast_2d(chi).T,len(pivot_points),axis = 1)
    pp_2d  = np.repeat(np.atleast_2d(locations),len(chi),axis = 0)
    idxs   = (np.abs(chi_2d - pp_2d)).argmin(axis = 0) 
    
    chi_pivs[1:-1] = chi[idxs]
    
    fspline = interpolate.CubicSpline(chi_pivs,x_pivs)
    x_bar = fspline(chi)
    
    return x_bar