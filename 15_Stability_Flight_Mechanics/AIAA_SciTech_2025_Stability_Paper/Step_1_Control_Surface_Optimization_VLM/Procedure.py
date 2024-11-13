# Procedure.py 
# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------     
import RCAIDE
from RCAIDE.Framework.Core import Units, Data
from RCAIDE.Library.Methods.Geometry.Planform              import wing_segmented_planform   
from RCAIDE.Framework.Analyses.Process                     import Process    
from RCAIDE.Framework.Mission.Common                       import Results, State 
from RCAIDE.Library.Methods.Stability.Common               import compute_dynamic_flight_modes   
from RCAIDE.Library.Methods.Weights.Moment_of_Inertia      import compute_aircraft_moment_of_inertia
from RCAIDE.Library.Methods.Weights.Center_of_Gravity      import compute_vehicle_center_of_gravity
from RCAIDE.Library.Methods.Aerodynamics.Vortex_Lattice_Method.evaluate_VLM import evaluate_no_surrogate
from RCAIDE.Library.Mission.Common.Update  import orientations
from RCAIDE.Library.Mission.Common.Unpack_Unknowns import orientation

# Routines  
import Missions 

import numpy as np
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
    summary     = nexus.summary 
    vehicle     = nexus.vehicle_configurations.stick_fixed_cruise 
    g           = 9.81   
    m           = vehicle.mass_properties.max_takeoff 
    S           = vehicle.reference_area
    atmosphere  = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data   = atmosphere.compute_values(altitude = nexus.missions['stick_fixed_cruise'].segments['cruise'].altitude )       
      
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # Stick Fixed 
    # ------------------------------------------------------------------------------------------------------------------------       
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions                                         
    stick_fixed_conditions                                             = Results()
    stick_fixed_conditions.freestream.density[:,0]                     = atmo_data.density[0,0]
    stick_fixed_conditions.freestream.gravity[:,0]                     = g          
    stick_fixed_conditions.freestream.speed_of_sound[:,0]              = atmo_data.speed_of_sound[0,0] 
    stick_fixed_conditions.freestream.dynamic_viscosity[:,0]           = atmo_data.dynamic_viscosity[0,0] 
    stick_fixed_conditions.freestream.temperature[:,0]                 = atmo_data.temperature[0,0] 
    stick_fixed_conditions.freestream.pressure[:,0]                    = atmo_data.pressure[0,0]   
    stick_fixed_conditions.freestream.altitude[:,0]                    = nexus.missions['stick_fixed_cruise'].segments['cruise'].altitude     
    stick_fixed_conditions.freestream.velocity[:,0]                    = nexus.missions['stick_fixed_cruise'].segments['cruise'].air_speed 
    stick_fixed_conditions.frames.inertial.velocity_vector[:,0]        = nexus.missions['stick_fixed_cruise'].segments['cruise'].air_speed 
    stick_fixed_conditions.freestream.mach_number                      = stick_fixed_conditions.freestream.velocity/stick_fixed_conditions.freestream.speed_of_sound
    stick_fixed_conditions.freestream.dynamic_pressure                 = 0.5 * stick_fixed_conditions.freestream.density *  (stick_fixed_conditions.freestream.velocity ** 2)
    stick_fixed_conditions.freestream.reynolds_number                  = stick_fixed_conditions.freestream.density * stick_fixed_conditions.freestream.velocity / stick_fixed_conditions.freestream.dynamic_viscosity  
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses                                      
    stability_stick_fixed                                       = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method()   
    stability_stick_fixed.settings.number_of_spanwise_vortices  = 20
    stability_stick_fixed.settings.number_of_chordwise_vortices = 4
    stability_stick_fixed.settings.discretize_control_surfaces  = False
    stability_stick_fixed.vehicle                               = vehicle 
    stability_stick_fixed.settings.use_surrogate = False 
    stability_stick_fixed.initialize()
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run VLM
    stick_fixed_segment =  nexus.missions.stick_fixed_cruise.segments['cruise']
    stick_fixed_segment.conditions = stick_fixed_conditions
    stick_fixed_segment.state.conditions = stick_fixed_conditions
    orientation(stick_fixed_segment)
    orientations(stick_fixed_segment)
    evaluate_no_surrogate(stick_fixed_segment,stability_stick_fixed.settings,vehicle) 
    compute_dynamic_flight_modes(stick_fixed_segment,stability_stick_fixed.settings,vehicle)

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------     
    summary.CD              = stick_fixed_conditions.static_stability.coefficients.drag[0,0]  
    summary.CM_residual     = abs(stick_fixed_conditions.static_stability.coefficients.M[0,0])
    summary.spiral_criteria = stick_fixed_conditions.static_stability.spiral_criteria[0,0]
    NP                      = stick_fixed_conditions.static_stability.neutral_point[0,0]
    cg                      = vehicle.mass_properties.center_of_gravity[0][0]
    MAC                     = vehicle.wings.main_wing.chords.mean_aerodynamic
    summary.static_margin   = (NP - cg)/MAC
    summary.CM_alpha        = stick_fixed_conditions.static_stability.derivatives.CM_alpha[0,0] 

    CL_stick_fixed                   = stick_fixed_conditions.aerodynamics.coefficients.lift.total
    CL_stick_fixed_required          = m*g/(S*stick_fixed_conditions.freestream.dynamic_pressure)    
    summary.CL_stick_fixed_residual  = abs(CL_stick_fixed[0, 0] - CL_stick_fixed_required[0, 0])    
 
    if np.count_nonzero(vehicle.mass_properties.moments_of_inertia.tensor) > 0:  
        summary.phugoid_damping_ratio       = stick_fixed_conditions.dynamic_stability.LongModes.phugoidDamping[0,0] 
        summary.short_period_damping_ratio  = stick_fixed_conditions.dynamic_stability.LongModes.shortPeriodDamping[0,0] 
        summary.dutch_roll_frequency        = stick_fixed_conditions.dynamic_stability.LatModes.dutchRollFreqHz[0,0]* (2 * np.pi)  
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
    summary     = nexus.summary  
    g           = 9.81 
    vehicle     = nexus.vehicle_configurations.elevator_sizing 
    m           = vehicle.mass_properties.max_takeoff
    S           = vehicle.reference_area 
    V_trim      = nexus.cruise_velocity   
    V_max       = nexus.missions['elevator_sizing'].segments['cruise'].air_speed 
    atmosphere  = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data   = atmosphere.compute_values(altitude = nexus.missions['elevator_sizing'].segments['cruise'].altitude )

    
    # ------------------------------------------------------------------------------------------------------------------------  
    # Pull Up Maneuver 
    # ------------------------------------------------------------------------------------------------------------------------       
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions        
    pull_up_conditions                                       = Results()
    pull_up_conditions.freestream.density                    = np.array([[atmo_data.density[0,0]]])
    pull_up_conditions.freestream.gravity                    = np.array([[g]])           
    pull_up_conditions.freestream.speed_of_sound             = np.array([[atmo_data.speed_of_sound[0,0]]]) 
    pull_up_conditions.freestream.dynamic_pressure           = 0.5 * pull_up_conditions.freestream.density *  (pull_up_conditions.freestream.velocity ** 2)     
    pull_up_conditions.freestream.velocity                   = np.array([[V_max]])
    pull_up_conditions.frames.inertial.velocity_vector[:,0]  = V_max
    pull_up_conditions.freestream.mach_number                = pull_up_conditions.freestream.velocity/pull_up_conditions.freestream.speed_of_sound 
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses                                      
    stability_pull_up                                       = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method()   
    stability_pull_up.settings.number_of_spanwise_vortices  = 40
    stability_pull_up.settings.number_of_chordwise_vortices = 10
    stability_pull_up.settings.discretize_control_surfaces  = False
    stability_pull_up.vehicle                               = vehicle 
    stability_pull_up.settings.use_surrogate = False 
    stability_pull_up.initialize()
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run VLM
    pull_up_segment =  nexus.missions.elevator_sizing.segments['cruise']
    pull_up_segment.conditions = pull_up_conditions
    pull_up_segment.state.conditions = pull_up_conditions
    orientation(pull_up_segment)
    orientations(pull_up_segment)
    evaluate_no_surrogate(pull_up_segment,stability_pull_up.settings,vehicle) 
    compute_dynamic_flight_modes(pull_up_segment,stability_pull_up.settings,vehicle) 
    
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
    stability_push_over                                       = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method()   
    stability_push_over.settings.number_of_spanwise_vortices  = 40
    stability_push_over.settings.number_of_chordwise_vortices = 10
    stability_push_over.settings.discretize_control_surfaces  = False
    stability_push_over.vehicle                               = vehicle 
    stability_push_over.settings.use_surrogate = False 
    stability_push_over.initialize()

    # ------------------------------------------------------------------------------------------------------------------------  
    # Run VLM
    push_over_segment =  nexus.missions.elevator_sizing.segments['cruise']
    push_over_segment.conditions = push_over_conditions
    push_over_segment.state.conditions = push_over_conditions
    orientation(push_over_segment)
    orientations(push_over_segment)
    evaluate_no_surrogate(push_over_segment,stability_push_over.settings,vehicle) 
    compute_dynamic_flight_modes(push_over_segment,stability_push_over.settings,vehicle)  

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------  
    AoA_pull                          = pull_up_conditions.aerodynamics.angles.alpha[0,0]
    AoA_push                          = push_over_conditions.aerodynamics.angles.alpha[0,0]    
    
    # compute control surface area 
    control_surfaces = ['elevator'] 
    total_control_surface_area = compute_control_surface_areas(control_surfaces,vehicle)  
    summary.elevator_surface_area =  total_control_surface_area
    
    # compute lift coefficient and residual from target CLs 
    CL_pull_man                  = pull_up_conditions.aerodynamics.coefficients.lift.total
    CL_push_man                  = push_over_conditions.aerodynamics.coefficients.lift.total
    CL_pull_man_required         = vehicle.maxiumum_load_factor*m*g/(S*pull_up_conditions.freestream.dynamic_pressure) 
    CL_push_man_required         = vehicle.minimum_load_factor*m*g/(S*push_over_conditions.freestream.dynamic_pressure)
    summary.CL_pull_residual     = abs(CL_pull_man - CL_pull_man_required)  
    summary.CL_push_residual     = abs(CL_push_man - CL_push_man_required)  
    
    # compute elevator deflections 
    elevator_pull_deflection          = pull_up_conditions.control_surfaces.elevator.deflection  
    elevator_push_deflection          = push_over_conditions.control_surfaces.elevator.deflection 
    summary.elevator_pull_deflection  = abs(elevator_pull_deflection)  
    summary.elevator_push_deflection  = abs(elevator_push_deflection)  

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
    summary      = nexus.summary  
    g            = 9.81   
    m            = vehicle.mass_properties.max_takeoff
    S            = vehicle.reference_area 
    vehicle      = nexus.vehicle_configurations.aileron_rudder_sizing 
    #CL_trim      = vehicle.trim_cl 
    V_crosswind  = vehicle.crosswind_velocity 
    atmosphere   = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmo_data    = atmosphere.compute_values(altitude = nexus.missions['aileron_sizing'].segments['cruise'].altitude )
    

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
    roll_conditions.freestream.velocity                           = nexus.missions['aileron_sizing'].segments['cruise'].air_speed 
    roll_conditions.frames.inertial.velocity_vector[:, 0]         = nexus.missions['aileron_sizing'].segments['cruise'].air_speed 
    roll_conditions.freestream.mach_number                        = roll_conditions.freestream.velocity/roll_conditions.freestream.speed_of_sound
    roll_conditions.static_stability.coefficients.roll            = np.array([[ 0.07]]) 
    roll_conditions.static_stability.coefficients.pitch           = np.array([[0.0]]) 
    roll_conditions.aerodynamics.angles.beta                      = np.array([[0.0]])
    
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses                                      
    stability_roll                                       = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method()   
    stability_roll.settings.number_of_spanwise_vortices  = 40
    stability_roll.settings.number_of_chordwise_vortices = 10
    stability_roll.settings.discretize_control_surfaces  = False
    stability_roll.vehicle                               = vehicle 
    stability_roll.settings.use_surrogate = False 
    stability_roll.initialize()
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run VLM
    roll_segment =  nexus.missions.turn_criteria.segments['cruise']
    roll_segment.conditions = roll_conditions
    roll_segment.state.conditions = roll_conditions
    orientation(roll_segment)
    orientations(roll_segment)
    evaluate_no_surrogate(roll_segment,stability_roll.settings,vehicle) 
    compute_dynamic_flight_modes(roll_segment,stability_roll.settings,vehicle)
    
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
    crosswind_conditions.static_stability.coefficients.roll             = np.array([[0.0]])  
    crosswind_conditions.static_stability.coefficients.pitch            = np.array([[0.0]])
    crosswind_conditions.aerodynamics.angles.beta                       = np.array([[np.tan(V_crosswind/nexus.missions['aileron_sizing'].segments['cruise'].air_speed) ]])  
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses                                      
    crosswind_maneuver                                       = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method()   
    crosswind_maneuver.settings.number_of_spanwise_vortices  = 40
    crosswind_maneuver.settings.number_of_chordwise_vortices = 10
    crosswind_maneuver.settings.discretize_control_surfaces  = False
    crosswind_maneuver.vehicle                               = vehicle 
    crosswind_maneuver.settings.use_surrogate = False 
    crosswind_maneuver.initialize()
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run VLM
    crosswind_segment =  nexus.missions.turn_criteria.segments['cruise']
    crosswind_segment.conditions = roll_conditions
    crosswind_segment.state.conditions = roll_conditions
    orientation(crosswind_segment)
    orientations(crosswind_segment)
    evaluate_no_surrogate(crosswind_segment,crosswind_maneuver.settings,vehicle) 
    compute_dynamic_flight_modes(crosswind_segment,crosswind_maneuver.settings,vehicle)
    

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------  

    # ------------------------------------------------------------------------------------------------------------------------  
    # Aileron 
    # compute lift coefficient and residual from target CLs 
    CL_roll_man                   = roll_conditions.aerodynamics.coefficients.lift.total
    CL_crosswind_man              = crosswind_conditions.aerodynamics.coefficients.lift.total
    CL_roll_required              = vehicle.maxiumum_load_factor*m*g/(S*roll_conditions.freestream.dynamic_pressure) 
    CL_crosswind_required         = vehicle.minimum_load_factor*m*g/(S*crosswind_conditions.freestream.dynamic_pressure)
    summary.CL_roll_residual      = abs(CL_roll_man - CL_roll_required)  
    summary.CL_crosswind_residual = abs(CL_crosswind_man - CL_crosswind_required)  
    
    # compute aileron deflections 
    aileron_roll_deflection               = roll_conditions.control_surfaces.aileron.deflection  
    aileron_crosswind_deflection          = crosswind_conditions.control_surfaces.aileron.deflection 
    summary.aileron_roll_deflection       = abs(aileron_roll_deflection)  
    summary.aileron_crosswind_deflection  = abs(aileron_crosswind_deflection) 

    # ------------------------------------------------------------------------------------------------------------------------  
    # Rudder     
    if vehicle.rudder_flag: 
        rudder_roll_deflection               = roll_conditions.control_surfaces.rudder.deflection  
        rudder_crosswind_deflection          = crosswind_conditions.control_surfaces.rudder.deflection 
        summary.rudder_roll_deflection       = abs(rudder_roll_deflection)  
        summary.rudder_crosswind_deflection  = abs(rudder_crosswind_deflection)  
    else:
        rudder_roll_deflection              = 0 
        summary.rudder_roll_deflection      = 0  
        rudder_crosswind_deflection         = 0  
        summary.rudder_crosswind_deflection = 0 
        
    # compute control surface area 
    control_surfaces = ['aileron','rudder'] 
    total_control_surface_area = compute_control_surface_areas(control_surfaces,vehicle)   
    summary.aileron_rudder_surface_area =  total_control_surface_area

    print("Total Rudder Aileron Surface Area : " + str(summary.aileron_rudder_surface_area)) 
    print("Aileron Roll Defl                 : " + str(aileron_roll_deflection)) 
    print("Rudder Roll Defl                  : " + str(rudder_roll_deflection))  
    print("Aileron Crosswind Defl            : " + str(aileron_crosswind_deflection)) 
    print("Rudder  Crosswind Defl            : " + str(rudder_crosswind_deflection )) 
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
    # initialize analyses                                      
    stability_no_flap                                       = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method()   
    stability_no_flap.settings.number_of_spanwise_vortices  = 40
    stability_no_flap.settings.number_of_chordwise_vortices = 10
    stability_no_flap.settings.discretize_control_surfaces  = False
    vehicle.wings.main_wing.control_surfaces.flap.deflection = 0.0 
    stability_no_flap.vehicle                               = vehicle 
    stability_no_flap.settings.use_surrogate = False 
    stability_no_flap.initialize()
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run VLM
    no_flap_segment =  nexus.missions.turn_criteria.segments['cruise']
    no_flap_segment.conditions = no_flap_conditions
    no_flap_segment.state.conditions = no_flap_conditions
    orientation(no_flap_segment)
    orientations(no_flap_segment)
    evaluate_no_surrogate(no_flap_segment,stability_no_flap.settings,vehicle) 
    compute_dynamic_flight_modes(no_flap_segment,stability_no_flap.settings,vehicle)

 
    # ------------------------------------------------------------------------------------------------------------------------  
    # Flaps  Deployed 
    # ------------------------------------------------------------------------------------------------------------------------          
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize run conditions
    w_flap_conditions = deepcopy(no_flap_conditions)
    w_flap_conditions.control_surfaces.aileron.deflection  = 40*Units.degrees    
      
    # ------------------------------------------------------------------------------------------------------------------------  
    # initialize analyses                                      
    stability_w_flap                                       = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method()   
    stability_w_flap.settings.number_of_spanwise_vortices  = 40
    stability_w_flap.settings.number_of_chordwise_vortices = 10
    stability_w_flap.settings.discretize_control_surfaces  = False
    vehicle.wings.main_wing.control_surfaces.flap.deflection = 40*Units.degrees
    stability_w_flap.vehicle                               = vehicle 
    stability_w_flap.settings.use_surrogate = False 
    stability_w_flap.initialize()
     
    # ------------------------------------------------------------------------------------------------------------------------  
    # Run VLM
    w_flap_segment =  nexus.missions.flap_sizing.segments['cruise']
    w_flap_segment.conditions = w_flap_conditions
    w_flap_segment.state.conditions = w_flap_conditions
    orientation(w_flap_segment)
    orientations(w_flap_segment)
    evaluate_no_surrogate(w_flap_segment,stability_w_flap.settings,vehicle) 
    compute_dynamic_flight_modes(w_flap_segment,stability_w_flap.settings,vehicle)
    

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------       
    # critera    
    CL_12_deg_no_flap    = no_flap_conditions.aerodynamics.coefficients.lift.total[0,0]  
    CL_12_deg_flap       = w_flap_conditions.aerodynamics.coefficients.lift.total[0,0]
    flap_criteria        = (CL_12_deg_flap-CL_12_deg_no_flap) - 0.95*(CL_12_deg_flap-CL_12_deg_no_flap) 
    
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