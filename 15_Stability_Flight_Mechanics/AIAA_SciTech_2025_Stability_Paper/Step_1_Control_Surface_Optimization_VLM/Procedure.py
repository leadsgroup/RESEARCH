# Procedure.py 
# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------     
import RCAIDE
from RCAIDE.Framework.Core import Units 
from RCAIDE.Library.Methods.Geometry.Planform              import wing_segmented_planform   
from RCAIDE.Framework.Analyses.Process                     import Process     
from RCAIDE.Library.Methods.Weights.Moment_of_Inertia      import compute_aircraft_moment_of_inertia
from RCAIDE.Library.Methods.Weights.Center_of_Gravity      import compute_vehicle_center_of_gravity 

# Routines  
import Missions 

import numpy as np
from copy import deepcopy

# ############################################################################################################################################################       
# STICK FIXED MISSION 
# ############################################################################################################################################################    

def stick_fixed_stability_and_drag_procedure(): 
    procedure                 = Process()
    procedure.modify_vehicle  = modify_stick_fixed_vehicle
 
    procedure.missions                   = Process()
    procedure.missions.design_mission    = run_stick_fixed_mission     
    procedure.post_process               = longitudinal_static_stability_and_drag_post_process   
        
    return procedure
    
def run_stick_fixed_mission(nexus): 
    results                    = nexus.results
    results.stick_fixed_cruise = nexus.missions.stick_fixed_cruise.evaluate()
    return nexus


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
    vehicle.mass_properties.moments_of_inertia.tensor[0, 1] = 0
    vehicle.mass_properties.moments_of_inertia.tensor[0, 2] = 0
    vehicle.mass_properties.moments_of_inertia.tensor[1, 0] = 0
    vehicle.mass_properties.moments_of_inertia.tensor[1, 2] = 0
    vehicle.mass_properties.moments_of_inertia.tensor[2, 0] = 0
    vehicle.mass_properties.moments_of_inertia.tensor[2, 1] = 0    

    # Update Mission
    AoA_guess      = nexus.missions['stick_fixed_cruise'].segments['cruise'].angle_of_attack 
    Phi_guess      = nexus.missions['stick_fixed_cruise'].segments['cruise'].bank_angle 
    nexus.missions = Missions.stick_fixed_stability_setup(nexus.analyses,weight_analysis.vehicle,nexus.cruise_velocity, nexus.cruise_altitude,AoA_guess,Phi_guess)      
    
    # diff the new data
    vehicle.store_diff() 
    
    return nexus

def longitudinal_static_stability_and_drag_post_process(nexus): 
    '''
    This function analyses and post processes the aircraft at cruise conditions. 
    The objective of is to minimize the drag  of a trimmed aircraft 
    '''
    summary         = nexus.summary 
    vehicle         = nexus.vehicle_configurations.stick_fixed_cruise 
    segment_results = nexus.results.stick_fixed_cruise.segments.cruise
    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------     
    summary.CD              = segment_results.conditions.static_stability.coefficients.drag[0,0]  
    summary.CM_residual     = abs( segment_results.conditions.static_stability.coefficients.M[0,0])
    summary.spiral_criteria = segment_results.conditions.static_stability.spiral_criteria[0,0]
    NP                      = segment_results.conditions.static_stability.neutral_point[0,0]
    cg                      = vehicle.mass_properties.center_of_gravity[0][0]
    MAC                     = vehicle.wings.main_wing.chords.mean_aerodynamic
    summary.static_margin   = (NP - cg)/MAC
    summary.CM_alpha        = segment_results.conditions.static_stability.derivatives.CM_alpha[0,0]
    AoA                     = segment_results.conditions.aerodynamics.angles.alpha[0,0] / Units.degree
 
    summary.F_x_residual      =  segment_results.state.residuals.force_x[0,0] 
    summary.F_y_residual      =  segment_results.state.residuals.force_y[0,0]   
    summary.F_z_residual      =  segment_results.state.residuals.force_z[0,0] 
    summary.M_x_residual      =  segment_results.state.residuals.moment_x[0,0] 
    summary.M_y_residual      =  segment_results.state.residuals.moment_y[0,0] 
    summary.M_z_residual      =  segment_results.state.residuals.moment_z[0,0]  
 
    if np.count_nonzero(vehicle.mass_properties.moments_of_inertia.tensor) > 0:  
        summary.phugoid_damping_ratio       = segment_results.conditions.dynamic_stability.LongModes.phugoidDamping[0,0] 
        summary.short_period_damping_ratio  = segment_results.conditions.dynamic_stability.LongModes.shortPeriodDamping[0,0] 
        summary.dutch_roll_frequency        = segment_results.conditions.dynamic_stability.LatModes.dutchRollFreqHz[0,0]* (2 * np.pi)  
        summary.dutch_roll_damping_ratio    = segment_results.conditions.dynamic_stability.LatModes.dutchRollDamping[0,0]
        summary.spiral_doubling_time        = segment_results.conditions.dynamic_stability.LatModes.spiralTimeDoubleHalf[0,0] 
        print("Angle of Attack            : " + str(AoA))
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



# ############################################################################################################################################################   
# ELEVATOR SIZING
# ############################################################################################################################################################   

def elevator_sizing_setup(): 
    procedure = Process()  
    procedure.missions                   = Process()
    procedure.missions.design_mission    = run_elevator_sizing_mission     
    procedure.post_process               = elevator_sizing_post_process   
    return procedure
 

def run_elevator_sizing_mission(nexus): 
    results                           = nexus.results
    results.elevator_sizing_pull_up   = nexus.missions.elevator_sizing_pull_up.evaluate()
    results.elevator_sizing_push_over = nexus.missions.elevator_sizing_push_over.evaluate()
    return nexus
 

def elevator_sizing_post_process(nexus): 
    '''
    This function analyses and post processes the aircraft at the flight conditions required to size
    the elevator. These conditions are:
    1) Stick pull maneuver with a load factor of 3.0
    2) Stick push maneuver with a load factor of -1
    ''' 
    summary    = nexus.summary   
    vehicle    = nexus.vehicle_configurations.elevator_sizing_pull_up 
    pull_up    = nexus.results.elevator_sizing_pull_up.segments.cruise 
    push_over  = nexus.results.elevator_sizing_push_over.segments.cruise

    g                     = 9.81      
    S                     = vehicle.reference_area 
    m                     = vehicle.mass_properties.max_takeoff  
    q                     = push_over.conditions.freestream.dynamic_pressure[0, 0]   
    CL_pull_man_required  = vehicle.maxiumum_load_factor*m*g/(S*q)  
    CL_push_man_required  = vehicle.minimum_load_factor*m*g/(S*q)
            

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------  
    AoA_pull                          = pull_up.conditions.aerodynamics.angles.alpha[0,0]
    AoA_push                          = push_over.conditions.aerodynamics.angles.alpha[0,0]    
    
    # compute control surface area 
    control_surfaces = ['elevator'] 
    total_control_surface_area    = compute_control_surface_areas(control_surfaces,vehicle)  
    summary.elevator_surface_area =  total_control_surface_area
    
    # compute lift coefficient and residual from target CLs 
    CL_pull_man                  = pull_up.conditions.aerodynamics.coefficients.lift.total[0,0]
    CL_push_man                  = push_over.conditions.aerodynamics.coefficients.lift.total[0,0] 
    summary.CL_pull_residual     = abs(CL_pull_man - CL_pull_man_required)  
    summary.CL_push_residual     = abs(CL_push_man - CL_push_man_required)  
   
    summary.pull_up_CZ_residual     =  pull_up.conditions.static_stability.coefficients.Z[0,0] # pull_up.state.residuals.force_z[0,0]  
    summary.pull_up_CM_residual     =  pull_up.conditions.static_stability.coefficients.M[0,0] # pull_up.state.residuals.moment_y[0,0]    
    summary.push_over_CZ_residual   =  push_over.conditions.static_stability.coefficients.Z[0,0] # push_over.state.residuals.force_z[0,0]  
    summary.push_over_CM_residual   =  push_over.conditions.static_stability.coefficients.M[0,0] # push_over.state.residuals.moment_y[0,0]
      
    # compute elevator deflections 
    elevator_pull_deflection          = pull_up.conditions.control_surfaces.elevator.deflection[0,0]  
    elevator_push_deflection          = push_over.conditions.control_surfaces.elevator.deflection[0,0] 
    summary.elevator_pull_deflection  = abs(elevator_pull_deflection)  
    summary.elevator_push_deflection  = abs(elevator_push_deflection)  

    print("Elevator Area      : " + str(summary.elevator_surface_area))
    print("Aircraft CL Pull   : " + str(CL_pull_man))
    print("CL Pull Residual   : " + str(summary.CL_pull_residual))
    print("Aircraft AoA Pull  : " + str(AoA_pull/Units.degree))
    print("Elevator Pull Defl.: " + str(elevator_pull_deflection/Units.degree)) 
    print("Aircraft CL Push   : " + str(CL_push_man))
    print("CL Push Residual   : " + str(summary.CL_push_residual))
    print("Aircraft AoA Push  : " + str(AoA_push/Units.degree))
    print("Aircraft Push Defl.: " + str(elevator_push_deflection/Units.degree)) 
    print("\n\n")     
         
    return nexus
 


# ############################################################################################################################################################    
# AILERON AND RUDDER SIZING
# ############################################################################################################################################################  

def aileron_rudder_sizing_setup(): 
    procedure = Process()  
    procedure.missions                   = Process()
    procedure.missions.design_mission    = run_aileron_rudder_sizing_mission     
    procedure.post_process               = aileron_rudder_sizing_post_process   
    return procedure   
  
 
def run_aileron_rudder_sizing_mission(nexus): 
    results                     = nexus.results 
    results.roll_maneuver       = nexus.missions.roll_maneuver.evaluate() 
    results.crosswind_maneuver  = nexus.missions.crosswind_maneuver.evaluate()  
    return nexus


 
def aileron_rudder_sizing_post_process(nexus):  
    '''
    This function analyses and post processes the aircraft at the flight conditions required to size
    the aileron and rudder. These conditions are:
    1) A controlled roll at a  rate of 0.07 https://www.conachenuavs.com/downloads/ProAdvice%203%20-%20AILERON%20SIZING.pdf
    2) Trimmed flight in a 20 knot crosswind
    ''' 
    summary               = nexus.summary  
    g                     = 9.81   
    vehicle               = nexus.vehicle_configurations.aileron_rudder_roll_sizing  
    m                     = vehicle.mass_properties.max_takeoff
    S                     = vehicle.reference_area    
    roll                  = nexus.results.roll_maneuver.segments.cruise 
    crosswind             = nexus.results.crosswind_maneuver.segments.cruise
    V                     = roll.conditions.freestream.velocity[0, 0]  
    span                  = vehicle.wings.main_wing.spans.projected
    desired_roll_rate     = 0.07 
     

    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------  

    # ------------------------------------------------------------------------------------------------------------------------  
    # Aileron 
    # compute lift coefficient and residual from target CLs 
    CL_roll_man                   = roll.conditions.aerodynamics.coefficients.lift.total[0,0] 
    CL_p                          = roll.conditions.static_stability.derivatives.CL_p[0,0] 
    CL_delta_a                    = roll.conditions.static_stability.derivatives.CL_a[0,0] 
    CL_crosswind_man              = crosswind.conditions.aerodynamics.coefficients.lift.total[0,0]
    CL_roll_required              = vehicle.maxiumum_load_factor*m*g/(S*roll.conditions.freestream.dynamic_pressure[0,0] ) 
    CL_crosswind_required         = vehicle.minimum_load_factor*m*g/(S*crosswind.conditions.freestream.dynamic_pressure[0,0] )
    summary.CL_roll_residual      = abs(CL_roll_man - CL_roll_required)  
    summary.CL_crosswind_residual = abs(CL_crosswind_man - CL_crosswind_required)  
     
    
    # compute aileron deflections 
    aileron_roll_deflection               = roll.conditions.control_surfaces.aileron.deflection[0,0]   
    aileron_crosswind_deflection          = crosswind.conditions.control_surfaces.aileron.deflection[0,0]  
    summary.aileron_roll_deflection       = abs(aileron_roll_deflection)  
    summary.aileron_crosswind_deflection  = abs(aileron_crosswind_deflection)
 
    summary.crosswind_CY_residual    =  crosswind.conditions.static_stability.coefficients.Y[0,0] 
    summary.crosswind_CZ_residual    =  crosswind.conditions.static_stability.coefficients.Z[0,0]
    summary.crosswind_CL_residual    =  crosswind.conditions.static_stability.coefficients.L[0,0]
    summary.crosswind_CM_residual    =  crosswind.conditions.static_stability.coefficients.M[0,0]
    summary.crosswind_CN_residual    =  crosswind.conditions.static_stability.coefficients.N[0,0]   

    # ------------------------------------------------------------------------------------------------------------------------  
    # Rudder     
    if vehicle.rudder_flag: 
        rudder_roll_deflection               = roll.conditions.control_surfaces.rudder.deflection  
        rudder_crosswind_deflection          = crosswind.conditions.control_surfaces.rudder.deflection 
        summary.rudder_roll_deflection       = abs(rudder_roll_deflection)  
        summary.rudder_crosswind_deflection  = abs(rudder_crosswind_deflection)  
    else:
        rudder_roll_deflection              = 0 
        summary.rudder_roll_deflection      = 0  
        rudder_crosswind_deflection         = 0  
        summary.rudder_crosswind_deflection = 0 
         
    roll_rate =  - (CL_delta_a / CL_p) * aileron_roll_deflection * (2 * V * span)
        
    # compute control surface area 
    control_surfaces = ['aileron','rudder'] 
    total_control_surface_area = compute_control_surface_areas(control_surfaces,vehicle)   
    summary.aileron_rudder_surface_area =  total_control_surface_area
    summary.roll_rate_residual          = (roll_rate - desired_roll_rate)

    print("Total Rudder Aileron Surface Area : " + str(summary.aileron_rudder_surface_area)) 
    print("Roll Rate                         : " + str(roll_rate)) 
    print("Aileron Roll Defl                 : " + str(aileron_roll_deflection/Units.degree)) 
    print("Rudder Roll Defl                  : " + str(rudder_roll_deflection/Units.degree))  
    print("Aileron Crosswind Defl            : " + str(aileron_crosswind_deflection/Units.degree)) 
    print("Rudder  Crosswind Defl            : " + str(rudder_crosswind_deflection/Units.degree )) 
    print("\n\n")     
  
    return nexus     
 

# ############################################################################################################################################################   
# FLAP SIZING
# ############################################################################################################################################################  

def flap_sizing_setup(): 
    procedure = Process()    
    procedure.missions                   = Process()
    procedure.missions.design_mission    = run_flap_sizing_mission    
    procedure.post_process               = flap_sizing_post_process 
    return procedure  
 
def run_flap_sizing_mission(nexus): 
    results                         = nexus.results 
    results.flap_sizing_flaps_up    = nexus.missions.flap_sizing_flaps_up.evaluate()  
    results.flap_sizing_flaps_down  = nexus.missions.flap_sizing_flaps_down.evaluate()  
    return nexus

  
def flap_sizing_post_process(nexus): 
    '''
    This function analyses and post processes the aircraft at the flight conditions required to size
    the flap. These conditions are:
    1) A comparison of clean and deployed flap at 12 deg. angle of attack
    ''' 
    summary   = nexus.summary   
    vehicle   = nexus.vehicle_configurations.flap_sizing_flaps_up   

    no_flap  = nexus.results.flap_sizing_flaps_up.segments.cruise 
    w_flap   = nexus.results.flap_sizing_flaps_down.segments.cruise 
    # ------------------------------------------------------------------------------------------------------------------------  
    # Post Process Results 
    # ------------------------------------------------------------------------------------------------------------------------       
    # critera    
    CL_12_deg_no_flap    = no_flap.conditions.aerodynamics.coefficients.lift.total[0,0]  
    CL_12_deg_flap       = w_flap.conditions.aerodynamics.coefficients.lift.total[0,0]
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

# ############################################################################################################################################################  
# SUPPLEMENTAL FUNCTIONS 
# ############################################################################################################################################################  
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