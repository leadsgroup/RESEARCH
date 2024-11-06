# Optimize.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------  
import RCAIDE
from RCAIDE.Framework.Core import Units, Data 
from RCAIDE.Framework.Optimization.Packages.scipy import scipy_setup
from RCAIDE.Framework.Optimization.Common   import Nexus
from RCAIDE.Library.Methods.Weights.Moment_of_Inertia.compute_aircraft_moment_of_inertia import compute_aircraft_moment_of_inertia
from RCAIDE.Library.Methods.Weights.Center_of_Gravity     import compute_vehicle_center_of_gravity

# python imports 
import numpy  as np
import pandas as pd
import time

# local imports 
import Vehicles_Single_Engine_Piston_Baseline
import Missions
import Procedure
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------        
#  Main function
# ----------------------------------------------------------------------  

def main(): 
    vehicle = Vehicles_Single_Engine_Piston_Baseline.vehicle_setup()
    
    # ----------------------------------------------------------------------        
    # Step 1: Weight and moment of inertia estimation
    # ----------------------------------------------------------------------    
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    weight_analysis.method                        = 'Raymer'
    weight_analysis.settings.use_max_fuel_weight  = True 
    weight_analysis.evaluate() 
    
    # Calculate CG location
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location                                   = vehicle.mass_properties.center_of_gravity
    
    # Calculate aircraft MOI
    compute_aircraft_moment_of_inertia(vehicle, CG_location)
    
    # ----------------------------------------------------------------------        
    # Step 2-3: Control surface sizing & Stability analysis
    # ----------------------------------------------------------------------     
    
    
    '''
    STICK FIXED (STATIC STABILITY AND DRAG OTIMIZATION
    '''
    ti = time.time()  
    solver_name       = 'SLSQP' 
    planform_optimization_problem = stick_fixed_stability_and_drag_optimization_setup()
    output_stick_fixed = scipy_setup.SciPy_Solve(planform_optimization_problem,solver=solver_name, sense_step = 1E-2, tolerance = 1E-3)  
    print (output_stick_fixed)    
    tf           = time.time()
    elapsed_time = round((tf-ti)/60,2)
    print('Stick Fixed Stability and Drag Otimization Simulation Time: ' + str(elapsed_time))  
    
    df = pd.DataFrame({
    'mw_root_twist': [output_stick_fixed[0]], 
    'mw_tip_twist': [output_stick_fixed[1]], 
    'mw_dihedral': [output_stick_fixed[2]], 
    'vt_span': [output_stick_fixed[3]], 
    'vt_taper': [output_stick_fixed[4]], 
    'c_g_x': [output_stick_fixed[5]],
    'CD': [planform_optimization_problem.summary['CD']],
    'CM_residual': [planform_optimization_problem.summary['CM_residual']],
    'spiral_criteria': [planform_optimization_problem.summary['spiral_criteria']],
    'static_margin': [planform_optimization_problem.summary['static_margin']],
    'CM_alpha': [planform_optimization_problem.summary['CM_alpha']],
    'phugoid_damping_ratio': [planform_optimization_problem.summary['phugoid_damping_ratio']],
    'short_period_damping_ratio': [planform_optimization_problem.summary['short_period_damping_ratio']],
    'dutch_roll_frequency': [planform_optimization_problem.summary['dutch_roll_frequency']],
    'dutch_roll_damping_ratio': [planform_optimization_problem.summary['dutch_roll_damping_ratio']],
    'spiral_doubling_time': [planform_optimization_problem.summary['spiral_doubling_time']]}) 
    
    df.to_csv('Output_Single_Engine_Piston_Baseline.csv')
    #np.savetxt(f'{save_filename}.txt', np.column_stack((output)), delimiter=',', header='mw_root_twist,mw_tip_twist,mw_dihedral,vt_span,vt_taper,c_g_x', comments='')  
        
    
    '''
    ELEVATOR SIZING
    '''      
    # define vehicle for elevator sizing 
    optimized_vehicle_v1                             = planform_optimization_problem.vehicle_configurations.stick_fixed_cruise 
    optimized_vehicle_v1.maximum_elevator_deflection = 30*Units.degrees 
    optimized_vehicle_v1.maxiumum_load_factor        = 3.0
    optimized_vehicle_v1.minimum_load_factor         = -1
    
    ti = time.time()   
    solver_name       = 'SLSQP'  
    elevator_sizing_optimization_problem = elevator_sizing_optimization_setup(optimized_vehicle_v1)
    output_elevator_sizing = scipy_setup.SciPy_Solve(elevator_sizing_optimization_problem,solver=solver_name, sense_step = 1E-2, tolerance = 1E-3) 
    print (output_elevator_sizing)     
    tf           = time.time()
    elapsed_time = round((tf-ti)/60,2)
    print('Elevator Sizing Simulation Time: ' + str(elapsed_time))   
     
    '''
    AILERON AND RUDDER SIZING
    '''      
    # define vehicle for aileron and rudder sizing
    optimized_vehicle    = elevator_sizing_optimization_problem.vehicle_configurations.elevator_sizing 
    optimized_vehicle.rudder_flag                       = True 
    optimized_vehicle.maximum_aileron_rudder_deflection = 30*Units.degrees 
    optimized_vehicle.crosswind_velocity                = 20 * Units.knots

    ti = time.time()   
    solver_name       = 'SLSQP'  
    aileron_rudder_sizing_optimization_problem = aileron_rudder_sizing_optimization_setup(optimized_vehicle)
    output_aileron_and_rudder_sizing = scipy_setup.SciPy_Solve(aileron_rudder_sizing_optimization_problem,solver=solver_name, sense_step = 1E-2, tolerance = 1E-3) 
    print (output_aileron_and_rudder_sizing)     
    tf           = time.time()
    elapsed_time = round((tf-ti)/60,2)
    print('Aileron and Rudder Sizing Simulation Time: ' + str(elapsed_time))   

    '''
    FLAP SIZING
    '''      
    # define vehicle for flap sizing     
    optimized_vehicle_v3 = aileron_rudder_sizing_optimization_problem.vehicle_configurations.aileron_rudder_sizing
    optimized_vehicle_v3.maximum_flap_deflection = 40*Units.degrees
    
    ti = time.time()   
    solver_name       = 'SLSQP'  
    flap_sizing_optimization_problem = flap_sizing_optimization_setup(optimized_vehicle_v3)
    output_flap_sizing = scipy_setup.SciPy_Solve(flap_sizing_optimization_problem,solver=solver_name, sense_step = 1E-2, tolerance = 1E-3) 
    print (output_flap_sizing)     
    tf           = time.time()
    elapsed_time = round((tf-ti)/60,2)
    print('Flap Sizing Simulation Time: ' + str(elapsed_time))   

    '''
    PRINT VEHICLE CONTROL SURFACES
    '''          
    optimized_vehicle_v4  = flap_sizing_optimization_problem.vehicle_configurations.flap_sizing 
    print_vehicle_control_surface_geoemtry(optimized_vehicle_v4)
    
    # ----------------------------------------------------------------------        
    # Step 4: Mission simulation
    # ----------------------------------------------------------------------     
    
    
    
    
    return
  
def stick_fixed_stability_and_drag_optimization_setup(): 
    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #                 [ tag                       , initial,  (lb , ub) , scaling , units ]  
    problem.inputs = np.array([          
                  #[ 'mw_taper'                    , 0.54  , 0.4  , 0.8   , 1.0  ,  1*Units.less],    
                  #[ 'mw_area'                     , 17.11 , 15.0 , 19.0  , 100. ,  1*Units.meter**2],    
                  #[ 'mw_AR'                       , 6.04  , 5.0  , 8.0   , 100  ,  1*Units.less],         
                  [ 'mw_root_twist'               , 3.0   , 0.0  , 5.0   , 1. ,  1*Units.degree], 
                  [ 'mw_tip_twist'                , -1.0  , -4.0 , 0.0   , 1.  ,  1*Units.degree],
                  [ 'mw_dihedral'                 , 7.5   , 0.0  , 8.0   , 1.  ,  1*Units.degree], 
                  [ 'vt_span'                     , 1.4816, 1.0  , 2     , 1.  ,  1*Units.meter],
                  [ 'vt_taper'                    , 0.4820, 0.2  , 0.8   , 1.  ,  1*Units.meter],  
                  #[ 'mw_sweep'                    , 15  , 0.0  , 20   , 1.  ,  1*Units.degree], 
                  #[ 'hs_AR'                       , 4.0   , 3.0  , 5.0   , 10.  ,  1*Units.less], 
                  #[ 'hs_area'                     , 4.0   , 3.0  , 5.0   , 10.  ,  1*Units.meter**2],   
                  #[ 'hs_root_twist'               , 0.0   , -3.0 , 3.0   , 10.  , 1*Units.degree],
                  #[ 'hs_tip_twist'                , 0.0   , -3.0 , 3.0   , 10.  , 1*Units.degree], 
                  [ 'c_g_x'                       , 2.239696797  , 1  , 3, 1    ,  1*Units.less],
                  
    ],dtype=object)   

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                 [  'CD'  ,  1.0  ,    1*Units.less] 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ] See http://everyspec.com/MIL-SPECS/MIL-SPECS-MIL-F/MIL-F-8785C_5295/
    # Level 1, Category B
    problem.constraints = np.array([
        [ 'CM_residual'               ,   '<' ,   1E-2  ,   1E-2  , 1*Units.less], # close to zero 2 works 
        [ 'static_margin'             ,   '>' ,   0.1   ,   0.1   , 1*Units.less],  # checked 
        [ 'CM_alpha'                  ,   '<' ,   0.0   ,   1.0   , 1*Units.less],  # checked 
        [ 'phugoid_damping_ratio'     ,   '>' ,   0.04  ,   1.0   , 1*Units.less],  # checked 
        [ 'short_period_damping_ratio',   '<' ,   2.0   ,   1.0   , 1*Units.less],  # checked 
        [ 'short_period_damping_ratio',   '>' ,   0.3   ,   1.0   , 1*Units.less], # checked    
        [ 'dutch_roll_frequency'      ,   '>' ,   0.4   ,   1.0   , 1*Units.less],  # checked   frequency in rad/sec
        [ 'dutch_roll_damping_ratio'  ,   '>' ,   0.08  ,   1.0   , 1*Units.less],  # checked   
        [ 'spiral_doubling_time'      ,   '>' ,   20.0  ,   1.0  , 1*Units.less],  # checked   
        [ 'spiral_criteria'           ,   '>' ,   1.0   ,   1.0   , 1*Units.less],  # checked   
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [ 
        [ 'CD'                                , 'summary.CD' ],
        [ 'CM_residual'                       , 'summary.CM_residual' ],  
        [ 'CM_alpha'                          , 'summary.CM_alpha' ],    
        [ 'static_margin'                     , 'summary.static_margin' ], 
        [ 'phugoid_damping_ratio'             , 'summary.phugoid_damping_ratio' ],  
        [ 'short_period_damping_ratio'        , 'summary.short_period_damping_ratio' ],  
        [ 'dutch_roll_frequency'              , 'summary.dutch_roll_frequency' ],  
        [ 'dutch_roll_damping_ratio'          , 'summary.dutch_roll_damping_ratio' ],  
        [ 'spiral_doubling_time'              , 'summary.spiral_doubling_time' ],   
        [ 'spiral_criteria'                   , 'summary.spiral_criteria' ],      
        #[ 'mw_area'                           , 'vehicle_configurations.*.wings.main_wing.areas.reference'],
        #[ 'mw_taper'                          , 'vehicle_configurations.*.wings.main_wing.taper'],
        #[ 'mw_AR'                             , 'vehicle_configurations.*.wings.main_wing.aspect_ratio'],         
        [ 'mw_root_twist'                     , 'vehicle_configurations.*.wings.main_wing.twists.root' ], 
        [ 'mw_tip_twist'                      , 'vehicle_configurations.*.wings.main_wing.twists.tip'  ],     
        [ 'mw_dihedral'                       , 'vehicle_configurations.*.wings.main_wing.dihedral'  ],     
        [ 'vt_span'                            , 'vehicle_configurations.*.wings.vertical_stabilizer.spans.projected'],
        [ 'vt_taper'                           , 'vehicle_configurations.*.wings.vertical_stabilizer.taper'],
        #[ 'mw_sweep'                          , 'vehicle_configurations.*.wings.main_wing.sweeps.quarter_chord'  ],         
        #[ 'hs_AR'                             , 'vehicle_configurations.*.wings.horizontal_stabilizer.aspect_ratio'],     
        #[ 'hs_area'                           , 'vehicle_configurations.*.wings.horizontal_stabilizer.areas.reference'],
        #[ 'hs_taper'                          , 'vehicle_configurations.*.wings.horizontal_stabilizer.taper'],  
        #[ 'hs_root_twist'                     , 'vehicle_configurations.*.wings.horizontal_stabilizer.twists.root' ], 
        #[ 'hs_tip_twist'                      , 'vehicle_configurations.*.wings.horizontal_stabilizer.twists.tip'  ], 
        #[ 'hs_dihedral'                       , 'vehicle_configurations.*.wings.horizontal_stabilizer.dihedral'  ],
        [ 'c_g_x'                             , 'vehicle_configurations.*.mass_properties.center_of_gravity[0][0]'  ], 
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles_Single_Engine_Piston_Baseline.stick_fixed_stability_setup()
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = None 
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.stick_fixed_stability_setup(nexus.analyses,nexus.vehicle_configurations.stick_fixed_cruise)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.stick_fixed_stability_and_drag_procedure()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus 

def elevator_sizing_optimization_setup(vehicle):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    problem.inputs = np.array([            
                  [ 'hs_elevator_chord_fraction' , 0.2    , 0.1  , 0.3   ,  1.0  ,  1*Units.less],
                  [ 'hs_elevator_span_frac_start', 0.25   , 0.05 , 0.5   ,  1.0 ,  1*Units.less], 
                  [ 'hs_elevator_span_frac_end'  , 0.75   , 0.6  , 0.95  ,  1.0  ,  1*Units.less],     
                  
    ],dtype=object)   

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'elevator_surface_area', 1. , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([ 
        [ 'elevator_push_deflection_residual'           ,   '>' ,  0  , 1.0   , 1*Units.less], 
        [ 'elevator_pull_deflection_residual'           ,   '>' ,  0  , 1.0   , 1*Units.less], 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [ 
        [ 'elevator_surface_area'             , 'summary.elevator_surface_area' ], 
        [ 'elevator_push_deflection_residual' , 'summary.elevator_push_deflection_residual' ],   
        [ 'elevator_pull_deflection_residual' , 'summary.elevator_pull_deflection_residual' ],     
        [ 'hs_elevator_chord_fraction'        , 'vehicle_configurations.*.wings.horizontal_stabilizer.control_surfaces.elevator.chord_fraction'],    
        [ 'hs_elevator_span_frac_start'       , 'vehicle_configurations.*.wings.horizontal_stabilizer.control_surfaces.elevator.span_fraction_start'],    
        [ 'hs_elevator_span_frac_end'         , 'vehicle_configurations.*.wings.horizontal_stabilizer.control_surfaces.elevator.span_fraction_end'],  
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles_Single_Engine_Piston_Baselineelevator_sizing_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = None 
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.elevator_sizing_setup(nexus.analyses,nexus.vehicle_configurations.elevator_sizing)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.elevator_sizing_setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  


 
 
def aileron_rudder_sizing_optimization_setup(vehicle):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    if vehicle.rudder_flag:
        problem.inputs = np.array([             
                      [ 'mw_aileron_chord_fraction'  , 0.2    , 0.15 , 0.3  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_start' , 0.75   , 0.55 , 0.8  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_end'   , 0.9    , 0.85 , 0.95 ,  1.0 ,  1*Units.less],  
                      [ 'vs_rudder_chord_fraction'   , 0.2    , 0.15 , 0.3  ,  1.0 ,  1*Units.less],
                      [ 'vs_rudder_span_frac_start'  , 0.25   , 0.05 , 0.35 ,  1.0 ,  1*Units.less],
                      [ 'vs_rudder_span_frac_end'    , 0.75   , 0.5  , 0.95 ,  1.0 ,  1*Units.less] ],dtype=object)   
    else:
        problem.inputs = np.array([             
                      [ 'mw_aileron_chord_fraction'  , 0.2    , 0.15 , 0.3  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_start' , 0.75   , 0.55 , 0.8  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_end'   , 0.9    , 0.85 , 0.95 ,  1.0 ,  1*Units.less],   
                      
        ],dtype=object)      

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'aileron_rudder_surface_area', 1. , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    if vehicle.rudder_flag:
        problem.constraints = np.array([
            [ 'aileron_roll_deflection_residual'       ,   '>' ,  0  ,  1.0    , 1*Units.less], 
            [ 'rudder_roll_deflection_residual'        ,   '>' ,  0  ,  1.0    , 1*Units.less],  
            [ 'aileron_crosswind_deflection_residual'  ,   '>' ,  0  ,  1.0    , 1*Units.less], 
            [ 'rudder_crosswind_deflection_residual'   ,   '>' ,  0  ,  1.0    , 1*Units.less],  
        ],dtype=object)
    else:
        problem.constraints = np.array([
            [ 'aileron_roll_deflection_residual'       ,   '>' ,  0  ,  1.0    , 1*Units.less],  
            [ 'aileron_crosswind_deflection_residual'  ,   '>' ,  0  ,  1.0    , 1*Units.less],  
        ],dtype=object)
        
        
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    if vehicle.rudder_flag:
        problem.aliases = [ 
            [ 'aileron_rudder_surface_area'            , 'summary.aileron_rudder_surface_area' ],  
            [ 'aileron_roll_deflection_residual'       , 'summary.aileron_roll_deflection_residual' ],  
            [ 'rudder_roll_deflection_residual'        , 'summary.rudder_roll_deflection_residual' ], 
            [ 'aileron_crosswind_deflection_residual'  , 'summary.aileron_crosswind_deflection_residual' ],      
            [ 'rudder_crosswind_deflection_residual'   , 'summary.rudder_crosswind_deflection_residual' ],         
            [ 'mw_aileron_chord_fraction'              , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.chord_fraction'],  
            [ 'mw_aileron_span_frac_start'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_start'],    
            [ 'mw_aileron_span_frac_end'               , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_end'],     
            [ 'vs_rudder_chord_fraction'               , 'vehicle_configurations.*.wings.vertical_stabilizer.control_surfaces.rudder.chord_fraction'],    
            [ 'vs_rudder_span_frac_start'              , 'vehicle_configurations.*.wings.vertical_stabilizer.control_surfaces.rudder.span_fraction_start'],    
            [ 'vs_rudder_span_frac_end'                , 'vehicle_configurations.*.wings.vertical_stabilizer.control_surfaces.rudder.span_fraction_end']]      
    else:
        problem.aliases = [ 
            [ 'aileron_rudder_surface_area'            , 'summary.aileron_rudder_surface_area' ],  
            [ 'aileron_roll_deflection_residual'       , 'summary.aileron_roll_deflection_residual' ],          
            [ 'aileron_crosswind_deflection_residual'  , 'summary.aileron_crosswind_deflection_residual' ],   
            [ 'mw_aileron_chord_fraction'              , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.chord_fraction'],  
            [ 'mw_aileron_span_frac_start'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_start'],    
            [ 'mw_aileron_span_frac_end'               , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_end']] 
        
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles_Single_Engine_Piston_Baseline.aileron_rudder_sizing_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = None
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.aileron_rudder_sizing_setup(nexus.analyses,nexus.vehicle_configurations.aileron_rudder_sizing)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.aileron_rudder_sizing_setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  


def flap_sizing_optimization_setup(optimized_vehicle):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    problem.inputs = np.array([           
                  [ 'mw_flap_chord_fraction'     , 0.2    , 0.15 , 0.4  ,  1.0 ,  1*Units.less],
                  [ 'mw_flap_span_frac_start'    , 0.2    , 0.05 , 0.25 ,  1.0 ,  1*Units.less],
                  [ 'mw_flap_span_frac_end'      , 0.4    , 0.3  , 0.5  ,  1.0 ,  1*Units.less],   
                  
    ],dtype=object)   

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'flap_surface_area', 1. , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([
        [ 'flap_criteria'    ,   '>' ,  0.   ,  1.0   , 1*Units.less],  
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [ 
        [ 'flap_surface_area'                 , 'summary.flap_surface_area' ], 
        [ 'flap_criteria'                     , 'summary.flap_criteria' ],   
        [ 'mw_flap_chord_fraction'            , 'vehicle_configurations.*.wings.main_wing.control_surfaces.flap.chord_fraction'],    
        [ 'mw_flap_span_frac_start'           , 'vehicle_configurations.*.wings.main_wing.control_surfaces.flap.span_fraction_start'],    
        [ 'mw_flap_span_frac_end'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.flap.span_fraction_end'],    
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles_Single_Engine_Piston_Baseline.flap_sizing_setup(optimized_vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = None
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.flap_sizing_setup(nexus.analyses,nexus.vehicle_configurations.flap_sizing)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.flap_sizing_setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  
    
def print_vehicle_control_surface_geoemtry(vehicle): 
   
    for wing in vehicle.wings:
        if 'control_surfaces' in wing:  
            for CS in wing.control_surfaces:  
                print('Wing                : ' + wing.tag)
                print('Control Surface     : ' + CS.tag)
                print('Span Fraction Start : ' + str(CS.span_fraction_start))
                print('Span Fraction End   : ' + str(CS.span_fraction_end)) 
                print('Chord Fraction      : ' + str(CS.chord_fraction)) 
                print("\n\n")     

    return

if __name__ == '__main__':
    main()
    plt.show()
