# Optimize.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------  
import RCAIDE
from   RCAIDE.Framework.Core                        import Units, Data 
from   RCAIDE.Framework.Optimization.Packages.scipy import scipy_setup
from   RCAIDE.Framework.Optimization.Common         import Nexus
from   save_data                      import save_data    
from RCAIDE import  load 
from RCAIDE import  save

# python imports 
import numpy as np
import time
import  pickle
from copy import  deepcopy

# local imports 
import Vehicles
import Missions
import Analyses
import Procedure
import matplotlib.pyplot as plt 
# ----------------------------------------------------------------------        
#   Run the whole thing
# ----------------------------------------------------------------------  
def size_control_surfaces(CG_bat_1, CG_bat_2, vehicle, cruise_velocity = 120 * Units['mph'], cruise_altitude= 5000*Units.feet):

    ti_0 = time.time()
    
    #'''
    #STICK FIXED (STATIC STABILITY AND DRAG OTIMIZATION
    #'''
    #ti   = time.time()   
    #planform_optimization_problem = stick_fixed_stability_and_drag_optimization_setup(vehicle,cruise_velocity,cruise_altitude)
    #output_stick_fixed = scipy_setup.SciPy_Solve(planform_optimization_problem,solver='SLSQP', sense_step = 1E-3, tolerance = 1E-3)   
    #print (output_stick_fixed)    
    #tf           = time.time()
    #elapsed_time_stick_fixed = round((tf-ti)/60,2)
    #print('Stick Fixed Stability and Drag Otimization Simulation Time: ' + str(elapsed_time_stick_fixed))
    
    
    #'''
    #ELEVATOR SIZING (7 mins)
    #'''
    #ti = time.time()    
    #optimized_vehicle                             = planform_optimization_problem.vehicle_configurations.stick_fixed_cruise 
    #save(optimized_vehicle, 'optimized_vehicle_stick_fixed', pickle_format=True) 
    #optimized_vehicle.maxiumum_load_factor        = 3.0
    #optimized_vehicle.minimum_load_factor         = -1 
    #elevator_sizing_optimization_problem = elevator_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude)
    #output_elevator_sizing = scipy_setup.SciPy_Solve(elevator_sizing_optimization_problem,solver='SLSQP', sense_step = 1E-2, tolerance = 1E-2) 
    #print (output_elevator_sizing)     
    #tf           = time.time()
    #elapsed_time_elevator_sizing = round((tf-ti)/60,2)
    #print('Elevator Sizing Simulation Time: ' + str(elapsed_time_elevator_sizing))
   
    
    #'''
    #AILERON AND RUDDER SIZING (22 mins)
    #'''       
    #ti = time.time()  
    #optimized_vehicle                    = elevator_sizing_optimization_problem.vehicle_configurations.elevator_sizing_pull_up # 
    #save(optimized_vehicle, 'optimized_vehicle_ele', pickle_format=True)  
    #optimized_vehicle.rudder_flag        = True  
    #optimized_vehicle.crosswind_velocity = 20 * Units.knots 
    #aileron_rudder_sizing_optimization_problem = aileron_rudder_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude)
    #output_aileron_and_rudder_sizing = scipy_setup.SciPy_Solve(aileron_rudder_sizing_optimization_problem,solver='SLSQP', sense_step = 1E-3, tolerance = 1E-3) 
    #print (output_aileron_and_rudder_sizing)     
    #tf           = time.time()
    #elapsed_time_aileron_and_rudder_sizing = round((tf-ti)/60,2)
    #print('Aileron and Rudder Sizing Simulation Time: ' + str(elapsed_time_aileron_and_rudder_sizing))   
    
    
    #'''
    #AILERON AND RUDDER SIZING OEI
    #'''     
    #ti = time.time()      
    #optimized_vehicle                              = aileron_rudder_sizing_optimization_problem.vehicle_configurations.aileron_rudder_roll_sizing  
    #save(optimized_vehicle, 'optimized_vehicle_ail_rud', pickle_format=True) 
    optimized_vehicle =  load('optimized_vehicle_ail_rud', pickle_format=True) 
    optimized_vehicle.rudder_flag                  = True    
    aileron_rudder_oei_sizing_optimization_problem = aileron_rudder_oei_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude)
    output_aileron_and_rudder_oei_sizing = scipy_setup.SciPy_Solve(aileron_rudder_oei_sizing_optimization_problem,solver='SLSQP', sense_step = 1E-2, tolerance = 1E-2) # sense_step = 1E-2, tolerance = 1E-2
    print (output_aileron_and_rudder_oei_sizing)     
    tf   = time.time()
    elapsed_time_aileron_and_rudder_oei_sizing = round((tf-ti)/60,2)
    print('Aileron and Rudder Sizing Simulation Time: ' + str(elapsed_time_aileron_and_rudder_oei_sizing))   
     
     
    '''
    FLAP SIZING
    ''' 
    ti = time.time()     
    optimized_vehicle = aileron_rudder_oei_sizing_optimization_problem.vehicle_configurations.aileron_rudder_oei_sizing 
    save(optimized_vehicle, 'optimized_vehicle_ail_rud_oei', pickle_format=True) 
    flap_sizing_optimization_problem = flap_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude)
    output_flap_sizing = scipy_setup.SciPy_Solve(flap_sizing_optimization_problem,solver='SLSQP', sense_step = 1E-1, tolerance = 1E-1) 
    print (output_flap_sizing)     
    tf   = time.time()
    elapsed_time_flap_sizing = round((tf-ti)/60,2)
    print('Flap Sizing Simulation Time: ' + str(elapsed_time_flap_sizing))   
    
     
    '''
    ONE ENGINE INOPERATIVE = HOVER 
    '''     
    ti   = time.time()   
    optimized_vehicle  = flap_sizing_optimization_problem.vehicle_configurations.flap_sizing  
    save(optimized_vehicle, 'optimized_vehicle_flap', pickle_format=True) 
    hover_oei_optimization_problem = hover_oei_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude)
    output_hover_oei = scipy_setup.SciPy_Solve(hover_oei_optimization_problem,solver='SLSQP', sense_step = 1E-5, tolerance = 1E-4)  # -3, -3 works 
    print (output_hover_oei)    
    tf   = time.time()
    elapsed_time_hover_oei = round((tf-ti)/60,2)
    print('Hover OEI Simulation Time: ' + str(elapsed_time_hover_oei))
    

    '''
    PRINT VEHICLE CONTROL SURFACES
    '''          
    optimized_vehicle  = hover_oei_optimization_problem.vehicle_configurations.hover_prop_rotor_oei  
    save(optimized_vehicle, 'optimized_vehicle_hover_oei', pickle_format=True) 
    print_vehicle_control_surface_geoemtry(optimized_vehicle)
     
    tf_0           = time.time()
    total_elapsed_time = round((tf_0-ti_0)/60,2)    
    print('Total Control Surface Sizing Time: ' + str(total_elapsed_time))    
    
    # THIS CAN BE REWRITEN, need to store if optimization works or not as well
    #Optimization_Data_2_CSV(CG_bat_1, 
                            #CG_bat_2,
                            #optimized_vehicle_v4.tag, 
                            #output_stick_fixed, 
                            #output_elevator_sizing, 
                            #output_aileron_and_rudder_sizing, 
                            #output_flap_sizing, 
                            #planform_optimization_problem, 
                            #elevator_sizing_optimization_problem, 
                            #aileron_rudder_sizing_optimization_problem, 
                            #flap_sizing_optimization_problem,
                            #elapsed_time_stick_fixed,
                            #elapsed_time_elevator_sizing,
                            #elapsed_time_aileron_and_rudder_sizing,
                            #elapsed_time_flap_sizing)
    
    # Save Vehicle!!!
    
    return
  
def stick_fixed_stability_and_drag_optimization_setup(vehicle,cruise_velocity,cruise_altitude): 
    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude     

    # -------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------- 
    problem.inputs = np.array([       

    #                [ tag             , initial  ,    (lb  ,    ub)  , scaling ,      units ]  
                  [ 'AoA'             , 5       , 0       ,  12     ,  10      ,  1*Units.degree],
                  [ 'mw_span'         , 11.82855 , 10      , 13      , 10      ,  1*Units.less],    
                  [ 'mw_AR'           , 8.95198  , 7       , 10      , 10.     ,  1*Units.meter**2],        
                  [ 'mw_root_twist'   , 4        , 3.0     , 5.0     , 1.      ,  1*Units.degree], 
                  [ 'mw_tip_twist'    , 0        , -1.0    , 1.0     , 1.      ,  1*Units.degree], 
                  [ 'vt_span'         , 1.4816   , 1.0     , 2       , 1.      ,  1*Units.meter],  
                  [ 'vt_AR'           , 1.8874   , 1       , 2       , 1.      ,  1*Units.meter**2],    
                  [ 'ht_span'         , 4.0      , 3.0     , 5.0     , 10.     ,  1*Units.less], 
                  [ 'ht_AR'           , 5.3333   , 3.0     , 6.0     , 10.     ,  1*Units.meter**2], 
                  [ 'ht_tip_twist'    , 0         , -2      , 4       , 1       ,  1*Units.degree],                   
                  
    ],dtype=object)       

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                 [  'CD'  ,  1   ,    1*Units.less] 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ] See http://everyspec.com/MIL-SPECS/MIL-SPECS-MIL-F/MIL-F-8785C_5295/
    # Level 1, Category B
    problem.constraints = np.array([
        [ 'F_z_residual'               ,   '<' ,   1E-3  ,   1E-3  , 1*Units.less], # close to zero 2 works 
        [ 'M_y_residual'               ,   '<' ,   1E-3  ,   1E-3  , 1*Units.less], # close to zero 2 works 
        [ 'static_margin'              ,   '>' ,   0.1   ,   0.1   , 1*Units.less],  # checked 
        [ 'CM_alpha'                   ,   '<' ,   0.0   ,   1.0   , 1*Units.less],  # checked 
        [ 'phugoid_damping_ratio'      ,   '>' ,   0.04  ,   1.0   , 1*Units.less],  # checked 
        [ 'short_period_damping_ratio' ,   '<' ,   2.0   ,   1.0   , 1*Units.less],  # checked 
        [ 'short_period_damping_ratio' ,   '>' ,   0.3   ,   1.0   , 1*Units.less], # checked    
        #[ 'dutch_roll_frequency'      ,   '>' ,   0.4   ,   1.0   , 1*Units.less],  # checked   frequency in rad/sec
        #[ 'dutch_roll_damping_ratio'  ,   '>' ,   0.08  ,   1.0   , 1*Units.less],  # checked   
        #[ 'spiral_doubling_time'      ,   '>' ,   20.0  ,   1.0   , 1*Units.less],  # checked   
        #[ 'spiral_criteria'           ,   '>' ,   1.0   ,   1.0   , 1*Units.less],  # checked   
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [ 
        [ 'F_z_residual'                      , 'summary.F_z_residual' ],
        [ 'CD'                                , 'summary.CD' ],
        [ 'M_y_residual'                      , 'summary.M_y_residual' ],  
        [ 'CM_alpha'                          , 'summary.CM_alpha' ],    
        [ 'static_margin'                     , 'summary.static_margin' ], 
        [ 'phugoid_damping_ratio'             , 'summary.phugoid_damping_ratio' ],  
        [ 'short_period_damping_ratio'        , 'summary.short_period_damping_ratio' ],  
        [ 'dutch_roll_frequency'              , 'summary.dutch_roll_frequency' ],  
        [ 'dutch_roll_damping_ratio'          , 'summary.dutch_roll_damping_ratio' ],  
        [ 'spiral_doubling_time'              , 'summary.spiral_doubling_time' ],   
        [ 'spiral_criteria'                   , 'summary.spiral_criteria' ],       
        [ 'AoA'                               , 'missions.stick_fixed_cruise.segments.cruise.angle_of_attack' ],
        [ 'mw_span'                           , 'vehicle_configurations.*.wings.main_wing.spans.projected'],
        [ 'mw_AR'                             , 'vehicle_configurations.*.wings.main_wing.aspect_ratio'],         
        [ 'mw_root_twist'                     , 'vehicle_configurations.*.wings.main_wing.twists.root' ], 
        [ 'mw_tip_twist'                      , 'vehicle_configurations.*.wings.main_wing.twists.tip'  ],     
        [ 'vt_span'                           , 'vehicle_configurations.*.wings.vertical_tail.spans.projected'],
        [ 'vt_AR'                             , 'vehicle_configurations.*.wings.vertical_tail.aspect_ratio'],         
        [ 'ht_AR'                             , 'vehicle_configurations.*.wings.horizontal_tail.aspect_ratio'],     
        [ 'ht_span'                           , 'vehicle_configurations.*.wings.horizontal_tail.spans.projected'], 
        [ 'ht_tip_twist'                      , 'vehicle_configurations.*.wings.horizontal_tail.twists.tip'], 
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.stick_fixed_stability_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.stick_fixed_stability_setup(nexus.analyses,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.stick_fixed_stability_and_drag_procedure()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus 

def elevator_sizing_optimization_setup(vehicle,cruise_velocity,cruise_altitude):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude    

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    problem.inputs = np.array([            
                  [ 'pull_up_AoA'                , 5     , -20  ,  20    ,  1.0 ,  1*Units.degree],
                  [ 'push_over_AoA'              , -5    , -20  ,  20    ,  1.0 ,  1*Units.degree],
                  [ 'elevator_push_deflection'   ,  10   , -30  ,  30    ,  1.0 ,  1*Units.degree],
                  [ 'elevator_pull_deflection'   , -10   , -30  ,  30    ,  1.0 ,  1*Units.degree], 
                  [ 'ht_elevator_chord_fraction' , 0.35  , 0.15 , 0.45   ,  1.0  ,  1*Units.less],
                  [ 'ht_elevator_span_frac_start', 0.1   , 0.05 , 0.5    ,  1.0  ,  1*Units.less], 
                  [ 'ht_elevator_span_frac_end'  , 0.9   , 0.6  , 0.95   ,  1.0  ,  1*Units.less],   
                  
    ],dtype=object)   

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'elevator_surface_area', 1E-2 , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([   
        [ 'CL_pull_residual'                  ,   '<' ,   1E-3  ,   1.0  , 1*Units.less],  
        [ 'CL_push_residual'                  ,   '<' ,   1E-3  ,   1.0  , 1*Units.less],   
        [ 'pull_up_CM_residual'               ,   '<' ,   1E-3  ,   1.0  , 1*Units.less],     
        [ 'push_over_CM_residual'             ,   '<' ,   1E-3  ,   1.0  , 1*Units.less],      
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [   
        [ 'pull_up_CM_residual'               , 'summary.pull_up_CM_residual' ],   
        [ 'push_over_CM_residual'             , 'summary.push_over_CM_residual' ],  
        [ 'CL_pull_residual'                  , 'summary.CL_pull_residual' ],  
        [ 'CL_push_residual'                  , 'summary.CL_push_residual' ],    
        [ 'elevator_surface_area'             , 'summary.elevator_surface_area' ], 
        [ 'pull_up_AoA'                       , 'missions.elevator_sizing_pull_up.segments.cruise.angle_of_attack' ], 
        [ 'push_over_AoA'                     , 'missions.elevator_sizing_push_over.segments.cruise.angle_of_attack' ], 
        [ 'elevator_push_deflection'          , 'vehicle_configurations.elevator_sizing_pull_up.wings.horizontal_tail.control_surfaces.elevator.deflection' ],   
        [ 'elevator_pull_deflection'          , 'vehicle_configurations.elevator_sizing_push_over.wings.horizontal_tail.control_surfaces.elevator.deflection' ],     
        [ 'ht_elevator_chord_fraction'        , 'vehicle_configurations.*.wings.horizontal_tail.control_surfaces.elevator.chord_fraction'],    
        [ 'ht_elevator_span_frac_start'       , 'vehicle_configurations.*.wings.horizontal_tail.control_surfaces.elevator.span_fraction_start'],    
        [ 'ht_elevator_span_frac_end'         , 'vehicle_configurations.*.wings.horizontal_tail.control_surfaces.elevator.span_fraction_end'],  
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.elevator_sizing_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.elevator_sizing_setup(nexus.analyses,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.elevator_sizing_procedure()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  
 
 
def aileron_rudder_sizing_optimization_setup(vehicle,cruise_velocity,cruise_altitude):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude    

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    if vehicle.rudder_flag:
        problem.inputs = np.array([         
                      [ 'crosswind_AoA'                ,  5  , -5   ,  20  ,  1.0  ,  1*Units.degree],
                      [ 'aileron_roll_deflection'      ,  0  , -30  , 30   ,  1.0  ,  1*Units.degree],   
                      [ 'aileron_crosswind_deflection' ,  0  , -30  , 30   ,  1.0  ,  1*Units.degree], 
                      [ 'rudder_crosswind_deflection'  ,  0  , -30  , 30   ,  1.0  ,  1*Units.degree],  
                      [ 'mw_aileron_chord_fraction'    , 0.2 , 0.15 , 0.3  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_start'   , 0.7 , 0.55 , 0.8  ,  1.0 ,  1*Units.less],
                      [ 'mw_aileron_span_frac_end'     , 0.9 , 0.85 , 0.95 ,  1.0 ,  1*Units.less],  
                      [ 'vs_rudder_chord_fraction'     , 0.4 , 0.15 , 0.5  ,  1.0 ,  1*Units.less],
                      [ 'vs_rudder_span_frac_start'    , 0.1 , 0.05 , 0.35 ,  1.0 ,  1*Units.less],
                      [ 'vs_rudder_span_frac_end'      , 0.9 , 0.5  , 0.95 ,  1.0 ,  1*Units.less],      
                      ],dtype=object)   
    else:
        problem.inputs = np.array([              
                      [ 'crosswind_AoA'                ,  5  , -5   ,  20   ,  1.0  ,  1*Units.degree], 
                      [ 'aileron_roll_deflection'      ,  0  , -30  , 30    ,  1.0  ,  1*Units.degree],    
                      [ 'aileron_crosswind_deflection' ,  5  , -30   , 30   ,  1.0  ,  1*Units.degree], 
                      [ 'mw_aileron_chord_fraction'    , 0.2  , 0.15 , 0.3  ,  1.0  ,  1*Units.less],
                      [ 'mw_aileron_span_frac_start'   , 0.7  , 0.55 , 0.8  ,  1.0  ,  1*Units.less],
                      [ 'mw_aileron_span_frac_end'     , 0.9  , 0.85 , 0.95 ,  1.0  ,  1*Units.less],   
                      
        ],dtype=object)      

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'aileron_rudder_surface_area', 1E-1 , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ] 
    problem.constraints = np.array([
        [ 'roll_rate_residual'             ,   '<' ,   5E-3 ,   1.0 , 1*Units.less],  # 1E-2 works 
        [ 'F_y_residual'          ,   '<' ,   5E-2 ,   1.0 , 1*Units.less],  # 1E-1 works 
        [ 'F_z_residual'          ,   '<' ,   5E-2 ,   1.0 , 1*Units.less],  # 1E-1 works  
        [ 'M_x_residual'          ,   '<' ,   5E-2 ,   1.0 , 1*Units.less],  # 1E-1 works 
        [ 'M_y_residual'          ,   '<' ,   5E-2 ,   1.0 , 1*Units.less],  # 1E-1 works 
        [ 'M_z_residual'          ,   '<' ,   5E-2 ,   1.0 , 1*Units.less],  # 1E-1 works 
    ],dtype=object)
        
        
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    if vehicle.rudder_flag:
        crosswind_maneuver_segment_name =  'missions.crosswind_maneuver.segments.cruise'  
        problem.aliases = [   
            [ 'crosswind_AoA'                          , crosswind_maneuver_segment_name + '.angle_of_attack' ],   
            [ 'aileron_roll_deflection'                , 'vehicle_configurations.aileron_rudder_roll_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],  
            [ 'aileron_crosswind_deflection'           , 'vehicle_configurations.aileron_rudder_crosswind_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],   
            [ 'rudder_crosswind_deflection'            , 'vehicle_configurations.aileron_rudder_crosswind_sizing.wings.vertical_tail.control_surfaces.rudder.deflection' ],   
            [ 'mw_aileron_chord_fraction'              , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.chord_fraction'],  
            [ 'mw_aileron_span_frac_start'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_start'],    
            [ 'mw_aileron_span_frac_end'               , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_end'],     
            [ 'vs_rudder_chord_fraction'               , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.chord_fraction'],    
            [ 'vs_rudder_span_frac_start'              , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.span_fraction_start'],    
            [ 'vs_rudder_span_frac_end'                , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.span_fraction_end'], 
            [ 'aileron_rudder_surface_area'            , 'summary.aileron_rudder_surface_area' ], 
            [ 'roll_rate_residual'                     , 'summary.roll_rate_residual' ],  
            [ 'F_y_residual'                  , 'summary.F_y_residual' ], 
            [ 'F_z_residual'                  , 'summary.F_z_residual' ],
            [ 'M_x_residual'                  , 'summary.M_x_residual' ], 
            [ 'M_y_residual'                  , 'summary.M_y_residual' ], 
            [ 'M_z_residual'                  , 'summary.M_z_residual' ], 
        ]      
    else:
        problem.aliases = [  
            [ 'F_y_residual'          , 'summary.F_y_residual' ], 
            [ 'F_z_residual'          , 'summary.F_z_residual' ],
            [ 'M_x_residual'          , 'summary.M_x_residual' ], 
            [ 'M_y_residual'          , 'summary.M_y_residual' ], 
            [ 'M_z_residual'          , 'summary.M_z_residual' ],  
            [ 'roll_rate_residual'             , 'summary.roll_rate_residual' ], 
            [ 'aileron_rudder_surface_area'    , 'summary.aileron_rudder_surface_area' ],  
            [ 'aileron_roll_deflection'        , 'vehicle_configurations.aileron_rudder_roll_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],          
            [ 'aileron_crosswind_deflection'   , 'vehicle_configurations.aileron_rudder_crosswind_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],   
            [ 'mw_aileron_chord_fraction'      , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.chord_fraction'],  
            [ 'mw_aileron_span_frac_start'     , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_start'],    
            [ 'mw_aileron_span_frac_end'       , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_end']] 
        
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.aileron_rudder_sizing_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # ------------------------------------------------------------------- 
    nexus.missions = Missions.aileron_rudder_sizing_setup(nexus.analyses,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.aileron_rudder_sizing_procedure()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  

def aileron_rudder_oei_sizing_optimization_setup(vehicle,cruise_velocity,cruise_altitude):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude    

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    if vehicle.rudder_flag:
        problem.inputs = np.array([          
                      [ 'aileron_oei_deflection'       ,  5  , -30  ,  30   ,  1.0  ,  1*Units.degree], 
                      [ 'rudder_oei_deflection'        ,  5  , -30  ,  30   ,  1.0  ,  1*Units.degree],
                      [ 'elevator_oei_deflection'      ,  5  , -30  ,  30   ,  1.0  ,  1*Units.degree],
                      [ 'OEI_AoA'                      ,  5  , -5   ,  20  ,  1.0  ,  1*Units.degree],
                      [ 'PR_OEI_eta_1'                 , 0.5 , 0    ,  1.0 ,  1.0  ,  1*Units.less],        
                      ],dtype=object)   
    else:
        problem.inputs = np.array([              
                      [ 'OEI_AoA'                      ,  5  , -5   ,  20   ,  1.0  ,  1*Units.degree], 
                      [ 'aileron_oei_deflection'       ,  0  , -30  , 30    ,  1.0  ,  1*Units.degree],     
                      [ 'mw_aileron_chord_fraction'    , 0.2  , 0.15 , 0.3  ,  1.0  ,  1*Units.less],
                      [ 'mw_aileron_span_frac_start'   , 0.7  , 0.55 , 0.8  ,  1.0  ,  1*Units.less],
                      [ 'mw_aileron_span_frac_end'     , 0.9  , 0.85 , 0.95 ,  1.0  ,  1*Units.less],   
                      
        ],dtype=object)      

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                  [ 'oei_objective', 1E-1 , 1*Units.kg],
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ] 
    problem.constraints = np.array([ 
        [ 'PR_OEI_CX_residual'             ,   '<' ,   1E-2 ,   1 , 1*Units.less],
        [ 'PR_OEI_CY_residual'             ,   '<' ,   1E-2 ,   1 , 1*Units.less],  
        [ 'PR_OEI_CZ_residual'             ,   '<' ,   1E-2 ,   1 , 1*Units.less],     
        [ 'PR_OEI_CL_residual'             ,   '<' ,   1E-2 ,   1 , 1*Units.less],   
        [ 'PR_OEI_CM_residual'             ,   '<' ,   1E-2 ,   1 , 1*Units.less],     
        [ 'PR_OEI_CN_residual'             ,   '<' ,   1E-2 ,   1 , 1*Units.less],  
    ],dtype=object)
        
        
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    if vehicle.rudder_flag: 
        cruis_oei_maneuver_segment_name =  'missions.cruise_oei.segments.cruise' 
        problem.aliases = [    
            [ 'aileron_oei_deflection'                 , 'vehicle_configurations.aileron_rudder_oei_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],   
            [ 'rudder_oei_deflection'                  , 'vehicle_configurations.aileron_rudder_oei_sizing.wings.vertical_tail.control_surfaces.rudder.deflection' ], 
            [ 'elevator_oei_deflection'                , 'vehicle_configurations.aileron_rudder_oei_sizing.wings.horizontal_tail.control_surfaces.elevator.deflection' ], 
            [ 'oei_objective'                          , 'summary.oei_objective' ],        
            [ 'OEI_AoA'                                ,  cruis_oei_maneuver_segment_name + '.angle_of_attack' ], 
            [ 'PR_OEI_eta_1'                           ,  cruis_oei_maneuver_segment_name + '.assigned_control_variables.throttle.initial_guess_values[0][0]'],    
            [ 'PR_OEI_CX_residual'                     , 'summary.F_x_residual'],  
            [ 'PR_OEI_CY_residual'                     , 'summary.F_x_residual'],  
            [ 'PR_OEI_CZ_residual'                     , 'summary.F_z_residual'],     
            [ 'PR_OEI_CL_residual'                     , 'summary.M_x_residual'],   
            [ 'PR_OEI_CM_residual'                     , 'summary.M_y_residual'],     
            [ 'PR_OEI_CN_residual'                     , 'summary.M_z_residual'],  
        ]      
    else:
        problem.aliases = [  
            [ 'PR_OEI_CX_residual'                     , 'summary.F_x_residual'],  
            [ 'PR_OEI_CY_residual'                     , 'summary.F_x_residual'],  
            [ 'PR_OEI_CZ_residual'                     , 'summary.F_z_residual'],     
            [ 'PR_OEI_CL_residual'                     , 'summary.M_x_residual'],   
            [ 'PR_OEI_CM_residual'                     , 'summary.M_y_residual'],     
            [ 'PR_OEI_CN_residual'                     , 'summary.M_z_residual'],   
            [ 'aileron_rudder_surface_area'            , 'summary.aileron_rudder_surface_area' ],  
            [ 'aileron_oei_deflection'                 , 'vehicle_configurations.aileron_rudder_oei_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],      
            [ 'mw_aileron_chord_fraction'              , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.chord_fraction'],  
            [ 'mw_aileron_span_frac_start'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_start'],    
            [ 'mw_aileron_span_frac_end'               , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_end']] 
        
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.aileron_rudder_oei_sizing_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # ------------------------------------------------------------------- 
    nexus.missions = Missions.aileron_rudder_oei_sizing_setup(nexus.analyses,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.aileron_rudder_oei_sizing_procedure()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus

def flap_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude):

    nexus = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude    

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,         (lb , ub)        , scaling , units ]  
    problem.inputs = np.array([           
                  [ 'mw_flap_chord_fraction'     , 0.2    , 0.15 , 0.4  ,  1.0 ,  1*Units.less],
                  [ 'mw_flap_span_frac_start'    , 0.2    , 0.05 , 0.25 ,  1.0 ,  1*Units.less],
                  [ 'mw_flap_span_frac_end'      , 0.5    , 0.3  , 0.6  ,  1.0 ,  1*Units.less],   
                  
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
    nexus.vehicle_configurations = Vehicles.flap_sizing_setup(optimized_vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.flap_sizing_setup(nexus.analyses,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.flap_sizing_procedure()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  

def hover_oei_optimization_setup(vehicle,cruise_velocity,cruise_altitude): 
    nexus   = Nexus()
    problem = Data()
    nexus.optimization_problem = problem
    
    nexus.cruise_velocity = cruise_velocity
    nexus.cruise_altitude = cruise_altitude     

    # -------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------- 
    problem.inputs = np.array([       

    #                [ tag             , initial  ,    (lb  ,    ub)  , scaling ,      units ]  
                  [ 'PR_OEI_eta_group_1'  , 0.5     , 0       ,  1.0     ,  1.0   ,  1*Units.less],   
                  [ 'PR_OEI_eta_group_2'  , 0.5     , 0       ,  1.0     ,  1.0   ,  1*Units.less],    
                  [ 'PR_OEI_eta_group_3'  , 0.5     , 0       ,  1.0     ,  1.0   ,  1*Units.less],    
                  [ 'PR_OEI_eta_group_4'  , 0.5     , 0       ,  1.0     ,  1.0   ,  1*Units.less], 
                  [ 'LR_OEI_eta_group_1'  , 0.5     , 0       ,  1.0     ,  1.0   ,  1*Units.less],   
                  [ 'LR_OEI_eta_group_2'  , 0.5     , 0       ,  1.0     ,  1.0   ,  1*Units.less],    
                  [ 'LR_OEI_eta_group_3'  , 0.5     , 0       ,  1.0     ,  1.0   ,  1*Units.less],    
                  [ 'LR_OEI_eta_group_4'  , 0.5     , 0       ,  1.0     ,  1.0   ,  1*Units.less],                         
                  
    ],dtype=object)       

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([ 
                                 [  'Power'  ,  1   ,    1*Units.less] 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ] See http://everyspec.com/MIL-SPECS/MIL-SPECS-MIL-F/MIL-F-8785C_5295/
    # Level 1, Category B
    problem.constraints = np.array([
        [ 'PR_OEI_F_z_residual'      ,   '<' ,   1E-3  ,   1E-3  , 1*Units.less], 
        [ 'PR_OEI_M_x_residual'      ,   '<' ,   1E-3  ,   1E-3  , 1*Units.less],   
        [ 'PR_OEI_M_y_residual'      ,   '<' ,   1E-3  ,   1E-3  , 1*Units.less],   
        [ 'LR_OEI_F_z_residual'      ,   '<' ,   1E-3  ,   1E-3  , 1*Units.less], 
        [ 'LR_OEI_M_x_residual'      ,   '<' ,   1E-3  ,   1E-3  , 1*Units.less],   
        [ 'LR_OEI_M_y_residual'      ,   '<' ,   1E-3  ,   1E-3  , 1*Units.less],  
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ] 
    problem.aliases = [ 
        [ 'PR_OEI_F_z_residual'     , 'summary.PR_OEI_F_z_residual' ],  
        [ 'PR_OEI_M_x_residual'     , 'summary.PR_OEI_M_x_residual' ], 
        [ 'PR_OEI_M_y_residual'     , 'summary.PR_OEI_M_y_residual' ], 
        [ 'PR_OEI_M_z_residual'     , 'summary.PR_OEI_M_z_residual' ],
        [ 'LR_OEI_F_z_residual'     , 'summary.LR_OEI_F_z_residual' ],  
        [ 'LR_OEI_M_x_residual'     , 'summary.LR_OEI_M_x_residual' ], 
        [ 'LR_OEI_M_y_residual'     , 'summary.LR_OEI_M_y_residual' ],            
        [ 'Power'                   , 'summary.Power'],  
        [ 'PR_OEI_eta_group_1'      , ['missions.crosswind_maneuver.segments.cruise.conditions.energy.prop_rotor_bus.[propulsor.tag].throttle']],   
        [ 'PR_OEI_eta_group_2'      , ''],    
        [ 'PR_OEI_eta_group_3'      , ''],    
        [ 'PR_OEI_eta_group_4'      , ''],         
        [ 'LR_OEI_eta_group_1'      , 'missions.crosswind_maneuver.segments.cruise.conditions.energy.cruise_bus.[propulsor.tag].throttle'],   
        [ 'LR_OEI_eta_group_2'      , ''],    
        [ 'LR_OEI_eta_group_3'      , ''],    
        [ 'LR_OEI_eta_group_4'      , ''],     
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.hover_oei_setup(vehicle)
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.analyses_setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.hover_oei_setup(nexus.analyses,nexus.cruise_velocity,nexus.cruise_altitude)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.hover_oei_procedure()
    
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
 
def save_aircraft_geometry(geometry,filename): 
    pickle_file  = filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(geometry, file) 
    return 


def load_aircraft_geometry(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results


if __name__ == '__main__':
    main()
    plt.show()
