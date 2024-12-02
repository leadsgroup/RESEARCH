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
import sys 
import os

# local imports 
import Vehicles
import Missions
import Analyses
import Procedure
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.split(os.path.split(os.path.split(sys.path[0])[0])[0])[0], 'Aircraft'))
from append_control_surfaces             import append_control_surfaces
from compute_control_surface_derivatives import compute_control_surface_derivatives 
from Stopped_Rotor.Stopped_Rotor                                        import vehicle_setup as SR_vehicle_setup   
from Tiltrotor.Tiltrotor                                                import vehicle_setup as TR_vehicle_setup   
from Tiltwing.Tiltwing                                                  import vehicle_setup as TW_vehicle_setup   
from Hexacopter.Hexacopter                                              import vehicle_setup as HC_vehicle_setup   
from Tilt_Stopped_Rotor.Tilt_Stopped_Rotor_Conv_Tail                    import vehicle_setup as TSR_vehicle_setup   
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------        
#   Run the whole thing
# ----------------------------------------------------------------------  
def main():
    
 
    aircraft_model  =  'TSR'
    cruise_velocity = 120  * Units['mph']
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
        
    case_vehicle  = deepcopy(vehicle)
    
    # delete control surfaces if they have been defined 
    for wing in case_vehicle.wings:
        for control_surface in wing.control_surfaces:
            del wing.control_surfaces[control_surface.tag]
    
    
    # prop rotor battery module (first module)
    #                 CG: X,    Y,  Z 
    CG_bat_1 = np.array([[0.25, 0., 0.],
                         [0.35, 0., 0.],
                         [0.45, 0., 0.]])
    
    # lift rotor battery modules
    #                 CG: X,    Y,  Z 
    CG_bat_2 = np.array([[4.0,  0., 0.],
                         [4.1,  0., 0.],
                         [4.2,  0., 0.]])   
 
    for i in range(len(CG_bat_1)):
        for j in range(len(CG_bat_2)):
            
            
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
       
        
            '''
            STICK FIXED (STATIC STABILITY AND DRAG OTIMIZATION
            '''
            ti   = time.time()   
            planform_optimization_problem = stick_fixed_stability_and_drag_optimization_setup(case_vehicle,cruise_velocity,cruise_altitude)
            output_stick_fixed = scipy_setup.SciPy_Solve(planform_optimization_problem,solver='SLSQP', sense_step = 1E-3, tolerance = 1E-3)   
            print (output_stick_fixed)    
            tf           = time.time()
            elapsed_time_stick_fixed = round((tf-ti)/60,2)
            print('Stick Fixed Stability and Drag Otimization Simulation Time: ' + str(elapsed_time_stick_fixed))
            
            
            '''
            RUN LOOP TO GET DERIVATIVE FUNCTIONS OF CONTROL SURFACES     
            '''
            vehicle =  append_control_surfaces(vehicle)
             
            
            '''
            RUN LOOP TO GET DERIVATIVE FUNCTIONS OF CONTROL SURFACES     
            '''
            derivatives = compute_control_surface_derivatives(vehicle)
            
            
            '''
            ELEVATOR SIZING (7 mins)
            '''
            ti = time.time()    
            optimized_vehicle                             = planform_optimization_problem.vehicle_configurations.stick_fixed_cruise 
            save(optimized_vehicle, 'optimized_vehicle_stick_fixed', pickle_format=True)  
            elevator_sizing_optimization_problem = elevator_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude,derivatives)
            output_elevator_sizing = scipy_setup.SciPy_Solve(elevator_sizing_optimization_problem,solver='SLSQP', sense_step = 1E-2, tolerance = 1E-2) 
            print (output_elevator_sizing)     
            tf           = time.time()
            elapsed_time_elevator_sizing = round((tf-ti)/60,2)
            print('Elevator Sizing Simulation Time: ' + str(elapsed_time_elevator_sizing))
           
            
            '''
            AILERON AND RUDDER SIZING (22 mins)
            '''       
            ti = time.time()  
            optimized_vehicle                    = elevator_sizing_optimization_problem.vehicle_configurations.elevator_sizing 
            save(optimized_vehicle, 'optimized_vehicle_ele', pickle_format=True)   
            optimized_vehicle.rudder_flag        = True  
            optimized_vehicle.crosswind_velocity = 20 * Units.knots 
            trim_angle_of_attack =  0.09328504 # output_elevator_sizing[0]
            aileron_rudder_sizing_optimization_problem = aileron_rudder_sizing_optimization_setup(optimized_vehicle,cruise_velocity,cruise_altitude, trim_angle_of_attack,derivatives)
            output_aileron_and_rudder_sizing = scipy_setup.SciPy_Solve(aileron_rudder_sizing_optimization_problem,solver='SLSQP', sense_step = 1E-3, tolerance = 1E-3) 
            print (output_aileron_and_rudder_sizing)     
            tf           = time.time()
            elapsed_time_aileron_and_rudder_sizing = round((tf-ti)/60,2)
            print('Aileron and Rudder Sizing Simulation Time: ' + str(elapsed_time_aileron_and_rudder_sizing))   
            
              
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
        [ 'pull_up_AoA'                       , 'summary.elevator_pull_up_angle_of_attack' ], 
        [ 'push_over_AoA'                     , 'summary.elevator_push_over_angle_of_attack' ], 
        [ 'elevator_push_deflection'          , 'summary.elevator_push_over_deflection' ],   
        [ 'elevator_pull_deflection'          , 'summary.elevator_pull_up_deflection' ],      
        [ 'ht_elevator_span_frac_start'       , 'summary.elevator_span_fraction_start'],    
        [ 'ht_elevator_span_frac_end'         , 'summary.elevator_span_fraction_end'],  
    ]      
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = None
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = None
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = None
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.elevator_sizing_procedure()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()     
    return nexus  
 
 
def aileron_rudder_sizing_optimization_setup(vehicle,cruise_velocity,cruise_altitude,trim_angle_of_attack):

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
                  [ 'aileron_roll_deflection'      ,  0  , -30  , 30   ,  1.0  ,  1*Units.degree],   
                  [ 'aileron_crosswind_deflection' ,  0  , -30  , 30   ,  1.0  ,  1*Units.degree], 
                  [ 'rudder_crosswind_deflection'  ,  0  , -30  , 30   ,  1.0  ,  1*Units.degree],  
                  [ 'PR_OEI_eta_1'                 , 0.5 ,  0.1 ,  1.0 ,  1.0  ,  1*Units.less], 
                  [ 'aileron_oei_deflection'       ,  0  , -30  ,  30   ,  1.0 , 1*Units.degree],  
                  [ 'rudder_oei_deflection'        ,  0  , -30  ,  30  ,  1.0 ,  1*Units.degree], 
                  [ 'mw_aileron_chord_fraction'    , 0.2 , 0.15 , 0.3  ,  1.0 ,  1*Units.less],
                  [ 'mw_aileron_span_frac_start'   , 0.7 , 0.55 , 0.8  ,  1.0 ,  1*Units.less],
                  [ 'mw_aileron_span_frac_end'     , 0.9 , 0.85 , 0.95 ,  1.0 ,  1*Units.less],  
                  [ 'vs_rudder_chord_fraction'     , 0.4 , 0.15 , 0.5  ,  1.0 ,  1*Units.less],
                  [ 'vs_rudder_span_frac_start'    , 0.1 , 0.05 , 0.35 ,  1.0 ,  1*Units.less],
                  [ 'vs_rudder_span_frac_end'      , 0.9 , 0.5  , 0.95 ,  1.0 ,  1*Units.less],      
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
        [ 'roll_rate_residual'    ,   '<' ,   1E-6 ,   1 , 1*Units.less],  # 1E-2 works 
        [ 'F_y_residual'          ,   '<' ,   1E-6 ,   1 , 1*Units.less],  # 1E-1 works   
        [ 'M_y_residual'          ,   '<' ,   1E-6 ,   1 , 1*Units.less],  # 1E-1 works 
        [ 'M_z_residual'          ,   '<' ,   1E-6 ,   1 , 1*Units.less],  # 1E-1 works 
        [ 'OEI_CY_residual'       ,   '<' ,   1E-6 ,   1 , 1*Units.less],    
        [ 'OEI_CL_residual'       ,   '<' ,   1E-6 ,   1 , 1*Units.less],    
        [ 'OEI_CN_residual'       ,   '<' ,   1E-6 ,   1 , 1*Units.less],  
    ],dtype=object)
        
        
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ]   
    cruise_oei_maneuver_segment_name =  'missions.cruise_oei.segments.cruise' 
    problem.aliases = [    
        [ 'OEI_eta_1'                              , cruise_oei_maneuver_segment_name + '.assigned_control_variables.throttle.initial_guess_values[0][0]'],  
        [ 'aileron_roll_deflection'                , 'vehicle_configurations.aileron_rudder_roll_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],  
        [ 'aileron_crosswind_deflection'           , 'vehicle_configurations.aileron_rudder_crosswind_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],   
        [ 'rudder_crosswind_deflection'            , 'vehicle_configurations.aileron_rudder_crosswind_sizing.wings.vertical_tail.control_surfaces.rudder.deflection' ], 
        [ 'aileron_oei_deflection'                 , 'vehicle_configurations.aileron_rudder_oei_sizing.wings.main_wing.control_surfaces.aileron.deflection' ],      
        [ 'rudder_oei_deflection'                  , 'vehicle_configurations.aileron_rudder_oei_sizing.wings.vertical_tail.control_surfaces.rudder.deflection' ],     
        [ 'mw_aileron_chord_fraction'              , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.chord_fraction'],  
        [ 'mw_aileron_span_frac_start'             , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_start'],    
        [ 'mw_aileron_span_frac_end'               , 'vehicle_configurations.*.wings.main_wing.control_surfaces.aileron.span_fraction_end'],     
        [ 'vs_rudder_chord_fraction'               , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.chord_fraction'],    
        [ 'vs_rudder_span_frac_start'              , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.span_fraction_start'],    
        [ 'vs_rudder_span_frac_end'                , 'vehicle_configurations.*.wings.vertical_tail.control_surfaces.rudder.span_fraction_end'], 
        [ 'aileron_rudder_surface_area'            , 'summary.aileron_rudder_surface_area' ], 
        [ 'roll_rate_residual'                     , 'summary.roll_rate_residual' ],  
        [ 'F_y_residual'                           , 'summary.F_y_residual' ],  
        [ 'M_y_residual'                           , 'summary.M_y_residual' ], 
        [ 'M_z_residual'                           , 'summary.M_z_residual' ],   
        [ 'OEI_CY_residual'                        , 'summary.OEI_F_y_residual'],    
        [ 'OEI_CL_residual'                        , 'summary.OEI_M_x_residual'],       
        [ 'OEI_CN_residual'                        , 'summary.OEI_M_z_residual'],             
        ]       
        
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
    nexus.missions = Missions.aileron_rudder_sizing_setup(nexus.analyses,nexus.cruise_velocity,nexus.cruise_altitude, angle_of_attack=trim_angle_of_attack)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.aileron_rudder_sizing_procedure()
    
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
