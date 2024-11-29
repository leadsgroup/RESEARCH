# Optimize.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------  
import RCAIDE
from   RCAIDE.Framework.Core                        import Units, Data 
from   RCAIDE.Framework.Optimization.Packages.scipy import scipy_setup
from   RCAIDE.Framework.Optimization.Common         import Nexus
from   save_data_stick_fixed          import save_data_stick_fixed 
from   RCAIDE import  load 
from   RCAIDE import  save

# python imports 
import numpy as np
import time
import pickle
from   copy  import  deepcopy

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
    
    '''
    STICK FIXED (STATIC STABILITY AND DRAG OTIMIZATION
    '''
    ti                            = time.time()   
    planform_optimization_problem = stick_fixed_stability_and_drag_optimization_setup(vehicle,cruise_velocity,cruise_altitude)
    output_stick_fixed            = scipy_setup.SciPy_Solve(planform_optimization_problem,solver='SLSQP', sense_step = 1E-3, tolerance = 1E-3)   
    print (output_stick_fixed)    
    tf                            = time.time()
    elapsed_time_stick_fixed = round((tf-ti)/60,2)
    print('Stick Fixed Stability and Drag Otimization Simulation Time: ' + str(elapsed_time_stick_fixed))
    optimized_vehicle                             = planform_optimization_problem.vehicle_configurations.stick_fixed_cruise   
    
    save_data_stick_fixed(CG_bat_1, 
                          CG_bat_2,
                          optimized_vehicle.tag, 
                          output_stick_fixed, 
                          planform_optimization_problem,
                          elapsed_time_stick_fixed
                         )    
    
    cg_x1   =  str(CG_bat_1[0]).replace('.', "")
    cg_y1   =  str(CG_bat_1[1]).replace('.', "")
    cg_z1   =  str(CG_bat_1[2]).replace('.', "")
    cg_x2   =  str(CG_bat_2[0]).replace('.', "")
    cg_y2   =  str(CG_bat_2[1]).replace('.', "")
    cg_z2   =  str(CG_bat_2[2]).replace('.', "")
    
    def new_filename(filename):
        return filename.replace('[', '').replace(']', '').replace(' ', '_').replace('.', '_')

    raw_file_name1 = f"{cg_x1}_{cg_y1}_{cg_z1}_{cg_x2}_{cg_y2}_{cg_z2}_Optimized_Vehicle"
    vehicle_file_name = new_filename(raw_file_name1)
    save_data(optimized_vehicle,vehicle_file_name)
    
    raw_file_name2 = f"{cg_x1}_{cg_y1}_{cg_z1}_{cg_x2}_{cg_y2}_{cg_z2}_Planform_Optimization_Problem"
    planform_optimization_problem_file_name = new_filename(raw_file_name2)
    save_data(planform_optimization_problem,planform_optimization_problem_file_name)  
    
    raw_file_name3 = f"{cg_x1}_{cg_y1}_{cg_z1}_{cg_x2}_{cg_y2}_{cg_z2}_Output_Stick_Fixed"
    output_stick_fixed_file_name = new_filename(raw_file_name3)
    save_data(output_stick_fixed,output_stick_fixed_file_name)      
    
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
 
def save_data(geometry,filename): 
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
