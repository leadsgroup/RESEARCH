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
