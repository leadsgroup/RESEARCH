# vlm_pertubation_test.py
# 
# Created: May 2024, M. Clarke
 
# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
import RCAIDE 
from   RCAIDE.Framework.Core import Units     
from   RCAIDE.Library.Plots  import *  

from   Optimize import size_control_surfaces

# python imports  
import pylab as plt 
import numpy as np
from   copy import deepcopy

# local imports 
import sys 
import os

sys.path.append(os.path.join(os.path.split(os.path.split(os.path.split(sys.path[0])[0])[0])[0], 'Aircraft'))

from   Stopped_Rotor.Stopped_Rotor                                        import vehicle_setup as SR_vehicle_setup   
from   Tiltrotor.Tiltrotor                                                import vehicle_setup as TR_vehicle_setup   
from   Tiltwing.Tiltwing                                                  import vehicle_setup as TW_vehicle_setup   
from   Hexacopter.Hexacopter                                              import vehicle_setup as HC_vehicle_setup   
from   Tilt_Stopped_Rotor.Tilt_Stopped_Rotor_Conv_Tail                    import vehicle_setup as TSR_vehicle_setup   
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():  
    
    aircraft_model  = 'TSR'
    cruise_velocity = 150  * Units['mph']
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
    #                 CG: X,    Y,    Z,  wing_flag 
    CG_bat_1 = np.array([[0.25, 0.,   0., True], # Change wing flag back to false
                         [0.35, 0.,   0., False],
                         [0.45, 0.,   0., False],
                         [1.9375, 0.25, 1.0, True],
                         [1.9375, 0.35, 1.0, True],
                         [1.9375, 0.45, 1.0, True]])
    
    # lift rotor battery modules
    #                 CG: X,    Y,    Z 
    CG_bat_2 = np.array([[4.0,  0.,   0.],
                         [4.1,  0.,   0.],
                         [4.2,  0.,   0.]])   
 
    for i in range(len(CG_bat_1)-3):
        for j in range(len(CG_bat_2)):
            
            
            # prop rotor battery modules      
            case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_1[i, 0:3]])
            case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_1[i,0] + case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.length, 
                                                                                                                 -CG_bat_1[i,1], 
                                                                                                                 CG_bat_1[i,2]])
         
            # If the prop rotor batteries are in the wing:
            if CG_bat_1[i, 3]:
                # rotate Battery 90 degrees
                length1 = case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.length
                width1 = case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.width
                
                case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.length = width1
                case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.width = length1       
                
                length2 = case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.length
                width2 = case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.width
                
                case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.length = width2
                case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.width = length2
                
                case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_1[i, 0],
                                                                                                                     CG_bat_1[i,1] + width2 / 2, 
                                                                                                                     CG_bat_1[i,2]                                                                                                                     ])
                case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_1[i,0], 
                                                                                                                     -CG_bat_1[i,1] - width2 / 2, 
                                                                                                                     CG_bat_1[i,2]])
                
            #case_vehicle.networks.electric.busses.prop_rotor_bus.initialize_bus_properties()    
            # lift rotor battery modules 

            case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_2[j,  0:3]])
            case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_2[j,0], 
                                                                                                                 CG_bat_2[j,1], 
                                                                                                                 CG_bat_2[j,2] + case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.height])
            size_control_surfaces(CG_bat_1[i], CG_bat_2[j], case_vehicle, cruise_velocity, cruise_altitude)
        
    for i in range(len(CG_bat_1)-3):
        for j in range(len(CG_bat_2)):
            
            
            # prop rotor battery modules      
            case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_1[i, 0:3]])
            case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_1[i,0] + case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.length, 
                                                                                                                 -CG_bat_1[i,1], 
                                                                                                                 CG_bat_1[i,2]])
            # If the prop rotor batteries are in the wing:
            if CG_bat_1[i, 3]:
                # rotate Battery 90 degrees
                case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_1.orientation_euler_angles = [-np.pi / 2, 0, 0]
                case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.nmc_module_2.orientation_euler_angles = [np.pi / 2, 0, 0]
                
            # lift rotor battery modules 
            case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_1.origin = np.array([CG_bat_2[j,  0:3]])
            case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.origin = np.array([CG_bat_2[i,0], 
                                                                                                                 CG_bat_2[i,1], 
                                                                                                                 CG_bat_2[i,2] + case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.nmc_module_2.height])
            size_control_surfaces(CG_bat_1[i], CG_bat_2[j], case_vehicle, cruise_velocity, cruise_altitude)    
    return 

if __name__ == '__main__': 
    main()    
    plt.show() 