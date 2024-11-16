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
    
    aircraft_model =  'TSR' 
    
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
            size_control_surfaces(CG_bat_1, CG_bat_2, case_vehicle)
        
    return 

if __name__ == '__main__': 
    main()    
    plt.show() 