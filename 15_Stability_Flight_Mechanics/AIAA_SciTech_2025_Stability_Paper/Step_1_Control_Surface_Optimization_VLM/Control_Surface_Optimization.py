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
    
    #                                     CG: X,    Y,  Z,  rx  ry  rz 
    configuration_CGbatt_MOIbatt = np.array([[2.6,  0., 0., 0.,   0., 0.],
                                             [2.7,  0., 0., 0.,   0., 0.],
                                             [2.8,  0., 0., 0.,   0., 0.],
                                             [2.9,  0., 0., 0.,   0., 0.],
                                             [3.0,  0., 0., 0.,   0., 0.],
                                             [2.6,  0., 0., 0.25, 0., 0.],
                                             [2.7,  0., 0., 0.25, 0., 0.],
                                             [2.8,  0., 0., 0.25, 0., 0.],
                                             [2.9,  0., 0., 0.25, 0., 0.],
                                             [3.0,  0., 0., 0.25, 0., 0.],
                                             [2.6,  0., 0., 0.5,  0., 0.],
                                             [2.7,  0., 0., 0.5,  0., 0.],
                                             [2.8,  0., 0., 0.5,  0., 0.],
                                             [2.9,  0., 0., 0.5,  0., 0.],
                                             [3.0,  0., 0., 0.5,  0., 0.],
                                             [2.88, 0., 0., 0.,   1., 0.],
                                             [2.88, 0., 0., 0.,   2., 0.],
                                             [2.88, 0., 0., 0.,   3., 0.]])
               
    for i in range(len(configuration_CGbatt_MOIbatt)):
        
        case_vehicle  = deepcopy(vehicle)
      
        case_vehicle.networks.electric.busses.prop_rotor_bus.battery_modules.bus_battery.origin = np.array([[configuration_CGbatt_MOIbatt[i,0]+configuration_CGbatt_MOIbatt[i,3], 
                                                                                                             configuration_CGbatt_MOIbatt[i,1]+configuration_CGbatt_MOIbatt[i,4], 
                                                                                                             configuration_CGbatt_MOIbatt[i,2]+configuration_CGbatt_MOIbatt[i,5]]])
        case_vehicle.networks.electric.busses.lift_rotor_bus.battery_modules.lift_bus_battery.origin = np.array([[configuration_CGbatt_MOIbatt[i,0]-configuration_CGbatt_MOIbatt[i,3], 
                                                                                                             configuration_CGbatt_MOIbatt[i,1]-configuration_CGbatt_MOIbatt[i,4], 
                                                                                                             configuration_CGbatt_MOIbatt[i,2]-configuration_CGbatt_MOIbatt[i,5]]])
                   
        size_control_surfaces(configuration_CGbatt_MOIbatt, case_vehicle)
    
    return 

if __name__ == '__main__': 
    main()    
    plt.show() 