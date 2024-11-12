# vlm_pertubation_test.py
# 
# Created: May 2024, M. Clarke
 
# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
import RCAIDE 
from   RCAIDE.Framework.Core import Units     
from   RCAIDE.Library.Plots  import *  

from Optimize import size_control_surfaces

# python imports  
import pylab as plt 

# local imports 
import sys 
import os

sys.path.append(os.path.join(os.path.split(os.path.split(os.path.split(sys.path[0])[0])[0])[0], 'Aircraft'))

from   Stopped_Rotor.Stopped_Rotor                     import vehicle_setup as SR_vehicle_setup   
from   Tiltrotor.Tiltrotor                             import vehicle_setup as TR_vehicle_setup   
from   Tiltwing.Tiltwing                               import vehicle_setup as TW_vehicle_setup   
from   Hexacopter.Hexacopter                           import vehicle_setup as HC_vehicle_setup   
from   Tilt_Stopped_Rotor.Tilt_Stopped_Rotor_Conv_Tail import vehicle_setup as TSR_vehicle_setup   
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():  
    aircraft_model =  'TSR' 
    
    if aircraft_model == 'SR':
        vehicle =  SR_vehicle_setup()
    if aircraft_model == 'TR':
        vehicle =  TR_vehicle_setup()
    if aircraft_model == 'TW':
        vehicle =  TW_vehicle_setup()
    if aircraft_model == 'HC':
        vehicle =  HC_vehicle_setup()
    if aircraft_model == 'TSR':
        vehicle =  TSR_vehicle_setup()
        
    for i in range(18):
        
        
        size_control_surfaces(vehicle)
    
    return 

if __name__ == '__main__': 
    main()    
    plt.show() 