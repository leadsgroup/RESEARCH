# vlm_pertubation_test.py
# 
# Created: May 2024, M. Clarke
 
# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
import RCAIDE 
from RCAIDE.Framework.Core import Units     
from RCAIDE.Library.Plots  import *  

# python imports  
import pylab as plt 


# local imports 
import sys 
import os

sys.path.append(os.path.join( os.path.split(os.path.split(sys.path[0])[0])[0], 'Aircraft'))
from Stopped_Rotor         import vehicle_setup as SR_vehicle_setup   
from Tiltrotor             import vehicle_setup as TR_vehicle_setup   
from Tiltwing              import vehicle_setup as TW_vehicle_setup   
from Hexacopter            import vehicle_setup as HC_vehicle_setup   
from Tilt_Stopped_Rotor    import vehicle_setup as TSR_vehicle_setup   
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():  
    aircraft_model =  'SR' 
    
    if aircraft_model == 'SR':
        vehicle =  SR_vehicle_setup()
    if aircraft_model == 'TR':
        vehicle =  TR_vehicle_setup()
    if aircraft_model == 'SR':
        vehicle =  TW_vehicle_setup()
    if aircraft_model == 'SR':
        vehicle =  HC_vehicle_setup()
    if aircraft_model == 'SR':
        vehicle =  TSR_vehicle_setup()
        
    size_control_surfaces(vehicle)
    
    return 

if __name__ == '__main__': 
    main()    
    plt.show() 