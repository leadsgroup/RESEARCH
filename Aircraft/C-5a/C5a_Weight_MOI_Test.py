# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units      
from RCAIDE.Library.Methods.Propulsors.Turbofan_Propulsor   import design_turbofan
from RCAIDE.Library.Methods.Geometry.Planform               import wing_planform, segment_properties
from RCAIDE.Library.Plots                 import *     
import Lockheed_C5a
from RCAIDE.Library.Methods.Weights.Correlation_Buildups import Common           as Common
from RCAIDE.Library.Methods.Stability.Moment_of_Inertia.calculate_aircraft_MOI import calculate_aircraft_MOI
from RCAIDE.Library.Methods.Stability.Center_of_Gravity     import compute_vehicle_center_of_gravity


# python imports 
import numpy as np  
from copy import deepcopy
import matplotlib.pyplot as plt  
import os

def main():
    vehicle = Lockheed_C5a.vehicle_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown
    # ------------------------------------------------------------------ 
    settings = None    
    weight = Common.compute_operating_empty_weight(vehicle, settings = settings, method_type = 'RCAIDE')
    print("Operating empty weight estimate for C-5a: "+str(weight))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location = vehicle.mass_properties.center_of_gravity
    print("C-5a CG location: "+str(CG_location))
    
    # ------------------------------------------------------------------
    #   Aircraft MOI
    # ------------------------------------------------------------------    
    MOI = calculate_aircraft_MOI(vehicle, CG_location)

    print(MOI)

    
    
if __name__ == '__main__': 
    main()
    plt.show()    