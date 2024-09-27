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

# python imports 
import numpy as np  
from copy import deepcopy
import matplotlib.pyplot as plt  
import os

def main():
    vehicle = Lockheed_C5a.vehicle_setup()
    
    settings = None    
    weight = Common.compute_operating_empty_weight(vehicle, settings = settings, method_type = 'RCAIDE')
    print("Operating empty weight estimate for C-5a: "+str(weight))


    
    
    
    
    
    
    
    
if __name__ == '__main__': 
    main()
    plt.show()    