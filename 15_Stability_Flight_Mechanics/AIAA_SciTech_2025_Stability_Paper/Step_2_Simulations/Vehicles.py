# Vehicle.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------     

import RCAIDE
from RCAIDE.Framework.Core                              import Units    
from RCAIDE.Library.Plots                               import *
import numpy as np
import os


################################################################################################################################################
# STICK FIXED 
################################################################################################################################################
def stick_fixed_stability_setup(vehicle):  
    configs  = stick_fixed_stability_configs_setup(vehicle) 
    return configs 
 

def stick_fixed_stability_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle) 
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'stick_fixed_cruise'
    configs.append(config) 
    return configs  
