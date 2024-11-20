# Vehicle.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------     

import RCAIDE
from RCAIDE.Framework.Core                              import Units    
from RCAIDE.Library.Plots                               import *
import numpy as np
import os

def stick_fixed_stability_setup(vehicle):  
    configs  = stick_fixed_stability_configs_setup(vehicle) 
    return configs 
 
 
def elevator_sizing_setup(vehicle):   
    ht_wing                        = vehicle.wings.horizontal_tail
    elevator                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                   = 'elevator'
    elevator.span_fraction_start   = 0.1
    elevator.span_fraction_end     = 0.9
    elevator.deflection            = 0.0  * Units.deg
    elevator.chord_fraction        = 0.3
    ht_wing.append_control_surface(elevator)     
    configs                        = elevator_sizing_configs_setup(vehicle) 
    return configs

def aileron_rudder_sizing_setup(vehicle):    
    mw_wing                       = vehicle.wings.main_wing 
    aileron                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Aileron()
    aileron.tag                   = 'aileron'
    aileron.span_fraction_start   = 0.7
    aileron.span_fraction_end     = 0.9 
    aileron.deflection            = 0.0 * Units.degrees
    aileron.chord_fraction        = 0.2
    mw_wing.append_control_surface(aileron) 
    
    if vehicle.rudder_flag:
        vt_wing                      = vehicle.wings.vertical_tail
        rudder                       = RCAIDE.Library.Components.Wings.Control_Surfaces.Rudder()
        rudder.tag                   = 'rudder'
        rudder.span_fraction_start   = 0.2
        rudder.span_fraction_end     = 0.8
        rudder.deflection            = 0.0  * Units.deg
        rudder.chord_fraction        = 0.2
        vt_wing.append_control_surface(rudder) 
    
    configs  = aileron_rudder_sizing_configs_setup(vehicle) 
    return configs 
 
def flap_sizing_setup(vehicle):   
    mw_wing                       = vehicle.wings.main_wing
    flap                          = RCAIDE.Library.Components.Wings.Control_Surfaces.Flap()
    flap.tag                      = 'flap'
    flap.span_fraction_start      = 0.2
    flap.span_fraction_end        = 0.5
    flap.deflection               = 0.0 * Units.degrees 
    flap.chord_fraction           = 0.20
    mw_wing.append_control_surface(flap)   
    configs                       = flap_sizing_configs_setup(vehicle) 
    return configs 

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def stick_fixed_stability_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle) 
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'stick_fixed_cruise'
    configs.append(config) 
    return configs  

def elevator_sizing_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container()
    
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)
    
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'elevator_sizing_pull_up'   
    configs.append(config)

    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'elevator_sizing_push_over'   
    configs.append(config)
    
    return configs

def aileron_rudder_sizing_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)  
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'aileron_rudder_roll_sizing'   
    configs.append(config)

    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'aileron_rudder_crosswind_sizing'   
    configs.append(config)       
    return configs  

def flap_sizing_configs_setup(vehicle): 
    configs     = RCAIDE.Library.Components.Configs.Config.Container()  
    base_config = RCAIDE.Library.Components.Configs.Config(vehicle)  
    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'flap_sizing_flaps_up'   
    configs.append(config)
    

    config      = RCAIDE.Library.Components.Configs.Config(base_config)
    config.tag  = 'flap_sizing_flaps_down'   
    configs.append(config)       
    return configs
