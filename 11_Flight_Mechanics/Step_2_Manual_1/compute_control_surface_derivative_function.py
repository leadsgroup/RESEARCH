
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from   RCAIDE.Framework.Core import Units,Data  
from   RCAIDE.Library.Plots  import *

import numpy as np
import os   
import pickle  
import pylab as plt 
from   copy import deepcopy
import sys 
import pandas as pd

 
def compute_control_surface_derivative_function(vehicle, velocity,altitude, control_surface_tags):
    'creates function of control surface delfection and angle'
    
    deflection_angles     =  np.linspace(-30, 30, 7)
    control_surface_size  =  np.linspace(0.1,0.9,10)
    angle_of_attack_range =  np.array([[0]])

    # find control surface
    for cs_tag in range(len(control_surface_tags)):
        for wing in  vehicle.wings:
            for control_surface in  wing.control_surfaces:
                if control_surface == cs_tag:
                    '''  modify control surface'''
                    
                    # define empty arrays to store surrogate raw data             
                    for def_ang in range(len(deflection_angles)):
                        for cs_surf_size in  range(len(control_surface_size)):
                            
                            
                            # modify control surface on aircraft 
                            
                            
                            
                            
                            # run aircraft at flight condition and store derivatives
                            #------------------------------------------------------------------------
                            # setup flight conditions
                            #------------------------------------------------------------------------   
                            atmosphere     = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
                            atmo_data      = atmosphere.compute_values(altitude)
                            P              = atmo_data.pressure 
                            T              = atmo_data.temperature 
                            rho            = atmo_data.density 
                            a              = atmo_data.speed_of_sound 
                            mu             = atmo_data.dynamic_viscosity
                               
                            # -----------------------------------------------------------------
                            # Evaluate Without Surrogate
                            # ----------------------------------------------------------------- 
                            ctrl_pts = len(angle_of_attack_range[:, 0] )
                            state                                         = RCAIDE.Framework.Mission.Common.State()
                            state.conditions                              = RCAIDE.Framework.Mission.Common.Results() 
                            state.conditions.freestream.density           = rho * np.ones_like(angle_of_attack_range)
                            state.conditions.freestream.dynamic_viscosity = mu  * np.ones_like(angle_of_attack_range)
                            state.conditions.freestream.temperature       = T   * np.ones_like(angle_of_attack_range)
                            state.conditions.freestream.pressure          = P   * np.ones_like(angle_of_attack_range)
                            state.conditions.aerodynamics.angles.alpha    = angle_of_attack_range  
                            state.conditions.aerodynamics.angles.beta     = angle_of_attack_range *0  
                            state.conditions.freestream.u                 = angle_of_attack_range *0       
                            state.conditions.freestream.v                 = angle_of_attack_range *0       
                            state.conditions.freestream.w                 = angle_of_attack_range *0       
                            state.conditions.static_stability.roll_rate   = angle_of_attack_range *0       
                            state.conditions.static_stability.pitch_rate  = angle_of_attack_range *0 
                            state.conditions.static_stability.yaw_rate    = angle_of_attack_range *0  
                            state.conditions.expand_rows(ctrl_pts)
                          
                         
                            state.analyses                                  =  Data()
                            aerodynamics                                    = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
                            aerodynamics.settings.use_surrogate             = False 
                            aerodynamics.vehicle                            = vehicle
                            aerodynamics.settings.model_fuselage            = True   
                            aerodynamics.initialize()
                            state.analyses.aerodynamics = aerodynamics 
                             
                            state.conditions.freestream.mach_number                 = velocity / a 
                            state.conditions.freestream.velocity                    = velocity
                            state.conditions.freestream.reynolds_number             = state.conditions.freestream.density * state.conditions.freestream.velocity / state.conditions.freestream.dynamic_viscosity 
                            state.conditions.frames.inertial.velocity_vector[:,0]   = velocity
                            
                         
                            # ---------------------------------------------------------------------------------------
                            # Evaluate With Surrogate
                            # ---------------------------------------------------------------------------------------  
                            _                 = state.analyses.aerodynamics.evaluate(state)            
                            
                            
                            # store raw data for surrogate
                            
                            
                            
                   
                    # reset control surface to zero deflection so it does not affect results
                    
                    # create surrogate and store onto aircraft control surface data structure
                    control_surface.deflection_and_size_surrogate =  ...
                 
    return 