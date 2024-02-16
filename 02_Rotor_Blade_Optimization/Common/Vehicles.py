# Vehicles.py
# 
# Created: Jun 2022, E. Botero
# Modified: 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
from RCAIDE.Components.Energy.Networks import Battery_Electric_Rotor
from RCAIDE.Methods.Propulsion import propeller_design 
from RCAIDE.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil.compute_airfoil_properties \
     import compute_airfoil_properties
from RCAIDE.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil.import_airfoil_geometry \
     import import_airfoil_geometry     
from RCAIDE.Analyses.Propulsion.Rotor_Wake_Fidelity_Zero import Rotor_Wake_Fidelity_Zero
import numpy as np

from RCAIDE.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil.import_airfoil_geometry  import import_airfoil_geometry
from RCAIDE.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil.compute_airfoil_properties import compute_airfoil_properties


# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

def setup(problem):
    
    base_vehicle = base_setup(problem)
    configs = configs_setup(base_vehicle)
    
    return configs

def base_setup(problem):

    #----------------------------------------------------------------------
    # Use propeller design to setup a rotor
    #----------------------------------------------------------------------  
    
    rs      = problem.rotor_settings
    
    vehicle = RCAIDE.Vehicle()
    net     = Battery_Electric_Rotor()
    
    prop = RCAIDE.Components.Energy.Converters.Propeller()
    prop.tag                      = 'prop_rotor'
    prop.orientation_euler_angles = rs.orientation_euler_angles
    prop.number_of_blades         = rs.number_of_blades
    prop.tip_radius        = R    = rs.tip_radius
    prop.hub_radius        = Rh   = rs.hub_radius
    prop.design_freestream_velocity      = rs.hover.design_freestream_velocity
    prop.design_angular_velocity         = rs.design_angular_velocity
    prop.design_altitude          = rs.hover.design_altitude
    prop.design_thrust            = rs.hover.design_thrust 
    prop.design_Cl                = rs.design_Cl
    prop.airfoil_polar_stations   = rs.airfoil_polar_stations
    prop.airfoil_geometry         = rs.airfoil_geometry
    prop.airfoil_polars           = rs.airfoil_polars
    prop.number_of_stations   = N = rs.number_of_stations    

    prop.Wake = Rotor_Wake_Fidelity_Zero()  

    prop.number_of_airfoil_section_points               = 500
    
    # rotor surrogates and geometry data 
    airfoil_geometry_data                               = import_airfoil_geometry(prop.airfoil_geometry[0],npoints = prop.number_of_airfoil_section_points) 
    airfoil_polar_data                                  = compute_airfoil_properties(airfoil_geometry_data, airfoil_polar_files= prop.airfoil_polars[0],use_pre_stall_data=True,linear_lift=True ) 
     
    prop.RE_data                                        = airfoil_polar_data.reynolds_numbers
    prop.aoa_data                                       = airfoil_polar_data.angle_of_attacks
    prop.airfoil_cl_surrogates                          = airfoil_polar_data.lift_coefficients
    prop.airfoil_cd_surrogates                          = airfoil_polar_data.drag_coefficients
    prop.airfoil_bl_aoa_data                            = airfoil_polar_data.boundary_layer_angle_of_attacks           
    prop.airfoil_bl_RE_data                             = airfoil_polar_data.boundary_layer_reynolds_numbers   
    prop.airfoil_bl_lower_surface_theta_surrogates      = airfoil_polar_data.boundary_layer_theta_lower_surface           
    prop.airfoil_bl_lower_surface_delta_surrogates      = airfoil_polar_data.boundary_layer_delta_lower_surface           
    prop.airfoil_bl_lower_surface_delta_star_surrogates = airfoil_polar_data.boundary_layer_delta_star_lower_surface          
    prop.airfoil_bl_lower_surface_Ue_surrogates         = airfoil_polar_data.boundary_layer_cf_lower_surface          
    prop.airfoil_bl_lower_surface_cf_surrogates         = airfoil_polar_data.boundary_layer_Ue_Vinf_lower_surface              
    prop.airfoil_bl_lower_surface_dp_dx_surrogates      = airfoil_polar_data.boundary_layer_dcp_dx_lower_surface          
    prop.airfoil_bl_upper_surface_theta_surrogates      = airfoil_polar_data.boundary_layer_theta_upper_surface           
    prop.airfoil_bl_upper_surface_delta_surrogates      = airfoil_polar_data.boundary_layer_delta_upper_surface           
    prop.airfoil_bl_upper_surface_delta_star_surrogates = airfoil_polar_data.boundary_layer_delta_star_upper_surface            
    prop.airfoil_bl_upper_surface_Ue_surrogates         = airfoil_polar_data.boundary_layer_cf_upper_surface          
    prop.airfoil_bl_upper_surface_cf_surrogates         = airfoil_polar_data.boundary_layer_Ue_Vinf_upper_surface
    prop.airfoil_bl_upper_surface_dp_dx_surrogates      = airfoil_polar_data.boundary_layer_dcp_dx_upper_surface
    prop.airfoil_bl_upper_surface_dp_dx_surrogates      = airfoil_polar_data.boundary_layer_dcp_dx_upper_surface    
     
    # simple rotor planform 
    chi0                                                = prop.hub_radius/prop.tip_radius  
    chi                                                 = np.linspace(chi0,1,prop.number_of_stations+1)  
    chi                                                 = chi[0:prop.number_of_stations]   
    prop.radius_distribution                            = chi*prop.tip_radius    
    prop.twist_distribution                             = np.linspace(np.pi/4,np.pi/6,prop.number_of_stations) 
    prop.chord_distribution                             = np.ones_like(prop.twist_distribution)*0.2
    prop.thickness_to_chord                             = np.ones_like(prop.twist_distribution)* airfoil_geometry_data.thickness_to_chord 
    prop.max_thickness_distribution                     = prop.thickness_to_chord *prop.chord_distribution 
    prop.mid_chord_alignment                            = np.zeros_like(prop.chord_distribution)         
    prop.airfoil_flag                                   = True      
    prop.airfoil_geometry_data                          = airfoil_geometry_data
    
    net.rotors.append(prop)
    vehicle.append_component(net) 

    return vehicle


def configs_setup(base_vehicle):
    
    configs                             = RCAIDE.Components.Configs.Config.Container()
    base_config                         = RCAIDE.Components.Configs.Config(base_vehicle)
    
    config                              = RCAIDE.Components.Configs.Config(base_config)
    config.tag                          = 'cruise'
    configs.append(config)
    
    config                              = RCAIDE.Components.Configs.Config(base_config)
    config.tag                          = 'hover' 
    config.networks.Battery_Electric_Rotor.rotors.prop_rotor.orientation_euler_angles = [0.,np.pi/2.,0.]
    configs.append(config)      

    
    config                              = RCAIDE.Components.Configs.Config(base_config)
    config.tag                          = 'OEI' 
    config.networks.Battery_Electric_Rotor.rotors.prop_rotor.orientation_euler_angles = [0.,np.pi/2.,0.]
    configs.append(config)          
    
    
    return configs
