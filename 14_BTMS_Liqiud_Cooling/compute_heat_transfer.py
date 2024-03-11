# @ingroup Methods-Energy-Sources-Battery-Common
# RCAIDE/Methods/Energy/Sources/Battery/Common/compute_heat_transfer.py
# 
# 
# Created:  Feb 2024, M. Clarke 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------
from RCAIDE.Core   import Units , Data 
import numpy as np  
 
# ----------------------------------------------------------------------------------------------------------------------
# compute_heat transfer
# ---------------------------------------------------------------------------------------------------------------------- 
## @ingroup Energy-Sources-Batteries-Common

def compute_heat_transfer(battery_conditions,battery,state,dt,i):
    
    # Battery Properties
    d_cell                      = battery.cell.diameter                                         
    A_surface_cell              = battery.cell.surface_area  
    Cp                          = battery.cell.specific_heat_capacity     
    k                           = battery.cell.axial_thermal_conductivity 
    emissivity_battery          = battery.cell.emissivity
    mass                        = battery.cell.mass 
    
    
    # Temperature 
    T_cell                      = battery_conditions.cell.temperature[i+1] 
    T_ambient                   = state.conditions.freestream.temperature[i,0] 
    
    # Heat Transfer properties
    sigma                       = 5.69e-8   #Stefan Boltzman Constant
    h                           = 2.5       #[W/m^2-k]
    emissivity_air              = 0.9
    
    # Heat Transfer due to conduction. 
    dQ_dt_cond                  = k*A_surface_cell*(T_cell-T_ambient)/(d_cell/2)
    
    # Heat Transfer dur to natural convention 
    dQ_dt_conv                  = h*A_surface_cell*(T_cell-T_ambient)
    
    # Heat Transfer due to radiation 
    dQ_dt_rad                   = sigma*A_surface_cell*((emissivity_battery*T_cell**4)-(emissivity_air*T_ambient**4)) # xheck the equation 
    
    dQ_dt                       = (dQ_dt_cond+dQ_dt_conv+dQ_dt_rad)
    
    dT_dt                       = dQ_dt/(mass*Cp)
     
    T_current                   = T_cell + dT_dt*dt[i] 
     
    
    battery_conditions.cell.temperature[i+1]   = T_current
    

    
    return