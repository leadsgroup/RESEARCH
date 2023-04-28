# Procedure.py
# 
# Created:  Mar 2023, M Clarke

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------      
import MARC
from MARC.Analyses.Process import Process  
import  Missions
import numpy as np 
from MARC.Core                                                              import Data , Units 
from MARC.Components.Energy.Networks.Battery_Electric_Rotor                 import Battery_Electric_Rotor
from MARC.Methods.Propulsion                                                import propeller_design 
from MARC.Methods.Power.Battery.Sizing                                      import initialize_from_circuit_configuration  
from MARC.Methods.Propulsion.electric_motor_sizing                          import size_optimal_motor
from MARC.Methods.Weights.Correlations.Propulsion                           import nasa_motor
from MARC.Methods.Geometry.Two_Dimensional.Planform                         import segment_properties
from MARC.Methods.Weights.Buildups.eVTOL.empty                              import empty  
from MARC.Methods.Geometry.Two_Dimensional.Planform.wing_segmented_planform import wing_segmented_planform
from MARC.Methods.Weights.Buildups.eVTOL.converge_evtol_weight              import converge_evtol_weight  
from MARC.Methods.Center_of_Gravity.compute_component_centers_of_gravity    import compute_component_centers_of_gravity

# ----------------------------------------------------------------------        
#   Setup
# ----------------------------------------------------------------------       

def setup():
    procedure = Process() 
    procedure.modify_vehicle_and_mission  = modify_vehicle_and_mission
    procedure.finalize                    = finalize 
    procedure.missions                    = Process()
    procedure.missions.base               = evaluate_mission 
    procedure.post_process                = post_process
        
    return procedure 

# ----------------------------------------------------------------------        
#   Evaluate Mission
# ----------------------------------------------------------------------    
    
def evaluate_mission(nexus): 
    # Evaluate the missions and save to results    
    mission         = nexus.missions.base
    results         = nexus.results
    results.base   = mission.evaluate()                                        
    return nexus
 
# ----------------------------------------------------------------------        
#  Modify Vehicle
# ----------------------------------------------------------------------    

def modify_vehicle_and_mission(nexus): 
    # this function actually does nothing but here is where you would modify the vehicle/mission
    # segments that rely on the mission  
    # Pull out the vehicle
    vehicle_opt = nexus.vehicle_configurations.base
   
    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------

    # mass properties
    wing = vehicle_opt.wings.main_wing
    wing.aspect_ratio    = wing.spans.projected**2. / wing.areas.reference 
    
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------
    
    
    # Fill out more segment properties automatically 
    wing_segmented_planform(wing, overwrite_reference = True ) 
    wing = segment_properties(wing)
    vehicle_opt.reference_area              = wing.areas.reference    
    net = vehicle_opt.networks.battery_electric_rotor
    
    # Component 3: Battery       
    bat                                                    = MARC.Components.Energy.Storages.Batteries.Constant_Mass.Lithium_Ion_LiNiMnCoO2_18650() 
    bat.pack.electrical_configuration.series               =  np.ceil(net.battery.pack.electrical_configuration.series)
    bat.pack.electrical_configuration.parallel             =  np.ceil(net.battery.pack.electrical_configuration.parallel)
    initialize_from_circuit_configuration(bat)  
    bat.module_config.number_of_modules                    = 14  
    bat.module.geometrtic_configuration.total              = bat.pack.electrical_configuration.total
    bat.module_config.voltage                              = bat.pack.max_voltage/bat.module_config.number_of_modules # assumes modules are connected in parallel, must be less than max_module_voltage (~50) /safety_factor (~ 1.5)  
    bat.module.geometrtic_configuration.normal_count       = bat.pack.electrical_configuration.series/bat.module_config.number_of_modules
    bat.module.geometrtic_configuration.parallel_count     = bat.pack.electrical_configuration.parallel/bat.module_config.number_of_modules     
    net.voltage                                            = bat.pack.max_voltage
     
    settings = Data()
    converge_evtol_weight(vehicle_opt,settings,contingency_factor = 1.1) 
    breakdown = empty(vehicle_opt,settings,contingency_factor     = 1.1 ) 
    
    vehicle_opt.weight_breakdown  = breakdown
    compute_component_centers_of_gravity(vehicle_opt)
    vehicle_opt.center_of_gravity()    
    # ----------------------------------------------------------------------
    # Update Mission 
    # ---------------------------------------------------------------------- 
    # Re-set nattery charge each optimization
    nexus.missions.base.segments.takeoff.battery_energy =  vehicle_opt.networks.battery_electric_rotor.battery.pack.max_energy 
    
    # diff the new data
    vehicle_opt.store_diff()
    
    return nexus 
 

# ----------------------------------------------------------------------
#   Finalizing Function
# ----------------------------------------------------------------------    

def finalize(nexus):
    
    nexus.analyses.finalize()   
    
    return nexus         

# ----------------------------------------------------------------------
#   Post Process Results to give back to the optimizer
# ----------------------------------------------------------------------   

def post_process(nexus):
    
    # Unpack data  
    res       = nexus.results.base
    final_DOC = 1-res.segments[-1].conditions.propulsion.battery.cell.state_of_charge[-1,0]
    
    summary         = nexus.summary     
    summary.Nothing = 0.0 
    summary.SOC_EOF = final_DOC           
  
    nexus.total_number_of_iterations +=1
    return nexus 