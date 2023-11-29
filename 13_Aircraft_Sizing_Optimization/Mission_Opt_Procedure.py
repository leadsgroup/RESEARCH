# Procedure.py
# 
# Created:  Mar 2023, M Clarke

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------      
import MARC
from MARC.Analyses.Process import Process  
import  Mission_Opt_Missions 
import  Config_Opt_Optimize
import numpy as np 
import MARC.Optimization.Package_Setups.scipy_setup as scipy_setup 
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
    procedure.modify_mission              = modify_mission
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

def modify_mission(nexus): 
    
    # run the configuration optimization 
    config_opt_problem = Config_Opt_Optimize.setup()  
    config_opt_outputs = scipy_setup.SciPy_Solve(config_opt_problem,solver='SLSQP', sense_step = 1E-2, tolerance = 1E-3)    
    
    
    # use outputs of config optimization to update vehicle of mission optimiztion  

    mission_opt_vehicle  = nexus.vehicle_configurations.base    
    mission_opt_vehicle.span = config_opt_outputs[0] 
    ... 
    
    # ----------------------------------------------------------------------
    # Update Mission 
    # ---------------------------------------------------------------------- 
    # Re-set nattery charge each optimization
    nexus.missions.base.segments.takeoff.battery_energy =  mission_opt_vehicle.networks.battery_electric_rotor.battery.pack.max_energy 
    
    # diff the new data
    mission_opt_vehicle.store_diff()
    
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