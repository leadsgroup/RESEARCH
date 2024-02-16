# Procedure.py
# 
# Created:  Mar 2023, M Clarke

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------      
from RCAIDE.Analyses.Process import Process  

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
    results.mission = mission.evaluate()                                        
    return nexus
 
# ----------------------------------------------------------------------        
#  Modify Vehicle
# ----------------------------------------------------------------------    

def modify_vehicle_and_mission(nexus): 
    # this function actually does nothing but here is where you would modify the vehicle/mission
    # segments that rely on the mission  
    return  
 

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
    res       = nexus.results.mission    
    final_SOC = res.segments[-1].conditions.propulsion.battery.pack.state_of_charge[-1]   
    
    summary         = nexus.summary     
    summary.Nothing = 0.0 
    summary.SOC_EOF = final_SOC           
  
    nexus.total_number_of_iterations +=1
    return nexus 