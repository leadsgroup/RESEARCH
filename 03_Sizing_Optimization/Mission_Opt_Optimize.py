#Mission_Optimize.py
# Will Kupiec
# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import MARC
assert MARC.__version__=='1.0.0', 'These tutorials only work with the MARC 1.0.0 release'
from MARC.Core import Units, Data
import numpy as np
import matplotlib.pyplot as plt
import Mission_Opt_Vehicle
import Mission_Opt_Analyses
import Mission_Opt_Missions
import Mission_Opt_Procedure
import Plot_Mission
from MARC.Optimization import Nexus 
import MARC.Optimization.Package_Setups.scipy_setup as scipy_setup 

# ----------------------------------------------------------------------        
#   Run the whole thing
# ----------------------------------------------------------------------  
def main():
    
    problem = setup()  
    output = scipy_setup.SciPy_Solve(problem,solver='SLSQP', sense_step = 1E-2, tolerance = 1E-3)
    print (output)     
    
    Plot_Mission.plot_mission(problem)
    
    return 
 

def setup(): 

    nexus   = Nexus()
    problem = Data()
    nexus.optimization_problem = problem

    # -------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------

    #   [ tag                   , initial,  (lb , ub)    , scaling , units ]  
    problem.inputs = np.array([
                      ['climbRate'  , 100, 75, 125, 100., 1*Units.less],
                      ['descentRate', 140, 100, 180, 100., 1*Units.less],
        ],dtype=object)    

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([
        [ 'SOC_end_of_flight', 1. , 1*Units.less ] 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([
        [ 'SOC_end_of_flight' , '>' , 0.2 ,  1.  , 1*Units.less], 
         
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ]

    problem.aliases = [
        [ 'climbRate'                 , 'missions.base.segments.climb.climb_rate'],
        [ 'descentRate'               , 'missions.base.segments.descent.climb_rate'],
        [ 'SOC_end_of_flight'         , 'summary.SOC_EOF'],
        [ 'Nothing'                   , 'summary.Nothing'],               
    ]   
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Mission_Opt_Vehicle.setup()
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Mission_Opt_Analyses.setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Mission_Opt_Missions.setup(nexus.analyses)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure =  Mission_Opt_Procedure.setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()    
    nexus.total_number_of_iterations = 0
    return nexus    

if __name__ == '__main__':
    main()
    plt.show()
