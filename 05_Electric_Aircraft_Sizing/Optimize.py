# Optimize.py
# Created:  Mar 2023, M. Clarke

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import MARC
assert MARC.__version__=='1.0.0', 'These tutorials only work with the MARC 1.0.0 release'
from MARC.Core import Units, Data
import numpy as np
import Vehicles
import Analyses
import Missions
import Procedure
import Plot_Mission 
from MARC.Optimization import Nexus 
import MARC.Optimization.Package_Setups.scipy_setup as scipy_setup 

# ----------------------------------------------------------------------        
#   Run the whole thing
# ----------------------------------------------------------------------  
def main():
    
    problem = setup()  
    output = scipy_setup.SciPy_Solve(problem,solver='SLSQP')
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
                  [ 'cruise_range' , 30   , 0.01 , 300.  ,   10. , 1*Units.nmi      ],  
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
        [ 'SOC_end_of_flight'  ,   '>'   , 0.2     ,   1.  , 1*Units.less], 
         
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ]

    problem.aliases = [
        [ 'cruise_range'              , 'missions.base.segments.cruise.distance'   ],  
        [ 'SOC_end_of_flight'         , 'summary.SOC_EOF'        ],
        [ 'Nothing'                   , 'summary.Nothing'    ],               
    ]    
    
    # -------------------------------------------------------------------
    #  Vehicles
    # -------------------------------------------------------------------
    nexus.vehicle_configurations = Vehicles.setup()
    
    # -------------------------------------------------------------------
    #  Analyses
    # -------------------------------------------------------------------
    nexus.analyses = Analyses.setup(nexus.vehicle_configurations)
    
    # -------------------------------------------------------------------
    #  Missions
    # -------------------------------------------------------------------
    nexus.missions = Missions.setup(nexus.analyses,nexus.vehicle_configurations.opt)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure =  Procedure.setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()    
    nexus.total_number_of_iterations = 0
    return nexus    

if __name__ == '__main__':
    main()
