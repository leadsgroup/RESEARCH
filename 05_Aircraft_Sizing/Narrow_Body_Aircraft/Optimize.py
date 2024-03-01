# Optimize.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    
 
from RCAIDE.Core import Units, Data 
from RCAIDE.Optimization.Common             import Nexus , carpet_plot
from RCAIDE.Optimization.Packages.scipy     import scipy_setup


import Vehicles
import Analyses
import Missions
import Procedure
import Plot_Mission 

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------        
#   Run the whole thing
# ----------------------------------------------------------------------  
def main():
    
    problem = setup()
    
    ## Base Input Values
    output = problem.objective()
    
    # Uncomment to view contours of the design space
    #variable_sweep(problem)
    
    # Uncomment for the first optimization
    output = scipy_setup.SciPy_Solve(problem,solver='SLSQP')
    print (output)    

    print('fuel burn = ', problem.summary.base_mission_fuelburn)
    print('fuel margin = ', problem.all_constraints())
    
    Plot_Mission.plot_mission(problem)
    
    return

# ----------------------------------------------------------------------        
#   Inputs, Objective, & Constraints
# ----------------------------------------------------------------------  

def setup():

    nexus  = Nexus()
    problem = Data()

    # -------------------------------------------------------------------
    # Inputs (Design Variables)
    # -------------------------------------------------------------------

    #   [ tag                   , initial,     lb , ub        , scaling , units ]
    problem.inputs = np.array([
        [ 'wing_area'           ,  92    ,    50. ,   130.    ,   100.  , 1*Units.meter**2],
        [ 'cruise_altitude'     ,   8    ,     6. ,    12.    ,   10.   , 1*Units.km],
    ],dtype=object)

    # -------------------------------------------------------------------
    # Objective
    # -------------------------------------------------------------------

    # [ tag, scaling, units ]
    problem.objective = np.array([
        [ 'fuel_burn', 10000, 1*Units.kg ]
    ],dtype=object)
    
    # -------------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------------
    
    # [ tag, sense, edge, scaling, units ]
    problem.constraints = np.array([
        [ 'design_range_fuel_margin' , '>', 0., 1E-1, 1*Units.less], #fuel margin defined here as fuel 
    ],dtype=object)
    
    # -------------------------------------------------------------------
    #  Aliases
    # -------------------------------------------------------------------
    
    # [ 'alias' , ['data.path1.name','data.path2.name'] ]

    problem.aliases = [
        [ 'wing_area'                        ,   ['vehicle_configurations.*.wings.main_wing.areas.reference',
                                                  'vehicle_configurations.*.reference_area'                    ]],
        [ 'cruise_altitude'                  , 'missions.base.segments.climb_5.altitude_end'                    ],
        [ 'fuel_burn'                        ,    'summary.base_mission_fuelburn'                               ],
        [ 'design_range_fuel_margin'         ,    'summary.max_zero_fuel_margin'                                ],
    ]    
    

    nexus.optimization_problem = problem    
    
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
    nexus.missions = Missions.setup(nexus.analyses)
    
    # -------------------------------------------------------------------
    #  Procedure
    # -------------------------------------------------------------------    
    nexus.procedure = Procedure.setup()
    
    # -------------------------------------------------------------------
    #  Summary
    # -------------------------------------------------------------------    
    nexus.summary = Data()    
    nexus.total_number_of_iterations = 0
    return nexus
    
def variable_sweep(problem):    
    number_of_points = 5
    outputs     = carpet_plot(problem, number_of_points, 0, 0)  #run carpet plot, suppressing default plots
    inputs      = outputs.inputs
    objective   = outputs.objective
    constraints = outputs.constraint_val
    plt.figure(0)
    CS   = plt.contourf(inputs[0,:],inputs[1,:], objective, 20, linewidths=2,cmap='jet')
    cbar = plt.colorbar(CS)
    
    cbar.ax.set_ylabel('fuel burn (kg)')
    CS_const = plt.contour(inputs[0,:],inputs[1,:], constraints[0,:,:],cmap='jet')
    plt.clabel(CS_const, inline=1, fontsize=10)
    cbar = plt.colorbar(CS_const)
    cbar.ax.set_ylabel('fuel margin')
    
    plt.xlabel('Wing Area (m^2)')
    plt.ylabel('Cruise Altitude (km)')
    
    plt.show(block=True)    
    
    return

if __name__ == '__main__':
    main()
    plt.show()
