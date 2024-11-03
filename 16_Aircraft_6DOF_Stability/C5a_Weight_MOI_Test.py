# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units
from RCAIDE.Library.Plots                 import *     
from RCAIDE.Library.Methods.Weights.Correlation_Buildups import Common
from RCAIDE.Library.Methods.Stability.Moment_of_Inertia.calculate_aircraft_MOI import calculate_aircraft_MOI
from RCAIDE.Library.Methods.Stability.Center_of_Gravity     import compute_vehicle_center_of_gravity


# python imports 
import numpy as np  
import matplotlib.pyplot as plt  
import sys
from RCAIDE.Library.Methods.Stability.Moment_of_Inertia import compute_cuboid_moment_of_inertia
sys.path.append('../Aircraft/C-5a')
import Lockheed_C5a

def main():
    vehicle = Lockheed_C5a.vehicle_setup()
    
    # update fuel weight to 60%
    vehicle.networks.fuel.fuel_lines.fuel_line.fuel_tanks.wing_fuel_tank.fuel.mass_properties.mass = 0.6 * vehicle.networks.fuel.fuel_lines.fuel_line.fuel_tanks.wing_fuel_tank.fuel.mass_properties.mass

    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    weight_analysis.method                        = 'Raymer'
    weight_analysis.settings.use_max_fuel_weight  = False  
    results                                       = weight_analysis.evaluate() 
    print("Operating empty weight estimate for C-5a: "+str(results))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location      = vehicle.mass_properties.center_of_gravity
    CG_location_true = np.array([[29.5, 0, 0.547]])
    print("C-5a CG location: "+str(CG_location))
    
    # ------------------------------------------------------------------
    #   Operating Aircraft MOI
    # ------------------------------------------------------------------    
    MOI, total_mass = calculate_aircraft_MOI(vehicle, CG_location)

    # ------------------------------------------------------------------
    #   Payload MOI
    # ------------------------------------------------------------------    
    Cargo_MOI, mass =  compute_cuboid_moment_of_inertia(CG_location, 99790*Units.kg, 36.0, 3.66, 3, 0, 0, 0, CG_location)
    MOI             += Cargo_MOI
    total_mass      += mass

    print(MOI)
    print("MOI Mass: " + str(total_mass))
    sft2     = 1.355817
    C5a_true = np.array([[27800000 , 0, 2460000], [0, 31800000, 0], [2460000, 0, 56200000]]) * sft2
    error    = (MOI - C5a_true) / C5a_true * 100
    print(error)
    
if __name__ == '__main__': 
    main()
    plt.show()    