# --------------------------------------
# Imports
# --------------------------------------
import RCAIDE
from   RCAIDE.Framework.Core import Units
import os 
import pickle  

def main():
    
    # ----------------------------------------------------------------------------
    # -1  load optimized_vehicle
    # ----------------------------------------------------------------------------    
    
    Optimized_Vehicle_1_pkl               = "025_00_00_40_00_00_Optimized_Vehicle"
    Optimized_Vehicle_1                   = load_results(Optimized_Vehicle_1_pkl)
    
    # ----------------------------------------------------------------------------
    # -2  create aircraft with any elevator 
    # ----------------------------------------------------------------------------     
    
    elevator                              = RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator()
    elevator.tag                          = 'elevator'
    elevator.span_fraction_start          = 0.1
    elevator.span_fraction_end            = 0.9
    elevator.deflection                   = 0.0  * Units.deg
    elevator.chord_fraction               = 0.35
    Optimized_Vehicle_1.wings.horizontal_tail.control_surfaces.append(elevator)     
    
    # ----------------------------------------------------------------------------
    # -3  run 1st mission: Pull up 
    # ---------------------------------------------------------------------------- 
    
    


    
    
    
    
    
    
    debug = 0
    
    return



def load_results(filename):  

    #current_dir = '/Users/aidanmolloy/Documents/LEADS/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations'
    current_dir = 'C:/Users/Matteo/Documents/UIUC/RESEARCH/15_Stability_Flight_Mechanics/AIAA_SciTech_2025_Stability_Paper/Step_1_Simulations'
    load_file = os.path.join(current_dir, filename + '.pkl')
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    
    return results

if __name__ == '__main__': 
    main()
    plt.show()