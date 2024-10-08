# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units ,Data  
from RCAIDE.Library.Plots  import *       
from RCAIDE.Library.Methods.Propulsors.Converters.Ducted_Fan.design_ducted_fan import design_ducted_fan
# python imports  
import matplotlib.pyplot as plt  
import os  

 
def  main():
     

    ospath                                = os.path.abspath(__file__)
    separator                             = os.path.sep
    rel_path                              = ospath.split()[0]  +  '..' + separator + '..' + separator
    
    
    ducted_fan                            = RCAIDE.Library.Components.Propulsors.Converters.Ducted_Fan()
    ducted_fan.tag                        = 'test_ducted_fan'
    ducted_fan.design_thrust              = 50 #
    ducted_fan.design_altitude            = 3000 * Units.feet
    ducted_fan.design_angular_velocity    = 8000 *Units.rpm
    ducted_fan.number_of_blades           = 11
    ducted_fan.number_of_radial_stations  = 20
    ducted_fan.design_freestream_velocity = 50  
    ducted_fan.design_reference_velocity  = 50  
    airfoil                               = RCAIDE.Library.Components.Airfoils.Airfoil()
    airfoil.tag                           = 'NACA_4412' 
    airfoil.coordinate_file               =  rel_path + 'Airfoils' + separator + 'NACA_4412.txt'    
    #ducted_fan.append_airfoil(airfoil)
    design_ducted_fan(ducted_fan)
    
    
    return 

if __name__ == '__main__': 
    main()
    plt.show()    