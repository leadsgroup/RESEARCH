# multi_airfoil_panel_code.py
 
# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
# RCAIDE Imports 
from RCAIDE.Core import Units
from RCAIDE.Methods.Aerodynamics.Airfoil_Panel_Method     import airfoil_analysis  
from RCAIDE.Methods.Geometry.Two_Dimensional.Airfoil      import import_airfoil_geometry
from RCAIDE.Visualization import * 
from RCAIDE.Methods.Geometry.Two_Dimensional.Airfoil      import  compute_naca_4series

# Python imports
import os 
import numpy as np
import matplotlib.pyplot as plt    
 

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():    
    # path import 
    ospath    = os.path.abspath(__file__)
    separator = os.path.sep 
    rel_path  = ospath.split('multi_airfoil_panel_code.py')[0] +  'Airfoils' + separator 
    
    # operating conditions
    nu              = 1.5*10**(-5)     # Air kinematic viscosity
    RPM             = 2400
    Omega           = RPM*Units.rpm
    Radius          = 1
    n_sections      = 11
    Radial_sections = np.linspace(0, Radius, n_sections)[1:]
    AoA             = np.atleast_2d(np.linspace(-5,10,16))*Units.degrees
    # set of Re for a blade given the Re at the tip
    Re_vals         = np.atleast_1d(Radial_sections*1e7)
    
    
    # airfoil geometry  
    airfoil_geometry_file_1  = rel_path + 'NACA_4412.txt' 
    
    # compute geometry points of two airfoils 
    airfoil_file_1       = '2412'

    airfoil_geometry_1   = compute_naca_4series(airfoil_file_1,npoints = 201) # number of points must always be (a*10^b) + 1
    
    airfoil_properties_1 = airfoil_analysis(airfoil_geometry_1,AoA,Re_vals)  
    plot_airfoil_surface_forces(airfoil_properties_1 )   
    plot_airfoil_boundary_layer_properties(airfoil_properties_1  )    
     
 
    return 
    
 

if __name__ == '__main__': 
    main() 
    plt.show()