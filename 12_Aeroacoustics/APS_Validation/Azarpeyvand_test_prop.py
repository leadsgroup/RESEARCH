# Azarpeyvand_test_prop.py
#
# Created:  Nov 2024, Niranjan Nanjappa

# Imports
import RCAIDE
from RCAIDE.Framework.Core import Units, Data  
from RCAIDE.Library.Methods.Geometry.Airfoil  import compute_airfoil_properties
from RCAIDE.Library.Methods.Geometry.Airfoil  import import_airfoil_geometry    
from scipy.interpolate import interp1d
import os
import numpy as np   

# design propeller                                       
def Azarpeyvand_test_prop():  
    prop                            = RCAIDE.Library.Components.Propulsors.Converters.Propeller()
    prop.inputs                     = Data()
    prop.inputs.pitch_command       = 0 
    prop.inputs.y_axis_rotation     = 0.
    prop.tag                        = 'Azarpeyvand_test_prop'  
    prop.tip_radius                 = 0.3048/2
    prop.hub_radius                 = prop.tip_radius*0.113139456717213
    prop.number_of_blades           = 2  
    r_R_data                        = np.array([0.113139457,0.16471,0.216280543,0.267851086,0.319421629,0.370992172,0.422562715,0.474133258,0.525703801,
                                                0.577274344,0.628844887,0.68041543,0.731985973,0.783556517,0.83512706,0.886697603,0.938268146,0.989838689])    
    t_c_data                        = np.ones(len(r_R_data))*0.117    
    c_data                          = np.array([16.45519126,20.45882952,21.08246808,21.4098497,21.79997386,22.04670533,22.08708736,21.89358428,21.4101269,
                                                20.71029621,19.8310804,18.74263889,17.48978606,16.08010468,14.4529353,12.58488109,10.49129543,6.878770192])*(1e-3)
    beta_data                       = np.array([52.27770172,52.21064838,44.32572311,36.91881267,31.00590752,26.37881448,22.69919181,19.93749725,17.82350885,
                                                16.18482705,14.918502,13.70006443,12.77316622,11.84349977,11.1099226,10.39560036,9.766041674,11.02957779])
    num_sec = 30          
    new_radius_distribution         = np.linspace(0.113139457,0.989838689 ,num_sec)
    func_twist_distribution         = interp1d(r_R_data, (beta_data)*Units.degrees , kind='cubic')
    func_chord_distribution         = interp1d(r_R_data, c_data , kind='cubic')
    func_radius_distribution        = interp1d(r_R_data, r_R_data*prop.tip_radius , kind='cubic')
    func_max_thickness_distribution = interp1d(r_R_data, t_c_data*c_data, kind='cubic')  
    prop.twist_distribution         = func_twist_distribution(new_radius_distribution)     
    prop.chord_distribution         = func_chord_distribution(new_radius_distribution)         
    prop.radius_distribution        = func_radius_distribution(new_radius_distribution)        
    prop.max_thickness_distribution = func_max_thickness_distribution(new_radius_distribution) 
    prop.thickness_to_chord         = prop.max_thickness_distribution/prop.chord_distribution 
    ospath    = os.path.abspath(__file__)
    separator = os.path.sep
    rel_path  = os.path.dirname(ospath) + separator  
    airfoil                          = RCAIDE.Library.Components.Airfoils.Airfoil()   
    airfoil.coordinate_file          = rel_path +'/Airfoils/Clark_y.txt'
    airfoil.polar_files              = [rel_path +'/Airfoils/Polars/Clark_y_polar_Re_50000.txt' ,
                                       rel_path +'/Airfoils/Polars/Clark_y_polar_Re_100000.txt',
                                       rel_path +'/Airfoils/Polars/Clark_y_polar_Re_200000.txt',
                                       rel_path +'/Airfoils/Polars/Clark_y_polar_Re_500000.txt',
                                       rel_path +'/Airfoils/Polars/Clark_y_polar_Re_1000000.txt']
    airfoil.geometry                 = import_airfoil_geometry(airfoil.coordinate_file,airfoil.number_of_points)
    airfoil.polars                   = compute_airfoil_properties(airfoil.geometry,airfoil.polar_files)
    prop.append_airfoil(airfoil) 
    prop.airfoil_polar_stations      = list(np.zeros(num_sec).astype(int))  
    prop.mid_chord_alignment         = np.zeros_like(prop.chord_distribution)  
        
    return prop


prop = Azarpeyvand_test_prop()