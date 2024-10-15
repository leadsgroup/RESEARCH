# Geopy comparison test.
# This serves to validate that the Geopy implementation agrees with the original library
#
# Created Oct. 2024, A. Molloy
# -------------------------------------------------------------------------

import  numpy as  np
from geopy.distance import geodesic as GD
from  RCAIDE.Framework.Analyses.Geodesics.Geodesics import Calculate_Distance

def  main():
    Geopy_Comparison()

def  Geopy_Comparison():
    lat_range = -40 + 3 / 20 * np.array(range(0, 20))
    lon_range = -69 + 3 / 20 * np.array(range(0, 20))
    
    for i in  range(len(lat_range)):
        for j in  range(len(lon_range)):
            test_point1 = np.array([lat_range[i], lon_range[j]])
            for k in  range(len(lat_range)):
                for l in  range(len(lon_range)):
                    test_point2 = np.array([lat_range[k], lon_range[l]])
                    True_result = GD(test_point1, test_point2).km
                    Test_result = Calculate_Distance(test_point1, test_point2)
                    error = True_result - float(Test_result)
                    if abs(error) > 0.00000001:
                        print("Error of "+str(error)+" km with points: "+str(test_point1)+" "+str(test_point2))
    print("Last point: "+str(test_point1)+" "+str(test_point2))

if __name__ == '__main__': 
    #with jax.disable_jit():
            
    main() 