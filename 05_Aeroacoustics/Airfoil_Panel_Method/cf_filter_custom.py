## @ingroup Methods-Aerodynamics-Airfoil_Panel_Method  
# RCAIDE/Methods/Aerodynamics/Airfoil_Panel_Method/cf_filter.py
# 
# 
# Created:  Jun 2024, Niranjan Nanjappa 

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------    

# RCAIDE imports  
from RCAIDE.Framework.Core import Data

# package imports  
import numpy as np 
from scipy.signal import lfilter 

# ----------------------------------------------------------------------------------------------------------------------
# cf_filter.py
# ----------------------------------------------------------------------------------------------------------------------
## @ingroup Methods-Aerodynamics-Airfoil_Panel_Method
def cf_filter_custom(ncpts,ncases,npanel,CF):
    
    n = npanel
    
    # number of points to be excluded from the filter
    n_f = n/10
    n_r = n/20
    CF_new = CF
    
    for i in range(ncpts):
        for j in range(ncases):
            CF_diff = np.diff(CF_new[:,j,i])
            for k in range(int(n_r),int(n/2-n_f)):
                if CF_diff[k]*CF_diff[k+1]<0:
                    CF_new[k,j,i] = 0.5*(CF_new[k-1,j,i]+CF_new[k+1,j,i])
    
    return CF_new