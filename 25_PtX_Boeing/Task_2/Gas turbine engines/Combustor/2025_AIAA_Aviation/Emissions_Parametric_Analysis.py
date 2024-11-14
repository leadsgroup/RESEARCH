# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from   RCAIDE.Framework.Core                                                             import Units   
from   RCAIDE.Library.Methods.Emissions.Chemical_Reactor_Network_Method.evaluate_cantera import evaluate_cantera 
from   RCAIDE.Library.Plots                                                              import *     
from   RCAIDE.Library.Components.Propulsors.Converters.Combustor                         import Combustor

# python imports 
import numpy             as np  
import matplotlib.pyplot as plt  
import os   

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
    # Parameters under investigation
    #fuel_types        = ['Jet_A1', 'Methane']
    comb_lengths      = np.array([0.08, 0.1, 0.11])
    comb_areas        = np.array([0.15, 0.3])
    temperatures      = np.array([700, 725, 750, 775, 800])
    pressures         = np.array([20000000, 30000000])
    mdots             = np.array([1, 5, 10])
    FARs              = np.array([0.06, 0.07])     
    
    EI_CO2            = np.zeros((len(comb_lengths), len(comb_areas), len(temperatures), len(pressures), len(mdots), len(FARs)))
    EI_CO             = np.zeros_like(EI_CO2)  
    EI_H2O            = np.zeros_like(EI_CO2)  
    EI_NO             = np.zeros_like(EI_CO2)  
    EI_NO2            = np.zeros_like(EI_CO2)     

    #for f in fuel_types:
    for l in comb_lengths:  
        for a in comb_areas:   
            for t in temperatures:   
                for p in pressures: 
                    for m in mdots:  
                        for far in FARs:      

                            Combustor.A_PZ = comb_areas[a]         # [m**2]    Primary Zone cross-sectional area     
                            Combustor.L_PZ = 0.2*comb_lengths[l]   # [m]       Primary Zone length     
                            Combustor.A_SZ = comb_areas[a]         # [m**2]    Secondary Zone cross-sectional area
                            Combustor.L_SZ = 0.8*comb_lengths[l]   # [m]       Secondary Zone length          
                            T              = temperatures[t]
                            P              = pressures[p]
                            mdot           = mdots[m]
                            FAR            = FARs[far]
                            
                            results = evaluate_cantera(Combustor,T,P,mdot,FAR)
                            
                            EI_CO2[l,a,t,p,m,far] = results.EI_CO2
                            EI_CO[l,a,t,p,m,far]  = results.EI_CO 
                            EI_H2O[l,a,t,p,m,far] = results.EI_H2O
                            EI_NO[l,a,t,p,m,far]  = results.EI_NO 
                            EI_NO2[l,a,t,p,m,far] = results.EI_NO2
    
    return    

if __name__ == '__main__': 
    main()
    plt.show()