# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from   RCAIDE.Framework.Core                                                             import Units   
from   RCAIDE.Library.Methods.Emissions.Chemical_Reactor_Network_Method.evaluate_cantera import evaluate_cantera 
from   RCAIDE.Library.Plots                                                              import *     

# python imports 
import numpy             as np  
import matplotlib.pyplot as plt  
import os   

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
    comb_lengths      = np.linspace(0.05, 0.14, 10)
    comb_areas        = np.linspace(0.15, 0.30, 10)
    temperatures      = np.linspace(700,  800,  10)
    pressures         = np.linspace(20e6, 30e6, 10)
    mdots             = np.linspace(1,    10,   10) 
    FARs              = np.linspace(0.06, 0.07, 10)
    
    EI_CO2            = np.zeros((len(comb_lengths), len(comb_areas), len(temperatures), len(pressures), len(mdots), len(FARs)))
    EI_CO             = np.zeros_like(EI_CO2)  
    EI_H2O            = np.zeros_like(EI_CO2)  
    EI_NO             = np.zeros_like(EI_CO2)  
    EI_NO2            = np.zeros_like(EI_CO2) 
    
    combustor         = RCAIDE.Library.Components.Propulsors.Converters.Combustor()
    combustor.fuel_data = RCAIDE.Library.Attributes.Propellants.Jet_A1()  

    #for f in fuel_types:
    for l in range(len(comb_lengths)):  
        for a in range(len(comb_areas)):   
            for t in range(len(temperatures)):   
                for p in range(len(pressures)): 
                    for m in range(len(mdots)):  
                        for far in range(len(FARs)):      

                            combustor.A_PZ = comb_areas[a]         # [m**2]    Primary Zone cross-sectional area     
                            combustor.L_PZ = 0.2*comb_lengths[l]   # [m]       Primary Zone length     
                            combustor.A_SZ = comb_areas[a]         # [m**2]    Secondary Zone cross-sectional area
                            combustor.L_SZ = 0.8*comb_lengths[l]   # [m]       Secondary Zone length          
                            T              = temperatures[t]
                            P              = pressures[p]
                            mdot           = mdots[m]
                            FAR            = FARs[far]
                            
                            results = evaluate_cantera(combustor,T,P,mdot,FAR)
                            
                            EI_CO2[l,a,t,p,m,far] = results.EI_CO2
                            EI_CO[l,a,t,p,m,far]  = results.EI_CO 
                            EI_H2O[l,a,t,p,m,far] = results.EI_H2O
                            EI_NO[l,a,t,p,m,far]  = results.EI_NO 
                            EI_NO2[l,a,t,p,m,far] = results.EI_NO2
                            
                            
    
    return    

if __name__ == '__main__': 
    main()
    plt.show()