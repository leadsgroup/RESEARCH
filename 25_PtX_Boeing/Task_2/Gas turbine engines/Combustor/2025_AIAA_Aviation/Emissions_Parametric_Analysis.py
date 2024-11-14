# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# RCAIDE imports 
import RCAIDE
from   RCAIDE.Framework.Core                                         import Units   
from   RCAIDE.Library.Methods.Propulsors.Turbofan_Propulsor          import design_turbofan
from   RCAIDE.Library.Methods.Emissions.Chemical_Reactor_Network_Method.evaluate_cantera import evaluate_cantera 
from   RCAIDE.Library.Methods.Geometry.Planform                      import segment_properties
from   RCAIDE.Library.Plots                                          import *     

# python imports 
import numpy as np  
from   copy import deepcopy
import matplotlib.pyplot as plt  
import os   

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():

    fuel_types        = ['Jet_A1', 'Methane']
    comb_lengths      = np.array([0.08, 0.1])
    comb_areas        = np.array([0.15, 0.3])
    temperatures      = np.array([700, 800])
    pressures         = np.array([20000000, 30000000])
    mdots             = np.array([1, 10])
    FARs              = np.array([0.06, 0.07])     
    
    CO2_total = np.zeros(5) 
    CO_total  = 0 * state.ones_row(1) 
    H2O_total = 0 * state.ones_row(1) 
    NO_total  = 0 * state.ones_row(1) 
    NO2_total = 0 * state.ones_row(1) 

    for f in fuel_types:
        for l in comb_lengths:  
            for a in comb_areas:   
                for t in temperatures:   
                    for p in pressures: 
                        for m in mdots:  
                            for far in FARs:      

                                combustor.A_PZ = comb_areas[a]         # [m**2]    Primary Zone cross-sectional area     
                                combustor.L_PZ = 0.2*comb_lengths[l]   # [m]       Primary Zone length     
                                combustor.A_SZ = comb_areas[a]         # [m**2]    Secondary Zone cross-sectional area
                                combustor.L_SZ = 0.8*comb_lengths[l]   # [m]       Secondary Zone length                                      
                            
                                # unpack component conditions
                                n_cp                 = state.numerics.number_of_control_points 
                                propulsor_conditions = state.conditions.energy[fuel_line.tag][propulsor.tag] 
                                combustor_conditions = propulsor_conditions[combustor.tag]  
                                
                                T    = temperatures[t]
                                P    = pressures[p]
                                mdot = mdots[m]
                                FAR  = FARs[far]

                                EI_CO2_comb   = 0 * state.ones_row(1)   
                                EI_CO_comb    = 0 * state.ones_row(1)  
                                EI_H2O_comb   = 0 * state.ones_row(1)  
                                EI_NO_comb    = 0 * state.ones_row(1)  
                                EI_NO2_comb   = 0 * state.ones_row(1)   
                                
                                results = evaluate_cantera(combustor,T,P,mdot,FAR)
                                
                                EI_CO2_comb[t_idx,0] = results.EI_CO2
                                EI_CO_comb[t_idx,0]  = results.EI_CO 
                                EI_H2O_comb[t_idx,0] = results.EI_H2O
                                EI_NO_comb[t_idx,0]  = results.EI_NO 
                                EI_NO2_comb[t_idx,0] = results.EI_NO2
                                
                                EI_CO2_prev = EI_CO2_comb 
                                EI_CO_prev  =  EI_CO_comb  
                                EI_H2O_prev = EI_H2O_comb 
                                EI_NO_prev  =  EI_NO_comb  
                                EI_NO2_prev = EI_NO2_comb 
                                    
                                CO2_total  += np.dot(I,mdot*EI_CO2_comb)
                                CO_total   += np.dot(I,mdot *EI_CO_comb )
                                H2O_total  += np.dot(I,mdot*EI_H2O_comb)
                                NO_total   += np.dot(I,mdot *EI_NO_comb ) 
                                NO2_total  += np.dot(I,mdot *EI_NO2_comb)
    
    emissions                 = Data()
    emissions.total           = Data()
    emissions.index           = Data() 
    emissions.total.CO2       = CO2_total  * combustor.fuel_data.global_warming_potential_100.CO2 
    emissions.total.H2O       = H2O_total  * combustor.fuel_data.global_warming_potential_100.H2O  
    emissions.total.NOx       = (NO_total + NO2_total) * combustor.fuel_data.global_warming_potential_100.NOx 
    emissions.index.CO2       = EI_CO2_comb
    emissions.index.CO        = EI_CO_comb 
    emissions.index.H2O       = EI_H2O_comb
    emissions.index.NO        = EI_NO_comb 
    emissions.index.NO2       = EI_NO2_comb 
    
    state.conditions.emissions =  emissions
    
    return    

if __name__ == '__main__': 
    main()
    plt.show()