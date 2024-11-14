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
import pickle

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():
    
 
    #parametric_analyses()
    read_results('Emissions')
 
    return    

def parametric_analyses():
    comb_lengths      = np.linspace(0.05, 0.14, 1)
    comb_areas        = np.linspace(0.15, 0.30, 1)
    temperatures      = np.linspace(700,  800,  1)
    pressures         = np.linspace(20e6, 30e6, 1)
    mdots             = np.linspace(1,    10,   1) 
    FARs              = np.linspace(0.06, 0.07, 1)

    EI_CO2            = np.zeros((len(comb_lengths), len(comb_areas), len(temperatures), len(pressures), len(mdots), len(FARs)))
    EI_CO             = np.zeros_like(EI_CO2)  
    EI_H2O            = np.zeros_like(EI_CO2)  
    #EI_NO             = np.zeros_like(EI_CO2)  
    #EI_NO2            = np.zeros_like(EI_CO2) 

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
                            #EI_NO[l,a,t,p,m,far]  = results.EI_NO 
                            #EI_NO2[l,a,t,p,m,far] = results.EI_NO2
                            
                            
    #emission_indices  = [EI_CO2,EI_CO,EI_H2O,EI_NO,EI_NO2]  
    emission_indices = {
    "EI_CO2": EI_CO2,
    "EI_CO": EI_CO,
    "EI_H2O": EI_H2O}

    filename = '=emissions'
    save_results(emission_indices, filename)


    return     
def read_results(file_name):
    emissions = load_results(file_name)
    print(emissions)
    return  


# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def save_results(results,filename):
   #  Pickle Backup Files
    current_dir = os.path.dirname(os.path.abspath(__file__))
    load_dir = os.path.join(current_dir)
    pickle_file = os.path.join(load_dir, filename + '.pkl')
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return

# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    current_dir = os.path.dirname(os.path.abspath(__file__))
    load_dir = os.path.join(current_dir)
    load_file = os.path.join(load_dir, filename+ '.pkl')
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results



if __name__ == '__main__': 
    main()
    plt.show()