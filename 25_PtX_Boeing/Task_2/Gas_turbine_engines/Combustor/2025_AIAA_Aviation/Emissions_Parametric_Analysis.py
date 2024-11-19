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
    
    comb_lengths      = np.linspace(0.05, 0.45,   5)
    comb_areas        = np.linspace(0.05, 0.45,   5)
    temperatures      = np.linspace(800,   900,   5)
    pressures         = np.linspace(4e6,   6e6,   5)
    mdots             = np.linspace(60,     90,   5) 
    FARs              = np.linspace(0.05,  0.7,   5)    
    
    parametric_analyses(comb_lengths, comb_areas, temperatures, pressures, mdots, FARs)         
    
    #emissions_L_A = read_results('emissions_L_A')
    #emissions_T_P = read_results('emissions_T_P')
    #emissions_M_F = read_results('emissions_M_F')
    
    #plot_results(emissions_L_A, emissions_T_P, emissions_M_F, comb_lengths, comb_areas, temperatures, pressures, mdots, FARs)
      
    return    

def parametric_analyses(comb_lengths, comb_areas, temperatures, pressures, mdots, FARs):

    EI_CO2_L_A          = np.zeros((len(comb_lengths), len(comb_areas)))
    EI_CO_L_A           = np.zeros_like(EI_CO2_L_A)  
    EI_H2O_L_A          = np.zeros_like(EI_CO2_L_A)  
    EI_CO2_T_P          = np.zeros((len(temperatures), len(pressures)))
    EI_CO_T_P           = np.zeros_like(EI_CO2_T_P)  
    EI_H2O_T_P          = np.zeros_like(EI_CO2_T_P)
    EI_CO2_M_F          = np.zeros((len(mdots), len(FARs)))
    EI_CO_M_F           = np.zeros_like(EI_CO2_M_F)  
    EI_H2O_M_F          = np.zeros_like(EI_CO2_M_F)    

    combustor           = RCAIDE.Library.Components.Propulsors.Converters.Combustor()
    #combustor.fuel_data = RCAIDE.Library.Attributes.Propellants.Jet_A1()
    combustor.fuel_data = RCAIDE.Library.Attributes.Propellants.Methane()     

    for l in range(len(comb_lengths)):  
        for a in range(len(comb_areas)):   

            combustor.A_PZ = comb_areas[a]         # [m**2]    Primary Zone cross-sectional area     
            combustor.L_PZ = 0.2*comb_lengths[l]   # [m]       Primary Zone length     
            combustor.A_SZ = comb_areas[a]         # [m**2]    Secondary Zone cross-sectional area
            combustor.L_SZ = 0.8*comb_lengths[l]   # [m]       Secondary Zone length          
            T              = 850
            P              = 3000000
            mdot           = 80
            FAR            = 0.0175
            
            print(f"Simulating L =", "%0.3f" % comb_lengths[l] , "m and A =", "%0.3f" % comb_areas[a] , "m^2.")    
            
            try:

                results = evaluate_cantera(combustor, T, P, mdot, FAR)
            
                EI_CO2_L_A[l,a] = results.EI_CO2
                EI_CO_L_A[l,a]  = results.EI_CO 
                EI_H2O_L_A[l,a] = results.EI_H2O
                
            except Exception as e:

                print(f"An error occurred while evaluating Cantera at L =", "%0.3f" % comb_lengths[l] , "m and A =", "%0.3f" % comb_areas[a] , "m^2.")                
            
            emission_indices_L_A = {
            "EI_CO2": EI_CO2_L_A,
            "EI_CO": EI_CO_L_A,
            "EI_H2O": EI_H2O_L_A}
        
            filename = 'emissions_L_A'
            save_results(emission_indices_L_A, filename)
      

    for t in range(len(temperatures)):   
        for p in range(len(pressures)):      

            combustor.A_PZ = 0.15             # [m**2]    Primary Zone cross-sectional area     
            combustor.L_PZ = 0.015            # [m]       Primary Zone length     
            combustor.A_SZ = 0.15             # [m**2]    Secondary Zone cross-sectional area
            combustor.L_SZ = 0.075            # [m]       Secondary Zone length          
            T              = temperatures[t]
            P              = pressures[p]
            mdot           = 80
            FAR            = 0.0175
            
            print(f"Simulating T =", "%0.3f" % temperatures[t] , "K and P =", "%0.3f" % pressures[p] , "Pa")
            
            try:

                results = evaluate_cantera(combustor, T, P, mdot, FAR)
            
                # Assign values to EI arrays
                EI_CO2_T_P[t, p] = results.EI_CO2
                EI_CO_T_P[t, p] = results.EI_CO
                EI_H2O_T_P[t, p] = results.EI_H2O
                
            except Exception as e:

                print(f"An error occurred while evaluating Cantera at T =", "%0.3f" % temperatures[t] , "K and P =", "%0.3f" % pressures[p] , "Pa")            
            
            emission_indices_T_P = {
            "EI_CO2": EI_CO2_T_P,
            "EI_CO": EI_CO_T_P,
            "EI_H2O": EI_H2O_T_P}
        
            filename = 'emissions_T_P'
            save_results(emission_indices_T_P, filename)
    
    
    for m in range(len(mdots)):  
        for far in range(len(FARs)):      
    
            combustor.A_PZ = 0.15    # [m**2]    Primary Zone cross-sectional area     
            combustor.L_PZ = 0.015   # [m]       Primary Zone length     
            combustor.A_SZ = 0.15    # [m**2]    Secondary Zone cross-sectional area
            combustor.L_SZ = 0.075   # [m]       Secondary Zone length          
            T              = 850
            P              = 3000000
            mdot           = mdots[m]
            FAR            = FARs[far] 
            
            print(f"Evaluating mdot_a =", "%0.3f" % mdots[m] , "kg/s and FAR =", "%0.3f" % FARs[far] , ".")
            
            try:

                results = evaluate_cantera(combustor, T, P, mdot, FAR)
            
                EI_CO2_M_F[m,far] = results.EI_CO2
                EI_CO_M_F[m,far]  = results.EI_CO 
                EI_H2O_M_F[m,far] = results.EI_H2O
                
            except Exception as e:

                print(f"An error occurred while evaluating Cantera at mdot_a =", "%0.3f" % mdots[m] , "kg/s and FAR =", "%0.3f" % FARs[far] , ".")                 
            
            emission_indices_M_F = {
            "EI_CO2": EI_CO2_M_F,
            "EI_CO": EI_CO_M_F,
            "EI_H2O": EI_H2O_M_F}
        
            filename = 'emissions_M_F'
            save_results(emission_indices_M_F, filename)

    return  




def read_results(file_name):
    emissions = load_results(file_name)
    print(emissions)
    return emissions 


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

def plot_results(emissions_L_A, emissions_T_P, emissions_M_F, comb_lengths, comb_areas, temperatures, pressures, mdots, FARs):
    
    plot_style()
    
    data = emissions_L_A['EI_CO2']
    
    comb_diameters    = 2 * np.sqrt(comb_areas / np.pi)
    X, Y              = np.meshgrid(comb_lengths, comb_diameters)
    
    fig    =  plt.figure(1)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_CO2')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Length [m]')
    plt.ylabel('Diameter [m]')
    
    data = emissions_L_A['EI_CO']
    
    fig    =  plt.figure(2)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_CO')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Length [m]')
    plt.ylabel('Diameter [m]')
    
    data = emissions_L_A['EI_H2O']
    
    fig    =  plt.figure(3)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_H2O')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Length [m]')
    plt.ylabel('Diameter [m]')    
    
    # ---------------------------------------------------------
    
    data = emissions_T_P['EI_CO2']
    
    X, Y              = np.meshgrid(temperatures, pressures)
    
    fig    =  plt.figure(4)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_CO2')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Temperature [K]')
    plt.ylabel('Pressure [Pa]')
    
    data = emissions_T_P['EI_CO']
    
    fig    =  plt.figure(5)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_CO')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Temperature [K]')
    plt.ylabel('Pressure [Pa]')
    
    data = emissions_T_P['EI_H2O']
    
    fig    =  plt.figure(6)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_H2O')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Temperature [K]')
    plt.ylabel('Pressure [Pa]')  
    
    # ---------------------------------------------------------
    
    data = emissions_M_F['EI_CO2']
    
    X, Y              = np.meshgrid(mdots, FARs)
    
    fig    =  plt.figure(7)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_CO2')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Air mass flow rate [kg/s]')
    plt.ylabel('Fuel to Air ratio [-]')
    
    data = emissions_M_F['EI_CO']
    
    fig    =  plt.figure(8)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_CO')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Air mass flow rate [kg/s]')
    plt.ylabel('Fuel to Air ratio [-]')
    
    data = emissions_M_F['EI_H2O']
    
    fig    =  plt.figure(9)
    contour = plt.contourf(X, Y, data, levels=50, cmap='viridis')
    color_bar = plt.colorbar(contour, label='EI_H2O')
    color_bar.formatter.useOffset = False
    color_bar.update_ticks()   
    plt.xlabel('Air mass flow rate [kg/s]')
    plt.ylabel('Fuel to Air ratio [-]')       
    
    return

def plot_style():
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {
        'axes.labelsize': 20,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'axes.titlesize': 18,
        'figure.dpi': 130
    }
    plt.rcParams.update(parameters)


if __name__ == '__main__': 
    main()
    plt.show()