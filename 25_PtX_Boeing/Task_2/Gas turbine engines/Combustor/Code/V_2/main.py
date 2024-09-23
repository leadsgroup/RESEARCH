import pandas               as pd
import matplotlib.pyplot    as plt
import time
import json

from src import Combustor

'''
This code is built to estimate the Emission Index of various species deriving
from the combustion of Jet-A fuel. The model is based on a Chemical Reactor
Network (CRN) built using Cantera.
'''

# Created:  Sep 2024, M. Guidotti

# References:
# [1]: Allaire, Douglas & Astronautics, Massachusetts. (2007). A physics-based 
#      emissions model for aircraft gas turbine combustors.
# [2]: Brink, L. F. J. (2020). Modeling the impact of fuel composition on aircraft 
#      engine NOx, CO and soot emissions. Masters thesis, Massachusetts Institute of
#      Technology.

def main():

    ti                      = time.time()                                           # [s]       Define the initial simulation time

    # ------------------------------------------------------------------------------
    # ------------------------------ Combustor Inputs ------------------------------
    # ------------------------------------------------------------------------------

    with open('combustor_input.json', 'r') as file:
        combustor_input                  = json.load(file)                                           # [-]       Import input parameters for the combustor
    
    combustor_input['m_dot_fuel_tot']    = combustor_input['m_dot_air_tot']*combustor_input['FAR']                                     # [kg/s]    Fuel mass flow going through all combustors
    combustor_input['m_dot_air']         = combustor_input['m_dot_air_tot']/combustor_input['N_comb']                                  # [kg/s]    Air mass flow inside each combustor, scaled inside each PSR to vary the Equivalence Ratio
    combustor_input['m_dot_fuel']        = combustor_input['m_dot_fuel_tot']/combustor_input['N_comb']      

    if combustor_input['high_fidelity_kin_mech']:
        with open('fuel_input_1.json', 'r') as file:
            combustor_input['dict_fuel'] = json.load(file)                                             # [-]       Fuel species and corresponding mole fractions for full fuel model
    else:
        with open('fuel_input_2.json', 'r') as file:
            combustor_input['dict_fuel'] = json.load(file)                                             # [-]       Fuel species and corresponding mole fractions for surrogate fuel model
    
    with open('air_input.json', 'r') as file:
        combustor_input['dict_oxy'] = json.load(file)                                                  # [-]       Air species and corresponding mole fractions    
    
    col_names               = ['Tout(K)', 'T_stag_out','P_stag_out', 'h_stag_out'] + ['X_' +str(sp) for sp in combustor_input['dict_fuel']['list_sp']] + ['Y_' +str(sp) for sp in combustor_input['dict_fuel']['list_sp']] + ['EI_' +str(sp) for sp in combustor_input['dict_fuel']['list_sp']] # [-]       Define output variables
    output                  = pd.DataFrame(columns=col_names)                                       # [-]       Assign output variables space to output

    for n in range(1):
        gas, EI, T_stag_out, P_stag_out, h_stag_out = Combustor(combustor_input) # [-]       Run combustor function

        sp_idx              = [gas.species_index(sp) for sp in combustor_input.list_sp]             # [-]       Retrieve the species index
        data_n              = [gas.T, T_stag_out, P_stag_out, h_stag_out] + list(gas.X[sp_idx]) + list(gas.Y[sp_idx]) + list(EI[sp_idx]) # [-]       Assign output variables
        output.loc[n]       = data_n                                                # [-]       Assign output variables to output

    print(output['EI_CO2'])                                                             # [-]       Print the value of EI_CO2
    print(output['EI_CO'])                                                              # [-]       Print the value of EI_CO
    print(output['EI_H2O'])                                                             # [-]       Print the value of EI_H2O

    if combustor_input['high_fidelity_kin_mech']:
        print(output['EI_NO'])                                                          # [-]       Print the value of EI_NO
        print(output['EI_NO2'])                                                         # [-]       Print the value of EI_NO2
        print(output['EI_CSOLID'])                                                      # [-]       Print the value of EI_CSOLID

    tf                      = time.time()                                           # [s]       Define the final simulation time
    output.elapsed_time     = round((tf-ti),2)                                      # [s]       Compute the total simulation time

    print('Simulation Time: ' + str(output.elapsed_time) + ' seconds per timestep')        # [-]       Print the value of total simulation time

    return

if __name__ == '__main__':
    main()
    plt.show()