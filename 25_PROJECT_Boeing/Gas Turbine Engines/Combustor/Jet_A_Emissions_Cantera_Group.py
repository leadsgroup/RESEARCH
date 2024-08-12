import cantera as ct
import numpy as np
import pandas as pd

#--------------------------------------------------------------------------------
# initial conditions

temp = 600 # K
patm = 25 # atm
equ_ratio = 0.7

list_tpfr = np.linspace(0.01,10,10)*1e-4 # psr residence time in sec
tpsr = 1e-3 

#dict_fuel = {'N-C12H26':0.6, 'A1CH3':0.2, 'A1':0.2}
dict_fuel = {'NC10H22':0.16449, 'NC12H26':0.34308, 'NC16H34':0.10335, 'IC8H18':0.08630, 'NC7H14':0.07945, 'C6H5C2H5': 0.07348, 'C6H5C4H9': 0.05812, 'C10H7CH3': 0.10972}      # [2] More accurate kinetic mechanism, slower simulation    
dict_oxy = {'O2':0.2095, 'N2':0.7809, 'AR':0.0093, 'CO2':0.0003}

#gas = ct.Solution('JetFuelSurrogate.yaml')
gas = ct.Solution('chem.yaml')

#--------------------------------------------------------------------------------

def combustor(tpfr):
    
    """ combustor simulation using a simple psr-pfr reactor network with varying pfr residence time """

    gas.TP = temp, patm*ct.one_atm
    gas.set_equivalence_ratio(equ_ratio, fuel = dict_fuel, oxidizer = dict_oxy )
        
    comp_fuel = list(dict_fuel.keys())
    Y_fuel = gas[comp_fuel].Y
    
    # psr (flame zone)
    
    upstream = ct.Reservoir(gas)
    downstream = ct.Reservoir(gas)
        
    gas.equilibrate('HP')
    psr = ct.IdealGasReactor(gas)
    
    func_mdot = lambda t: psr.mass/tpsr
    
    inlet = ct.MassFlowController(upstream, psr)
    inlet.mass_flow_rate = func_mdot
    outlet = ct.Valve(psr, downstream, K=100)
            
    sim_psr = ct.ReactorNet([psr])
        
    try:
        sim_psr.advance_to_steady_state()
    except RuntimeError:
        pass
    
    # pfr (burn-out zone)
    
    pfr = ct.IdealGasConstPressureReactor(gas)
    sim_pfr = ct.ReactorNet([pfr])
    
    try:
        sim_pfr.advance(tpfr)
    except RuntimeError:
        pass
        
    mdot = inlet.mass_flow_rate
    mdot_fuel = sum(mdot * Y_fuel)
    Emission_Index = gas.Y * mdot/mdot_fuel
    
    return (gas, Emission_Index)
            
#--------------------------------------------------------------------------------

list_sp = ['CO', 'CO2', 'H2O', 'NO', 'NO2', 'CSOLID']
col_names = ['tau(s)', 'Tout(K)'] + ['X_' +str(sp) for sp in list_sp] + ['Y_' +str(sp) for sp in list_sp] + ['EI_' +str(sp) for sp in list_sp]
df = pd.DataFrame(columns=col_names)

for (n, tau) in enumerate(list_tpfr):
    gas, EI = combustor(tau)
    sp_idx = [gas.species_index(sp) for sp in list_sp]
    data_n = [tau, gas.T] + list(gas.X[sp_idx]) + list(gas.Y[sp_idx]) + list(EI[sp_idx])
    df.loc[n] = data_n

df.to_csv('output.csv')    