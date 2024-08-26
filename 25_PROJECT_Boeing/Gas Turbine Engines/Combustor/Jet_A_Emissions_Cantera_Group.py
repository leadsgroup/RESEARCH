import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time 

def main():

    ti = time.time()    
        
    # initial conditions 
    temp      = np.array([600,600]) # K
    patm      = np.array([25,25]) # atm
    equ_ratio = np.array([0.7, 0.7]) 
    tpfr      = np.array([5 ,10]) * 1E-4 # np.linspace(1,10,10)*1e-4 # psr residence time in sec
    tpsr      = np.array([5e-3,5e-3])
    
    dict_fuel = {'N-C12H26':0.6, 'A1CH3':0.2, 'A1':0.2}
    #dict_fuel = {'NC10H22':0.16449, 'NC12H26':0.34308, 'NC16H34':0.10335, 'IC8H18':0.08630, 'NC7H14':0.07945, 'C6H5C2H5': 0.07348, 'C6H5C4H9': 0.05812, 'C10H7CH3': 0.10972}      # [2] More accurate kinetic mechanism, slower simulation    
    dict_oxy = {'O2':0.2095, 'N2':0.7809, 'AR':0.0093, 'CO2':0.0003}
    
    gas = ct.Solution('JetFuelSurrogate.yaml')
    #gas = ct.Solution('chem.yaml')  
    #-------------------------------------------------------------------------------- 
    list_sp = ['CO', 'CO2', 'H2O']
    #list_sp   = ['CO', 'CO2', 'H2O', 'NO', 'NO2', 'CSOLID']
    col_names = ['tau(s)', 'Tout(K)', 'EqRatio','Ent', 'Tstag', 'Pstag'] + ['X_' +str(sp) for sp in list_sp] + ['Y_' +str(sp) for sp in list_sp] + ['EI_' +str(sp) for sp in list_sp]
    df        = pd.DataFrame(columns=col_names)
    
    for n in range(len(tpfr)):
        gas, EI, phi, h, T_stag, P_stag = combustor(tpfr[n],temp[n],patm[n],equ_ratio[n],tpsr[n],dict_fuel,dict_oxy,gas)
        sp_idx = [gas.species_index(sp) for sp in list_sp]
        data_n = [tpfr[n], gas.T, phi, h, T_stag, P_stag] + list(gas.X[sp_idx]) + list(gas.Y[sp_idx]) + list(EI[sp_idx])
        df.loc[n] = data_n
    
    tf           = time.time()
    elapsed_time = round((tf-ti)/len(tpfr),2)
    print('Simulation Time: ' + str(elapsed_time) + ' seconds per timestep')   
        
    plot_emission(df,gas,equ_ratio)
    
    return 
 
def combustor(tau,temp,patm,equ_ratio,tpsr,dict_fuel, dict_oxy, gas):
    
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
    pfr     = ct.IdealGasConstPressureReactor(gas)
    sim_pfr = ct.ReactorNet([pfr])
    
    try:
        sim_pfr.advance(tau) # sim_pfr.advance(tpfr)
    except RuntimeError:
        pass
        
    mdot           = inlet.mass_flow_rate
    mdot_fuel      = sum(mdot * Y_fuel)
    Emission_Index = gas.Y * mdot/mdot_fuel
    
    density_out = pfr.thermo.density  # density of the gas in the combustor
    area_out = 1  # Assuming the area is 1 m^2 for simplification
    vel = mdot / (density_out * area_out)   
    a = gas.sound_speed
    M = vel/a
    cp = gas.cp_mass
    cv = gas.cv_mass
    gamma = cp/cv
    phi = gas.equivalence_ratio(fuel = dict_fuel, oxidizer = dict_oxy )
    h = gas.h
    T_stag = gas.T/(1 + ((gamma - 1)/2)*M**2)
    P_stag = gas.P*(T_stag/gas.T)**(gamma/(gamma - 1))
    
    return (gas, Emission_Index, phi, h, T_stag, P_stag) 

    
def plot_emission(df,gas,equ_ratio): 
    # Plot results
    f, ax1 = plt.subplots(3, 1, figsize=(16, 12))
    f.suptitle('Jet-A EI')
    subtitle = f'Equivalence ratio: {equ_ratio[0]}, Temperature: {gas.T:.1f} K, Pressure: {gas.P/ct.one_atm:.1f} atm,'
    plt.figtext(0.5, 0.925, subtitle, ha='center', fontsize=12)
    ax1[0].plot(df['tau(s)'], df['EI_CO2'], '.-', color='C0')
    ax1[0].axhline(y=3.16, color='r', linestyle='--')
    ax1[0].annotate('Typical EI value: 3.16', xy=(0.5, 3.16), xytext=(0.5, 3.14), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
    ax1[0].set_title('Emission Index CO2', color='C0')
    ax1[0].set_ylabel('EI [kg/kg]')
    ax1[1].plot(df['tau(s)'], df['EI_CO'], '.-', color='C1')
    ax1[1].axhline(y=0.05, color='r', linestyle='--')                                                                                                                 # https://ntrs.nasa.gov/api/citations/19750007129/downloads/19750007129.pdf, https://pdf.sciencedirectassets.com/271798/1-s2.0-S1352231015X00148/1-s2.0-S1352231015301722/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEFQaCXVzLWVhc3QtMSJHMEUCIQCbJLgpzlmLNnfVH4WG%2FuFB%2FqOu8xrYTmZDE%2FFdwY9CRQIgehzVJdCWAPcetTxMSNhAjOt304L1xdH5qkaf5ZtBqVQqswUITRAFGgwwNTkwMDM1NDY4NjUiDFPyX8EpoZkVNqV3lSqQBWYvEJUcWvUdCATtDtCfyxKBy4jWz2XQVw7nsF9S4bmfJqHSb4XTUr5b1SQzBkqZusEjar1Qi65i0CxHA8FSi2T4hdDbv399n2x9e3fLe%2FaI3MBt7ln%2BzuRZTJ7WH18TSE7v1I5rsejXhUloB3jFtkAn%2FWVF2o7uJnrpGPJgsIxMhYMhDk6OZb4FM%2BnhS%2Fb1FE2Q275EvoF2Dvce%2FlJ5srqnH38xoC2CoxxEWmQEDf9dZcQdaANorm2HkZCb2Wfe2o3tFnJruPAgxBOWMrSrt1B01CP3nfahP0R4svM03hcbnG0AFPkZs4ASaPE6cuJgIM0atyIF6QOndgrfRz5Itq5hSixBdjDtwy8AkDPbeDIT%2FTChYcRvJug46fRys56a9gB8w2IIKR5PbGreEDqD76O%2FSFm2%2FthB%2BZPEBGx1VGi7BPtHx%2FCCMxGPQToqx7M0XdDGWB3CLP7ae81xPu%2BtQo5B%2BnUPSVUzCExpO%2Fy6yISXRVeXOGFdbKUUAXFvn6EG0HPnRsPGke7B59RGvcycjRNT77rUD71wsBrrQ4sTmPpjLPrQsgg03z9a8QqWRNuUx4nkGyBCTjaxCCyeH81ZDZcTO1wpOya6AGbrjh61%2Fztmo2D8LiCkhR%2F9%2BuYkSusLq4JhFs%2BPh5%2BUtPA%2BwzZQGmKNCByIsJrvL5dEU8NeSEkwMp2Lc9AcHz8LZFi4FfXwUpnE9N%2FIdEmlUvjPVtmvdi8RICk6sLcOxMmHdSWUvS6aH1qegi1b5p2aul%2BieilyIVYcw533DSA5e9QxF6CXX%2B9lnwFRn2ijs1FLdrE1cBMxLM3hUqmUpcQfAYoX%2BHIxzsiCJjcFzEZmUAvQPO0RCdy1tbv9sElfG%2F1j35eQ3%2BeNMO%2Fj2bUGOrEBlJLvANktdqzp1uCuSH0eCO%2BAH7r5OuiXzQkXuVQRIgNZh7BOpRQUe17BOSmJ8zLSPk2%2BeXQboT8VMd4A6Mau1gAAEd9KCR4LnvOtxC0LbJzih3f0M%2FMaiygaPuRBP5uQ2BcDX4Sl%2Bq7YMu%2F0cLhr4%2BogeEEXJCHeqQy2aU8sxyJhgoc01%2FK4x7DrFG71HUQ6RfC5YTNFBn%2FDvs67d3KmUiuyhnNYzZTh7aGE8pzGonoe&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240809T210310Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY67WYP4PM%2F20240809%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=49f95b8789decc69cd2507fffa9b06d7f74a8dbcb638af573e4e905e76b5923a&hash=77bdba4d533e0c7c4359d88b85c90a90fb1a7e6fc1ef817544435a0d2e984525&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1352231015301722&tid=spdf-0c1e36cc-69f9-44bf-b0e4-bae0ff444eb9&sid=21a47c3f56535040da8be6c1225a593131dfgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=17155c0557070d5703&rr=8b0ab0fc59411cde&cc=us
    ax1[1].annotate('Typical EI value: 0.05', xy=(0.5, 0.05), xytext=(0.5, 0.03), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->')) 
    ax1[1].set_title('Emission Index CO', color='C1')
    ax1[1].set_ylabel('EI [kg/kg]')
    ax1[2].plot(df['tau(s)'], df['EI_H2O'], '.-', color='C2')
    ax1[2].set_xlabel('PFR residence time [s]')
    ax1[2].axhline(y=1.34, color='r', linestyle='--')
    ax1[2].annotate('Typical EI value: 1.34', xy=(0.5, 1.34), xytext=(0.5, 1.32), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
    ax1[2].set_title('Emission Index H2O', color='C2')
    ax1[2].set_ylabel('EI [kg/kg]')
    
    #f, ax1 = plt.subplots(3, 1, figsize=(16, 12))
    #f.suptitle('Jet-A EI')
    #subtitle = f'Equivalence ratio: {equ_ratio[0]}, Temperature: {gas.T:.1f} K, Pressure: {gas.P/ct.one_atm:.1f} atm,'
    #plt.figtext(0.5, 0.925, subtitle, ha='center', fontsize=12)
    #ax1[0].plot(df['tau(s)'], df['EI_NO2'], '.-', color='C0')
    #ax1[0].axhline(y=0.01, color='r', linestyle='--')
    #ax1[0].annotate('Typical EI value: 0.01', xy=(0.5, 0.01), xytext=(0.5, 0.008), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
    #ax1[0].set_title('Emission Index NO2', color='C0')
    #ax1[0].set_ylabel('EI [kg/kg]')
    #ax1[1].plot(df['tau(s)'], df['EI_NO'], '.-', color='C1')
    #ax1[1].axhline(y=0.01, color='r', linestyle='--')
    #ax1[1].annotate('Typical EI value: 0.01', xy=(0.5, 0.01), xytext=(0.5, 0.008), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
    #ax1[1].set_title('Emission Index NO', color='C1')
    #ax1[1].set_ylabel('EI [kg/kg]')
    #ax1[2].plot(df['tau(s)'], df['EI_CSOLID'], '.-', color='C2')
    #ax1[2].axhline(y=0.00004, color='r', linestyle='--')
    #ax1[2].annotate('Typical EI value: 0.00004', xy=(0.5, 0.00004), xytext=(0.5, -0.00196), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
    #ax1[2].set_xlabel('PFR residence time [s]')
    #ax1[2].set_title('Emission Index C-soot', color='C2')
    #ax1[2].set_ylabel('EI [kg/kg]')
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['X_CO2'], '.-', label='CO2 Mole Fraction', color='C0')
    plt.plot(df['tau(s)'], df['X_CO'], '.-', label='CO Mole Fraction', color='C1')
    plt.plot(df['tau(s)'], df['X_H2O'], '.-', label='H2O Mole Fraction', color='C2')
    #plt.plot(df['tau(s)'], df['X_NO2'], '.-', label='NO2 Mole Fraction', color='C3')
    #plt.plot(df['tau(s)'], df['X_NO'], '.-', label='NO Mole Fraction', color='C4')
    #plt.plot(df['tau(s)'], df['X_CSOLID'], '.-', label='Soot Mole Fraction', color='C5')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('Mole Fraction')
    plt.title('Mole Fraction of CO2, CO, H2O, NO2, NO and soot vs. PFR residence time')
    plt.legend()
    plt.grid(True)
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['Y_CO2'], '.-', label='CO2 Mass Fraction', color='C0')
    plt.plot(df['tau(s)'], df['Y_CO'], '.-', label='CO Mass Fraction', color='C1')
    plt.plot(df['tau(s)'], df['Y_H2O'], '.-', label='H2O Mass Fraction', color='C2')
    #plt.plot(df['tau(s)'], df['Y_NO2'], '.-', label='NO2 Mass Fraction', color='C3')
    #plt.plot(df['tau(s)'], df['Y_NO'], '.-', label='NO Mass Fraction', color='C4')
    #plt.plot(df['tau(s)'], df['Y_CSOLID'], '.-', label='Soot Mass Fraction', color='C5')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('Mass Fraction')
    plt.title('Mass Fraction of CO2, CO, H2O, NO2, NO and soot vs. PFR residence time')
    plt.legend()
    plt.grid(True)
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['Tout(K)'], '.-', label='Temperature', color='C0')
    plt.plot(df['tau(s)'], df['Tstag'], '.-', label='Stagnation Temperature', color='C1')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('Stagnation Temperature [K]')
    plt.title('Stagnation Temperature vs. PFR residence time')
    plt.legend()
    plt.grid(True)    
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['Pstag'], '.-', label='Stagnation Temperature', color='C0')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('Stagnation Pressure [Pa]')
    plt.title('Stagnation Pressure vs. PFR residence time')
    plt.legend()
    plt.grid(True) 
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['Ent'], '.-', label='Stagnation Temperature', color='C0')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('Enthalpy [J]')
    plt.title('Enthalpy vs. PFR residence time')
    plt.legend()
    plt.grid(True)  
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['EqRatio'], '.-', label='Stagnation Temperature', color='C0')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('Equivalence Ratio [-]')
    plt.title('Equivalence Ratio vs. PFR residence time')
    plt.legend()
    plt.grid(True)      
    
    plt.show() 
 
    return  

if __name__ == '__main__': 
    main()
    plt.show()