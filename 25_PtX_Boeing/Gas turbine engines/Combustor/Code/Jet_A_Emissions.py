import cantera           as ct
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
import time 

def main():

    ti        = time.time()
    
    #--------------------------------------------------------------------------------

    # Engine inputs 
    T_stag_0                = 700                                       # [K]  
    P_stag_0                = 2000000                                   # [Pa] 
    FAR                     = 0.02                                      # [-]
    FAR_TO                  = 0.0275                                    # [-]
    FAR_st                  = 0.068                                     # [-]
    mdot_fuel               = 0.8                                       # [kg/s]
    mdot_fuel_TO            = 1.2                                       # [kg/s]
    mdot_air                = mdot_fuel/FAR                             # [kg/s]
    mdot_air_TO             = mdot_fuel_TO/FAR_TO                       # [kg/s]     
    
    # Primary Zone inputs
    N_PZ                    = 21                                        # [-]
    V_PZ                    = 0.0023                                    # [m**3]           
    phi_PZ_des              = 1.77                                      # [-]
    S_PZ                    = 0.39                                      # [-]           
    LHV_input_fuel          = 1 #FIX
    LHV_model_fuel          = 1 #FIX
    F_SC                    = LHV_input_fuel/LHV_model_fuel             # [-]
    
    #Secondary Zone inputs
    A_SZ                    = 0.15                                      # [m**2]
    L_SZ                    = 0.075                                     # [m] 
    l_SA_SM                 = 0.55                                      # [-]          
    l_SA_FM                 = 0.055                                     # [-]    
    l_DA_start              = 0.95                                      # [-]          
    l_DA_end                = 0.1                                       # [-] 
    f_SM                    = 0.5                                       # [-]     
    phi_SZ_des              = 0.7                                       # [-]    
    
    dict_fuel = {'N-C12H26':0.6, 'A1CH3':0.2, 'A1':0.2}
    dict_oxy  = {'O2':0.2095,    'N2':0.7809, 'AR':0.0096}    

    #--------------------------------------------------------------------------------
    
    list_sp   = ['CO2', 'CO', 'H2O']
    col_names = ['Tout(K)', 'T_stag_out','P_stag_out', 'h_stag_out', 'FAR'] + ['X_' +str(sp) for sp in list_sp] + ['Y_' +str(sp) for sp in list_sp] + ['EI_' +str(sp) for sp in list_sp]
    df        = pd.DataFrame(columns=col_names)
    
    for n in range(1):
        gas, EI, T_stag_out, P_stag_out, h_stag_out, FAR = combustor(dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, mdot_fuel, mdot_fuel_TO, mdot_air, mdot_air_TO, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des, l_SA_SM, l_SA_FM, F_SC, A_SZ, L_SZ, l_DA_start, l_DA_end, f_SM)
        sp_idx = [gas.species_index(sp) for sp in list_sp]
        data_n = [gas.T, T_stag_out, P_stag_out, h_stag_out, FAR] + list(gas.X[sp_idx]) + list(gas.Y[sp_idx]) + list(EI[sp_idx])
        df.loc[n] = data_n
    
    tf           = time.time()
    elapsed_time = round((tf-ti),2)
    print('Simulation Time: ' + str(elapsed_time) + ' seconds per timestep')   
        
    #plot_emission(df,Fuel,equivalence_ratio)
    
    return 
 
def combustor(dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, mdot_fuel, mdot_fuel_TO, mdot_air, mdot_air_TO, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des, l_SA_SM, l_SA_FM, F_SC, A_SZ, L_SZ, l_DA_start, l_DA_end, f_SM):
      
    f_air_PZ              = (mdot_fuel_TO*F_SC)/(phi_PZ_des*mdot_air_TO*FAR_st)                       # fraction of total air present in the combustor that enters the primary zone
    phi_sign              = (mdot_fuel*F_SC)/(f_air_PZ*mdot_air*FAR_st)                               # mean equivalence ratio
    sigma_phi             = S_PZ*phi_sign                                                             # standard deviation
    V_PZ_i                = V_PZ/N_PZ
    phi                   = np.linspace(0, 2*phi_sign, N_PZ)   
    
    PSRs                  = []
    mass_flow_controllers = []
    f_phi                 = np.zeros(len(phi))
    Fuel                  = ct.Solution('JetFuelSurrogate.yaml')
    rho                   = Fuel.density
    upstream              = ct.Reservoir(Fuel, name='upstream')
    mixer_1               = ct.IdealGasReactor(Fuel, name='mixer 1')
    mixer_2               = ct.IdealGasReactor(Fuel, name='mixer 2')
    comp_fuel             = list(dict_fuel.keys())
    Y_fuel                = Fuel[comp_fuel].Y     
    
    # ----------------------------------------------------------------
    # ----------------------------- PSRs -----------------------------
    # ---------------------------------------------------------------- 
    
    for i in range(21):

        Delta_phi = np.abs(phi[0] - phi[1])
        f_phi[i] = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi[i] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
        Fuel.TP               = T_stag_0, P_stag_0
        Fuel.set_equivalence_ratio(phi[i], fuel=dict_fuel, oxidizer=dict_oxy)
        Fuel.equilibrate('HP')
        
        # Create PSR for this iteration
        PSR = ct.IdealGasReactor(Fuel, name=f'PSR {i+1}')
        PSR.volume = V_PZ_i
        PSRs.append(PSR)
    
        # Calculate mass flow rates
        m_dot_fuel_i = (mdot_fuel * f_phi[i])/N_PZ
        m_dot_air_i = (mdot_air * f_air_PZ)/N_PZ
        mass_flow_rate = m_dot_fuel_i + m_dot_air_i
            
        # Create MassFlowController
        mfc = ct.MassFlowController(upstream, PSR)
        mfc.mass_flow_rate = mass_flow_rate
        mass_flow_controllers.append(mfc)
    
        # Set up the time for the reactor to equilibrate
        rho = Fuel.density
        t_res_psr_i = V_PZ_i / (rho * mass_flow_rate)
    
        # Simulate the PSR
        sim_psr = ct.ReactorNet([PSR])
        sim_psr.advance(t_res_psr_i)
    
    # Set up mixing for all the reactors' outlets to go into the final mixer
    for PSR in PSRs:
        outlet = ct.MassFlowController(PSR, mixer_1)
    
    # ----------------------------------------------------------------
    # --------------------------- Mixing -----------------------------
    # ----------------------------------------------------------------     
    
    PFR_1 = ct.IdealGasConstPressureReactor(Fuel, name='PFR 1')
    PFR_2 = ct.IdealGasConstPressureReactor(Fuel, name='PFR 2')
    outlet_1 = ct.Valve(mixer_1, PFR_1, K=10.0) 
    outlet_2 = ct.Valve(mixer_1, PFR_2, K=10.0) 
    
    # Simulate mixing
    sim_mixer = ct.ReactorNet([mixer_1])
    sim_mixer.advance_to_steady_state()    
    
    ## ----------------------------------------------------------------
    ## ---------------------------- PFR #1 ----------------------------
    ## ----------------------------------------------------------------
    
    f_FM       = 1 - f_SM
    f_air_SA   = mdot_fuel_TO/(phi_SZ_des*FAR_st*mdot_air_TO)
    f_air_DA   = 1 - f_air_PZ - f_air_SA
    
    beta_SA_SM = (f_air_SA*f_SM*mdot_air)/(l_SA_SM * L_SZ)
    beta_SA_FM = (f_air_SA*f_FM*mdot_air)/(l_SA_FM * L_SZ)
    beta_DA    = (f_air_DA*mdot_air)/((l_DA_end - l_DA_start)*L_SZ)
    
    sim = ct.ReactorNet([PFR_1])
    sim.rtol = 1e-6  # Set the relative tolerance
    sim.atol = 1e-14  # Set the absolute tolerance  
    n_segments = 200  # Number of segments for the reactor
    dz = L_SZ / n_segments  # Step size 
    
    # Add air into the reactor
    air = ct.Solution('air.yaml')
    air.TPX = T_stag_0, P_stag_0 ,dict_oxy    
    air_reservoir = ct.Reservoir(air)  # Create a reservoir of air for adding into the reactor
    mfc_air = ct.MassFlowController(air_reservoir, PFR_1)    
    
    for i in range(n_segments):
        z = i * dz
    
        if 0 <= z <= l_SA_SM * L_SZ:
            beta_air_in = beta_SA_SM                                       # magnitude of incoming airflow at any point in the secondary zone
        elif l_DA_start * L_SZ <= z <= l_DA_end * L_SZ:
            beta_air_in = beta_DA
        else:
            beta_air_in = 0 
            
        air_mass_flow = max(0, beta_air_in * Fuel.density * dz)
        mfc_air.mass_flow_rate = air_mass_flow
        outlet = ct.MassFlowController(PFR_1, mixer_2)
        outlet.mass_flow_rate = mdot_air + mdot_fuel + air_mass_flow
        
        sim.advance(z + dz / 10)
        
        # Print results for this segment
        print(f"Position: {z:.2f} m, Temperature: {Fuel.T:.2f} K, Pressure: {Fuel.P:.2f} Pa")    
        
    ## ----------------------------------------------------------------
    ## ---------------------------- PFR #2 ----------------------------
    ## ----------------------------------------------------------------
    
    sim = ct.ReactorNet([PFR_2])
    sim.rtol = 1e-6  # Set the relative tolerance
    sim.atol = 1e-14  # Set the absolute tolerance   
    n_segments = 200  # Number of segments for the reactor
    dz = L_SZ / n_segments  # Step size
    
    # Add air into the reactor
    air = ct.Solution('air.yaml')
    air.TPX = T_stag_0, P_stag_0 ,dict_oxy    
    air_reservoir = ct.Reservoir(air)  # Create a reservoir of air for adding into the reactor
    mfc_air = ct.MassFlowController(air_reservoir, PFR_2)    
    
    for i in range(n_segments):
        z = i * dz
    
        if 0 <= z <= l_SA_FM * L_SZ:
            beta_air_in = beta_SA_FM
        elif l_DA_start * L_SZ <= z <= l_DA_end * L_SZ:
            beta_air_in = beta_DA
        else:
            beta_air_in = 0 
            
        air_mass_flow = max(0, beta_air_in * Fuel.density * dz)
        mfc_air.mass_flow_rate = air_mass_flow
        outlet = ct.MassFlowController(PFR_2, mixer_2)
        outlet.mass_flow_rate = mdot_air + mdot_fuel + air_mass_flow
        
        sim.advance(z + dz / 10)
        
        # Print results for this segment
        print(f"Position: {z:.2f} m, Temperature: {Fuel.T:.2f} K, Pressure: {Fuel.P:.2f} Pa")         

    # ----------------------------------------------------------------
    # --------------------------- Mixing -----------------------------
    # ----------------------------------------------------------------     
    
    PFR_3 = ct.IdealGasConstPressureReactor(Fuel, name='PFR 3')
    outlet = ct.Valve(mixer_2, PFR_3, K=10.0)
    # Simulate mixing
    sim_mixer = ct.ReactorNet([mixer_2])
    sim_mixer.advance_to_steady_state()  
    
    # ----------------------------------------------------------------
    # --------------------- Additional computations ------------------
    # ----------------------------------------------------------------    
          
    # Determine Emission Indices 
    Emission_Index = Fuel.Y * (mdot_fuel + mdot_air)/mdot_fuel 
    
    # Extract properties of combustor flow 
    a_out      = Fuel.sound_speed  # Speed of sound at PFR outlet
    rho_out    = Fuel.density # density of the Fuel in the combustor
    gamma      = Fuel.cp_mass / Fuel.cv_mass
    h          = Fuel.h # enthalpy
    vel_out    = (mdot_fuel + mdot_air) / (rho_out * A_SZ)  # Outlet velocity (m/s)  
    M_out      = vel_out / a_out  # Outlet Mach number
    
    phi        = Fuel.equivalence_ratio(fuel = dict_fuel, oxidizer = dict_oxy) 
    
    # Stagnation temperature 
    T_stag_out = Fuel.T * (1 + 0.5 * (gamma - 1) * (M_out)**2)
    
    # stagnation pressure 
    P_stag_out = Fuel.P * (1 + 0.5 * (gamma - 1) * (M_out)**2)**(gamma / (gamma - 1))
    
    # Stagnation enthalpy 
    h_stag_out = T_stag_out  * Fuel.cp_mass
    
    # Fuel-to-air ratio (FAR)
    FAR      = mdot_fuel / (mdot_air)   
    
    return (Fuel, Emission_Index, T_stag_out, P_stag_out, h_stag_out, FAR) 

    
def plot_emission(df,Fuel,equivalence_ratio): 
    # Plot results
    f, ax1 = plt.subplots(3, 1, figsize=(16, 12))
    f.suptitle('Jet-A EI')
    subtitle = f'Equivalence ratio: {equivalence_ratio[0]}, Temperature: {Fuel.T:.1f} K, Pressure: {Fuel.P/ct.one_atm:.1f} atm,'
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
    #subtitle = f'Equivalence ratio: {equivalence_ratio[0]}, Temperature: {Fuel.T:.1f} K, Pressure: {Fuel.P/ct.one_atm:.1f} atm,'
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
    plt.plot(df['tau(s)'], df['T_stag_out'], '.-', color='C0')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('T_stag_out [K]')
    plt.title('T_stag_out vs. PFR residence time')
    plt.grid(True)    
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['P_stag_out'], '.-', color='C0')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('P_stag_out [atm]')
    plt.title('P_stag_out vs. PFR residence time')
    plt.grid(True) 
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['h_stag_out'], '.-', color='C0')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('h_stag_out [J]')
    plt.title('h_stag_out vs. PFR residence time')
    plt.grid(True)  
    
    plt.figure(figsize=(16, 12))
    plt.plot(df['tau(s)'], df['FAR'], '.-', color='C0')
    plt.xlabel('PFR residence time [s]')
    plt.ylabel('FAR [-]')
    plt.title('FAR vs. PFR residence time')
    plt.grid(True)      
    
    plt.show() 
 
    return  

if __name__ == '__main__': 
    main()
    plt.show()