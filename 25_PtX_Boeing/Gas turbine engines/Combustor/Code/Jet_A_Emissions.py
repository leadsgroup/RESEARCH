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
    phi_PZ_des              = 1.7                                       # [-]
    S_PZ                    = 0.39                                      # [-]           
    phi_SZ_des              = 0.7                                       # [-]           
    l_SA_SM                 = 0.55                                      # [-]          
    l_SA_FM                 = 0.055                                     # [-]
    LHV_input_fuel          = 1 #FIX
    LHV_model_fuel          = 1 #FIX
    F_SC                    = LHV_input_fuel/LHV_model_fuel             # [-]
    
    #Secondary Zone inputs
    A_SZ                    = 0.15                                      # [m**2]
    L_SZ                    = 0.075                                     # [m]           
    l_DA_start              = 0.95                                      # [-]          
    l_DA_end                = 0.1                                       # [-] 
    f_SM                    = 0.5                                       # [-]
    
    dict_fuel = {'N-C12H26':0.6, 'A1CH3':0.2, 'A1':0.2} # Less accurate model (no NOx), faster
    #dict_fuel = {'NC10H22':0.16449, 'NC12H26':0.34308, 'NC16H34':0.10335, 'IC8H18':0.08630, 'NC7H14':0.07945, 'C6H5C2H5': 0.07348, 'C6H5C4H9': 0.05812, 'C10H7CH3': 0.10972}      # More accurate model (NOx), slower   
    dict_oxy = {'O2':0.2095, 'N2':0.7809, 'AR':0.0093, 'CO2':0.0003}    

    #--------------------------------------------------------------------------------
    
    list_sp = ['CO2', 'CO', 'H2O']                          # Less accurate model (no NOx), faster
    #list_sp   = ['CO2', 'CO', 'H2O', 'NO', 'NO2', 'CSOLID'] # More accurate model (NOx), slower
    col_names = ['tau(s)', 'Tout(K)', 'T_stag_out','P_stag_out', 'h_stag_out', 'FAR'] + ['X_' +str(sp) for sp in list_sp] + ['Y_' +str(sp) for sp in list_sp] + ['EI_' +str(sp) for sp in list_sp]
    df        = pd.DataFrame(columns=col_names)
    
    Fuel, EI, T_stag_out, P_stag_out, h_stag_out, FAR = combustor(dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, mdot_fuel, mdot_fuel_TO, mdot_air, mdot_air_TO, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des, l_SA_SM, l_SA_FM, F_SC, A_SZ, L_SZ, l_DA_start, l_DA_end, f_SM)
    #sp_idx = [Fuel.species_index(sp) for sp in list_sp]
    #data_n = [residence_time_pfr[n], Fuel.T, T_stag_out, P_stag_out, h_stag_out, FAR] + list(Fuel.X[sp_idx]) + list(Fuel.Y[sp_idx]) + list(EI[sp_idx])
    #df.loc[n] = data_n
    #print(EI[sp_idx])
    
    tf           = time.time()
    elapsed_time = round((tf-ti)/len(residence_time_pfr),2)
    print('Simulation Time: ' + str(elapsed_time) + ' seconds per timestep')   
        
    #plot_emission(df,Fuel,equivalence_ratio)
    
    return 
 
def combustor(dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, mdot_fuel, mdot_fuel_TO, mdot_air, mdot_air_TO, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des, l_SA_SM, l_SA_FM, F_SC, A_SZ, L_SZ, l_DA_start, l_DA_end, f_SM):
        
    Fuel = ct.Solution('JetFuelSurrogate.yaml')
    Fuel.TP = T_stag_0, P_stag_0
    
    sigma_phi = S_PZ*phi_PZ_des
    f_air_PZ  = (mdot_fuel_TO*F_SC)/(phi_PZ_des*mdot_air_TO*FAR_st)
    phi_sign       = (mdot_fuel*F_SC)/(f_air_PZ*mdot_air*FAR_st)
    V_PZ_i                   = V_PZ/N_PZ
    
    upstream   = ct.Reservoir(Fuel)
    mixer      = ct.IdealGasReactor(Fuel)    
    
    for i in N_PZ:
        f_phi[i]   = (1/(np.sqrt(2*np.pi)*sigma_phi)) * np.exp((-(phi_sign[i] - phi_PZ_des)**2)/(2*sigma_phi**2))*Delta_phi

        Fuel.set_equivalence_ratio(phi[i], fuel = dict_fuel, oxidizer = dict_oxy)
        Fuel.equilibrate('HP')
        PSR[i]     = ct.IdealGasReactor(Fuel)
        PSR[i].volume               = V_PZ_i
        rho[i]                      = Fuel.density
        t_res_psr_i              = V_PZ_i/(rho_i*(m_dot_fuel_i + m_dot_air_i))
        func_mdot  = lambda t: PSR[i].mass/t_res_psr_i
        inlet                = ct.MassFlowController(upstream, PSR[[i]])
        inlet.mass_flow_rate = func_mdot*f_phi[i]
        outlet  = ct.Valve(psr, mixer, K=100) 
        sim_psr = ct.ReactorNet([PSR[i]])
        sim_psr.advance(t_res_psr_i)         
    
    
   
    
    
    comp_fuel                = list(dict_fuel.keys())
    Y_fuel                   = Fuel[comp_fuel].Y    
    
    Air_composition          = 'O2:0.21, N2:0.78, AR:0.01'
    mdot_air                 = 14*0.7
    mdot_air_d               = 14*0.3
    T_injection_air          = temp
    P_injection_air          = press  
    Air                      = ct.Solution('air.yaml')
    Air.TPX                  = T_injection_air, P_injection_air, Air_composition   
    
    dilution_air_1           = Air
    dilution_air_1.TPX       = T_injection_air, P_injection_air, Air_composition 
    
    dilution_air_2           = Air
    dilution_air_2.TPX       = T_injection_air, P_injection_air, Air_composition 
    
    dilution_air_3           = Air
    dilution_air_3.TPX       = T_injection_air, P_injection_air, Air_composition     
    
    Fuel_composition         = 'N-C12H26:0.6, A1CH3:0.2, A1:0.2'
    mdot_fuel                = 0.65
    T_injection_fuel         = 300
    P_injection_fuel         = ct.one_atm 
    Fuel                     = ct.Solution('JetFuelSurrogate.yaml')  
    Fuel.TPX                 = T_injection_fuel, P_injection_fuel, Fuel_composition 
    
    mdot_tot                 = mdot_air + mdot_air_d + mdot_fuel    
    
    # ----------------------------------------------------------------
    # ----------------------- CRN Components -------------------------
    # ----------------------------------------------------------------
    
    res_air                  = ct.Reservoir(Air, name='air')
    res_fuel                 = ct.Reservoir(Fuel, name='fuel')
    mixer_1                  = ct.IdealGasReactor(Fuel, name='mixer 1')
    PSR                      = ct.IdealGasReactor(Fuel, name='PSR') 
    res_d_1                  = ct.Reservoir(dilution_air_1, name='dilution air 1')
    mixer_2                  = ct.IdealGasReactor(Fuel, name='mixer 2') 
    PFR_1                    = ct.IdealGasConstPressureReactor(Fuel, name='PFR 1')
    res_d_2                  = ct.Reservoir(dilution_air_2, name='dilution air 2')
    mixer_3                  = ct.IdealGasReactor(Fuel, name='mixer 3') 
    PFR_2                    = ct.IdealGasConstPressureReactor(Fuel, name='PFR 2')   
    res_d_3                  = ct.Reservoir(dilution_air_3, name='dilution air 3')
    mixer_4                  = ct.IdealGasReactor(Fuel, name='mixer 4')    
    PFR_3                    = ct.IdealGasConstPressureReactor(Fuel, name='PFR 2') 
    
    # ----------------------------------------------------------------
    # ---------------------- Air + Fuel Mixing -----------------------
    # ---------------------------------------------------------------- 
    
    mfc1_1                   = ct.MassFlowController(res_air, mixer_1, mdot = mdot_air)
    mfc2_1                   = ct.MassFlowController(res_fuel, mixer_1, mdot = mdot_fuel)    
    outlet_1                 = ct.Valve(mixer_1, PSR, K=10.0)
    sim_mixer_1              = ct.ReactorNet([mixer_1])    
    sim_mixer_1.advance_to_steady_state()
    #sim_mixer_1.advance(residence_time_pfr)     
    
    # ----------------------------------------------------------------
    # ----------------------------- PSRs -----------------------------
    # ----------------------------------------------------------------    
    
    PSR.volume               = 1
    #Fuel.TP                  = temp, press
    Fuel.equilibrate('HP')
    #Fuel.set_equivalence_ratio(equivalence_ratio, fuel = dict_fuel, oxidizer = dict_oxy )  
    outlet_2                 = ct.Valve(PSR, mixer_2, K=100) 
    sim_psr                  = ct.ReactorNet([PSR])
    #sim_psr.advance_to_steady_state()
    sim_psr.advance(residence_time_psr)
    
    # ----------------------------------------------------------------
    # ----------------- Addition of dilution air #1 ------------------
    # ----------------------------------------------------------------
    
    mfc1_2                   = ct.MassFlowController(res_d_1, mixer_2, mdot=mdot_air_d/3)   
    outlet_3                 = ct.Valve(mixer_2, PFR_1, K=10.0)
    sim_mixer_2              = ct.ReactorNet([mixer_2])    
    sim_mixer_2.advance_to_steady_state()
    #sim_mixer_2.advance(residence_time_pfr)
    
    # ----------------------------------------------------------------
    # ---------------------------- PFR #1 ----------------------------
    # ----------------------------------------------------------------
    
    outlet_4                 = ct.Valve(PFR_1, mixer_3, K=100)
    sim_pfr_1                = ct.ReactorNet([PFR_1])
    #sim_pfr_1.advance_to_steady_state() 
    sim_pfr_1.advance(residence_time_pfr)
    
    # ----------------------------------------------------------------
    # ----------------- Addition of dilution air #2 ------------------
    # ---------------------------------------------------------------- 
    
    mfc1_3                   = ct.MassFlowController(res_d_2, mixer_3, mdot=mdot_air_d/3)   
    outlet_5                 = ct.Valve(mixer_3, PFR_2, K=10.0)
    sim_mixer_3              = ct.ReactorNet([mixer_3])    
    sim_mixer_3.advance_to_steady_state()
    #sim_mixer_3.advance(residence_time_pfr)    
    
    # ----------------------------------------------------------------
    # ---------------------------- PFR #2 ----------------------------
    # ----------------------------------------------------------------
    
    outlet_6                 = ct.Valve(PFR_2, mixer_4, K=100)
    sim_pfr_2                = ct.ReactorNet([PFR_2])
    #sim_pfr_2.advance_to_steady_state() 
    sim_pfr_2.advance(residence_time_pfr)  
    
    # ----------------------------------------------------------------
    # ------------------- Addition of cooling air --------------------
    # ---------------------------------------------------------------- 
    
    mfc1_4                   = ct.MassFlowController(res_d_3, mixer_4, mdot=mdot_air_d/3)  
    outlet_7                 = ct.Valve(mixer_4, PFR_3, K=10.0)
    sim_mixer_4              = ct.ReactorNet([mixer_4])    
    sim_mixer_4.advance_to_steady_state()
    #sim_mixer_4.advance(residence_time_pfr)    
    
    # ----------------------------------------------------------------
    # ---------------------------- PFR #3 ----------------------------
    # ----------------------------------------------------------------
    
    sim_pfr_3                = ct.ReactorNet([PFR_3])
    #sim_pfr_3.advance_to_steady_state() 
    sim_pfr_3.advance(residence_time_pfr)     
    
    # ----------------------------------------------------------------
    # --------------------- Additional computations ------------------
    # ----------------------------------------------------------------    
          
    # Determine Emission Indices 
    Emission_Index = Fuel.Y * mdot_tot/mdot_fuel 
    
    # Extract properties of combustor flow 
    a_out      = Fuel.sound_speed  # Speed of sound at PFR outlet
    rho_out    = Fuel.density # density of the Fuel in the combustor
    gamma      = Fuel.cp_mass / Fuel.cv_mass
    h          = Fuel.h # enthalpy
    vel_out    = mdot / (rho_out * area_out)  # Outlet velocity (m/s)  
    M_out      = vel_out / a_out  # Outlet Mach number
    
    phi        = Fuel.equivalence_ratio(fuel = dict_fuel, oxidizer = dict_oxy) 
    
    # Stagnation temperature 
    T_stag_out = Fuel.T * (1 + 0.5 * (gamma - 1) * (M_out)**2)
    
    # stagnation pressure 
    P_stag_out = Fuel.P * (1 + 0.5 * (gamma - 1) * (M_out)**2)**(gamma / (gamma - 1))
    
    # Stagnation enthalpy 
    h_stag_out = T_stag_out  * Fuel.cp_mass
    
    # Fuel-to-air ratio (FAR)
    FAR      = mdot_fuel / (mdot_air + mdot_air_d)   
    
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