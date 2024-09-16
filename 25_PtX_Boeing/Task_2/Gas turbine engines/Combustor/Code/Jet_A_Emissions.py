import cantera           as ct
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
import time 

def main():

    ti        = time.time()
    
    # ----------------------------------------------------------------
    # ----------------------- Combustor Inputs -----------------------
    # ---------------------------------------------------------------- 

    # General Inputs 
    T_stag_0                = 700                                       # [K]    Stagnation Temperature entering all combustors
    P_stag_0                = 1500000                                   # [Pa]   Stagnation Pressure entering all combustors
    FAR                     = 0.02                                      # [-]    Fuel-to-Air ratio
    FAR_TO                  = 0.0275                                    # [-]    Fuel-to-Air ratio during TO 
    FAR_st                  = 0.068                                     # [-]    Stoichiometric Fuel-to-Air ratio
    m_dot_air_tot           = 40                                       # [kg/s] Air mass flow going through all combustors
    m_dot_air_TO_tot        = 44                                       # [kg/s] Air mass flow going through all combustors during TO 
    m_dot_fuel_tot          = m_dot_air_tot*FAR                        # [kg/s] Fuel mass flow going through all combustors
    m_dot_fuel_TO_tot       = m_dot_air_TO_tot*FAR_TO                  # [kg/s] Fuel mass flow going through all combustors during TO 
    N_comb                  = 9                                         # [-]    Number of can-annular combustors
    m_dot_air_id            = m_dot_air_tot/N_comb                     # [kg/s] Ideal Air mass flow inside each combustor, scaled inside each PSR to vary the Equivalence Ratio
    m_dot_air_TO            = m_dot_air_TO_tot/N_comb                  # [kg/s] Air mass flow inside each combustor during TO
    m_dot_fuel              = m_dot_fuel_tot/N_comb                    # [kg/s] Fuel mass flow inside each combustor
    m_dot_fuel_TO           = m_dot_fuel_TO_tot/N_comb                 # [kg/s] Fuel mass flow inside each combustor during TO    
    
    # Primary Zone Inputs
    N_PZ                    = 8                                         # [-]    Number of PSR (EVEN)
    V_PZ                    = 0.0023                                    # [m**3] Volume of the Primary Zone in a SINGLE combustor, must be split into the different PSRs       
    phi_PZ_des              = 1.77                                      # [-]    Design Equivalence Ratio
    S_PZ                    = 0.39                                      # [-]    Mixing parameter, used to define the Equivalence Ratio standard deviation         
    LHV_input_fuel          = 1                                         # [-]    FIX VALUE
    LHV_model_fuel          = 1                                         # [-]    FIX VALUE
    F_SC                    = LHV_input_fuel/LHV_model_fuel             # [-]    Fuel scaler, used to define the fraction of total air present in the combustor that enters the Primary Zone
    
    #Secondary Zone Inputs
    A_SZ                    = 0.15                                      # [m**2] Secondary Zone cross-sectional area
    L_SZ                    = 0.075                                     # [m]    Secondary Zone length
    t_mix                   = 0.001                                     # [s]    Mixing residence time
    n_segments              = 200                                       # [-]    Number of segments for each PFR
    f_SM                    = 0.5                                       # [-]     
    phi_SZ_des              = 0.7                                       # [-]
    
    Fuel                    = ct.Solution('JetFuelSurrogate.yaml')
    Air                     = ct.Solution('Air.yaml')
    dict_fuel               = {'N-C12H26':0.6, 'A1CH3':0.2, 'A1':0.2}
    dict_oxy                = {'O2':0.2095,    'N2':0.7809, 'AR':0.0096}    

    #--------------------------------------------------------------------------------
    
    list_sp   = ['CO2', 'CO', 'H2O']
    col_names = ['Tout(K)', 'T_stag_out','P_stag_out', 'h_stag_out', 'FAR'] + ['X_' +str(sp) for sp in list_sp] + ['Y_' +str(sp) for sp in list_sp] + ['EI_' +str(sp) for sp in list_sp]
    df        = pd.DataFrame(columns=col_names)
    
    for n in range(1):
        gas, EI, T_stag_out, P_stag_out, h_stag_out, FAR = combustor(Fuel, Air, dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, m_dot_fuel, m_dot_fuel_TO, m_dot_air_id, m_dot_air_TO, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des, F_SC, A_SZ, L_SZ, f_SM, t_mix)
        sp_idx = [gas.species_index(sp) for sp in list_sp]
        data_n = [gas.T, T_stag_out, P_stag_out, h_stag_out, FAR] + list(gas.X[sp_idx]) + list(gas.Y[sp_idx]) + list(EI[sp_idx])
        df.loc[n] = data_n
        
    print(df['EI_CO2'])
    print(df['EI_CO'])
    print(df['EI_H2O'])        
    
    tf           = time.time()
    elapsed_time = round((tf-ti),2)
    print('Simulation Time: ' + str(elapsed_time) + ' seconds per timestep')   
    
    return 
 
def combustor(Fuel, Air, dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, m_dot_fuel, m_dot_fuel_TO, m_dot_air_id, m_dot_air_TO, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des, F_SC, A_SZ, L_SZ, f_SM, t_mix):
    
    # ----------------------------------------------------------------
    # ---------------------- Initial Parameters ----------------------
    # ----------------------------------------------------------------    
      
    f_air_PZ               = (m_dot_fuel_TO*F_SC)/(phi_PZ_des*m_dot_air_TO*FAR_st) # Fraction of total air present in the combustor that enters the Primary Zone
    f_air_SZ               = 1 - f_air_PZ                                          # Fraction of total air present in the combustor that enters the Secondary Zone  
    m_dot_air_PZ           = f_air_PZ*m_dot_air_id                                 # Air mass flow going through the Primary Zone
    m_dot_air_SZ           = (f_air_SZ*m_dot_air_id)/3                             # Air mass flow going through each dilution air inlet (3 inlets)
    phi_sign               = ((m_dot_fuel*F_SC)/m_dot_air_PZ)/(FAR_st)             # Mean Equivalence Ratio
    sigma_phi              = S_PZ*phi_sign                                         # Standard deviation of the Equivalence Ratio    
    m_dot_air_PSR          = m_dot_air_PZ/N_PZ                                     # Air mass flow going through each PSR
    m_dot_fuel_PSR         = m_dot_fuel/N_PZ                                       # Fuel mass flow going through each PSR
    V_PZ_PSR               = V_PZ/N_PZ                                             # Volume of each PSR
    phi_PSR                = np.linspace(0, 2*phi_sign, N_PZ)                      # Distribution of Equivalence Ratio through the PSRs
    Delta_phi              = np.abs(phi_PSR[0] - phi_PSR[1])                       # Difference between two subsequent Equivalence Ratios
    comp_fuel              = list(dict_fuel.keys())                                # Fuel components
    
    # ----------------------------------------------------------------
    # ---------------------------- PSR #1 ----------------------------
    # ---------------------------------------------------------------- 
    
    f_fuel_PZ_1            = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[0] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
    Fuel_1                 = Fuel
    Fuel_1.TP              = T_stag_0, P_stag_0
    Fuel_1.set_equivalence_ratio(phi_PSR[0], fuel=dict_fuel, oxidizer=dict_oxy)
    Fuel_1.equilibrate('HP')   
    rho_1                  = Fuel_1.density
    m_dot_fuel_1           = m_dot_fuel_PSR * f_fuel_PZ_1
    mass_flow_rate_1       = m_dot_fuel_1 + m_dot_air_PSR  
    upstream_1             = ct.Reservoir(Fuel_1) 
    mixer_12               = ct.IdealGasReactor(Fuel_1)
    PSR_1                  = ct.IdealGasReactor(Fuel_1)
    PSR_1.volume           = V_PZ_PSR
    inlet_1                = ct.MassFlowController(upstream_1, PSR_1)
    inlet_1.mass_flow_rate = mass_flow_rate_1    
    outlet_1               = ct.MassFlowController(PSR_1, mixer_12, mdot=mass_flow_rate_1)                                                   
    t_res_PSR_1            = (rho_1 * V_PZ_PSR) / (mass_flow_rate_1)
    sim_PSR_1              = ct.ReactorNet([PSR_1])
    sim_PSR_1.advance(t_res_PSR_1)
    #sim_psr.advance_to_steady_state()
    Y_fuel_1               = Fuel_1[comp_fuel].Y
    
    # ----------------------------------------------------------------
    # ---------------------------- PSR #2 ----------------------------
    # ---------------------------------------------------------------- 
    
    f_fuel_PZ_2            = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[1] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
    Fuel_2                 = Fuel
    Fuel_2.TP              = T_stag_0, P_stag_0
    Fuel_2.set_equivalence_ratio(phi_PSR[1], fuel=dict_fuel, oxidizer=dict_oxy)
    Fuel_2.equilibrate('HP')   
    rho_2                  = Fuel_2.density
    m_dot_fuel_2           = m_dot_fuel_PSR * f_fuel_PZ_2
    mass_flow_rate_2       = m_dot_fuel_2 + m_dot_air_PSR  
    upstream_2             = ct.Reservoir(Fuel_2) 
    PSR_2                  = ct.IdealGasReactor(Fuel_2)
    PSR_2.volume           = V_PZ_PSR
    inlet_2                = ct.MassFlowController(upstream_2, PSR_2)
    inlet_2.mass_flow_rate = mass_flow_rate_2    
    outlet_2               = ct.MassFlowController(PSR_2, mixer_12, mdot=mass_flow_rate_2)                                                           
    t_res_PSR_2            = (rho_2 * V_PZ_PSR) / (mass_flow_rate_2)
    sim_PSR_2              = ct.ReactorNet([PSR_2])
    sim_PSR_2.advance(t_res_PSR_2)
    #sim_psr.advance_to_steady_state()
    Y_fuel_2               = Fuel_2[comp_fuel].Y 
    
    # ----------------------------------------------------------------
    # ---------------------------- PSR #3 ----------------------------
    # ---------------------------------------------------------------- 
    
    f_fuel_PZ_3            = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[2] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
    Fuel_3                 = Fuel
    Fuel_3.TP              = T_stag_0, P_stag_0
    Fuel_3.set_equivalence_ratio(phi_PSR[2], fuel=dict_fuel, oxidizer=dict_oxy)
    Fuel_3.equilibrate('HP')   
    rho_3                  = Fuel_3.density
    m_dot_fuel_3           = m_dot_fuel_PSR * f_fuel_PZ_3
    mass_flow_rate_3       = m_dot_fuel_3 + m_dot_air_PSR  
    upstream_3             = ct.Reservoir(Fuel_3) 
    mixer_34               = ct.IdealGasReactor(Fuel_3)
    PSR_3                  = ct.IdealGasReactor(Fuel_3)
    PSR_3.volume           = V_PZ_PSR
    inlet_3                = ct.MassFlowController(upstream_3, PSR_3)
    inlet_3.mass_flow_rate = mass_flow_rate_3    
    outlet_3               = ct.MassFlowController(PSR_3, mixer_34, mdot=mass_flow_rate_3)                                                          
    t_res_PSR_3            = (rho_3 * V_PZ_PSR) / (mass_flow_rate_3)
    sim_PSR_3              = ct.ReactorNet([PSR_3])
    sim_PSR_3.advance(t_res_PSR_3)
    #sim_psr.advance_to_steady_state()
    Y_fuel_3               = Fuel_3[comp_fuel].Y    
    
    # ----------------------------------------------------------------
    # ---------------------------- PSR #4 ----------------------------
    # ---------------------------------------------------------------- 
    
    f_fuel_PZ_4            = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[3] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
    Fuel_4                 = Fuel
    Fuel_4.TP              = T_stag_0, P_stag_0
    Fuel_4.set_equivalence_ratio(phi_PSR[3], fuel=dict_fuel, oxidizer=dict_oxy)
    Fuel_4.equilibrate('HP')   
    rho_4                  = Fuel_4.density
    m_dot_fuel_4           = m_dot_fuel_PSR * f_fuel_PZ_4
    mass_flow_rate_4       = m_dot_fuel_4 + m_dot_air_PSR  
    upstream_4             = ct.Reservoir(Fuel_4) 
    PSR_4                  = ct.IdealGasReactor(Fuel_4)
    PSR_4.volume           = V_PZ_PSR
    inlet_4                = ct.MassFlowController(upstream_4, PSR_4)
    inlet_4.mass_flow_rate = mass_flow_rate_4    
    outlet_4               = ct.MassFlowController(PSR_4, mixer_34, mdot=mass_flow_rate_4)                                                          
    t_res_PSR_4            = (rho_4 * V_PZ_PSR) / (mass_flow_rate_4)
    sim_PSR_4              = ct.ReactorNet([PSR_4])
    sim_PSR_4.advance(t_res_PSR_4)
    #sim_psr.advance_to_steady_state()
    Y_fuel_4               = Fuel_4[comp_fuel].Y
    
    # ----------------------------------------------------------------
    # ---------------------------- PSR #5 ----------------------------
    # ---------------------------------------------------------------- 
    
    f_fuel_PZ_5            = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[4] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
    Fuel_5                 = Fuel
    Fuel_5.TP              = T_stag_0, P_stag_0
    Fuel_5.set_equivalence_ratio(phi_PSR[4], fuel=dict_fuel, oxidizer=dict_oxy)
    Fuel_5.equilibrate('HP')   
    rho_5                  = Fuel_5.density
    m_dot_fuel_5           = m_dot_fuel_PSR * f_fuel_PZ_5
    mass_flow_rate_5       = m_dot_fuel_5 + m_dot_air_PSR  
    upstream_5             = ct.Reservoir(Fuel_5) 
    mixer_56               = ct.IdealGasReactor(Fuel_5)
    PSR_5                  = ct.IdealGasReactor(Fuel_5)
    PSR_5.volume           = V_PZ_PSR
    inlet_5                = ct.MassFlowController(upstream_5, PSR_5)
    inlet_5.mass_flow_rate = mass_flow_rate_5    
    outlet_5               = ct.MassFlowController(PSR_5, mixer_56, mdot=mass_flow_rate_5)                                                          
    t_res_PSR_5            = (rho_5 * V_PZ_PSR) / (mass_flow_rate_5)
    sim_PSR_5              = ct.ReactorNet([PSR_5])
    sim_PSR_5.advance(t_res_PSR_5)
    #sim_psr.advance_to_steady_state()
    Y_fuel_5               = Fuel_5[comp_fuel].Y    
    
    # ----------------------------------------------------------------
    # ---------------------------- PSR #6 ----------------------------
    # ---------------------------------------------------------------- 
    
    f_fuel_PZ_6            = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[5] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
    Fuel_6                 = Fuel
    Fuel_6.TP              = T_stag_0, P_stag_0
    Fuel_6.set_equivalence_ratio(phi_PSR[5], fuel=dict_fuel, oxidizer=dict_oxy)
    Fuel_6.equilibrate('HP')   
    rho_6                  = Fuel_6.density
    m_dot_fuel_6           = m_dot_fuel_PSR * f_fuel_PZ_6
    mass_flow_rate_6       = m_dot_fuel_6 + m_dot_air_PSR  
    upstream_6             = ct.Reservoir(Fuel_6) 
    PSR_6                  = ct.IdealGasReactor(Fuel_6)
    PSR_6.volume           = V_PZ_PSR
    inlet_6                = ct.MassFlowController(upstream_6, PSR_6)
    inlet_6.mass_flow_rate = mass_flow_rate_6    
    outlet_6               = ct.MassFlowController(PSR_6, mixer_56, mdot=mass_flow_rate_6)                                                          
    t_res_PSR_6            = (rho_6 * V_PZ_PSR) / (mass_flow_rate_6)
    sim_PSR_6              = ct.ReactorNet([PSR_6])
    sim_PSR_6.advance(t_res_PSR_6)
    #sim_psr.advance_to_steady_state()
    Y_fuel_6               = Fuel_6[comp_fuel].Y 
    
    # ----------------------------------------------------------------
    # ---------------------------- PSR #7 ----------------------------
    # ---------------------------------------------------------------- 
    
    f_fuel_PZ_7            = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[6] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
    Fuel_7                 = Fuel
    Fuel_7.TP              = T_stag_0, P_stag_0
    Fuel_7.set_equivalence_ratio(phi_PSR[6], fuel=dict_fuel, oxidizer=dict_oxy)
    Fuel_7.equilibrate('HP')   
    rho_7                  = Fuel_7.density
    m_dot_fuel_7           = m_dot_fuel_PSR * f_fuel_PZ_7
    mass_flow_rate_7       = m_dot_fuel_7 + m_dot_air_PSR  
    upstream_7             = ct.Reservoir(Fuel_7) 
    mixer_78               = ct.IdealGasReactor(Fuel_7)
    PSR_7                  = ct.IdealGasReactor(Fuel_7)
    PSR_7.volume           = V_PZ_PSR
    inlet_7                = ct.MassFlowController(upstream_7, PSR_7)
    inlet_7.mass_flow_rate = mass_flow_rate_7    
    outlet_7               = ct.MassFlowController(PSR_7, mixer_78, mdot=mass_flow_rate_7)                                                          
    t_res_PSR_7            = (rho_7 * V_PZ_PSR) / (mass_flow_rate_7)
    sim_PSR_7              = ct.ReactorNet([PSR_7])
    sim_PSR_7.advance(t_res_PSR_7)
    #sim_psr.advance_to_steady_state()
    Y_fuel_7               = Fuel_7[comp_fuel].Y
    
    # ----------------------------------------------------------------
    # ---------------------------- PSR #8 ----------------------------
    # ---------------------------------------------------------------- 
    
    f_fuel_PZ_8            = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[7] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi  # Fraction of mass flow entering reactor i at equivalence ratio phi_i
    Fuel_8                 = Fuel
    Fuel_8.TP              = T_stag_0, P_stag_0
    Fuel_8.set_equivalence_ratio(phi_PSR[7], fuel=dict_fuel, oxidizer=dict_oxy)
    Fuel_8.equilibrate('HP')   
    rho_8                  = Fuel_8.density
    m_dot_fuel_8           = m_dot_fuel_PSR * f_fuel_PZ_8
    mass_flow_rate_8       = m_dot_fuel_8 + m_dot_air_PSR  
    upstream_8             = ct.Reservoir(Fuel_8) 
    PSR_8                  = ct.IdealGasReactor(Fuel_8)
    PSR_8.volume           = V_PZ_PSR
    inlet_8                = ct.MassFlowController(upstream_8, PSR_8)
    inlet_8.mass_flow_rate = mass_flow_rate_8    
    outlet_8               = ct.MassFlowController(PSR_8, mixer_78, mdot=mass_flow_rate_8)                                                          
    t_res_PSR_8            = (rho_8 * V_PZ_PSR) / (mass_flow_rate_8)
    sim_PSR_8              = ct.ReactorNet([PSR_8])
    sim_PSR_8.advance(t_res_PSR_8)
    #sim_psr.advance_to_steady_state()
    Y_fuel_8               = Fuel_8[comp_fuel].Y
    
    # ----------------------------------------------------------------
    # -------------------------- Mixing 1-2 --------------------------
    # ----------------------------------------------------------------     
    
    mixer_1234             = ct.IdealGasReactor(Fuel_1)
    outlet_12              = ct.MassFlowController(mixer_12, mixer_1234, mdot = (mass_flow_rate_1 + mass_flow_rate_2))  
    sim_mixer_12           = ct.ReactorNet([mixer_12])
    sim_mixer_12.advance_to_steady_state()   
    
    # ----------------------------------------------------------------
    # -------------------------- Mixing 3-4 --------------------------
    # ----------------------------------------------------------------     
    
    outlet_34              = ct.MassFlowController(mixer_34, mixer_1234, mdot = (mass_flow_rate_3 + mass_flow_rate_4))  
    sim_mixer_34           = ct.ReactorNet([mixer_34])
    sim_mixer_34.advance_to_steady_state()   
    
    # ----------------------------------------------------------------
    # -------------------------- Mixing 5-6 --------------------------
    # ----------------------------------------------------------------     
    
    mixer_5678             = ct.IdealGasReactor(Fuel_5)
    outlet_56              = ct.MassFlowController(mixer_56, mixer_5678, mdot = (mass_flow_rate_5 + mass_flow_rate_6))  
    sim_mixer_56           = ct.ReactorNet([mixer_56])
    sim_mixer_56.advance_to_steady_state()   
    
    # ----------------------------------------------------------------
    # -------------------------- Mixing 7-8 --------------------------
    # ----------------------------------------------------------------     
    
    outlet_78              = ct.MassFlowController(mixer_78, mixer_5678, mdot = (mass_flow_rate_7 + mass_flow_rate_8))  
    sim_mixer_78           = ct.ReactorNet([mixer_78])
    sim_mixer_78.advance_to_steady_state()
    
    # ----------------------------------------------------------------
    # ------------------------ Mixing 1-2-3-4 ------------------------
    # ----------------------------------------------------------------     
    
    mixer_12345678         = ct.IdealGasReactor(Fuel_1)
    outlet_1234            = ct.MassFlowController(mixer_1234, mixer_12345678, mdot = (mass_flow_rate_1 + mass_flow_rate_2 + mass_flow_rate_3 + mass_flow_rate_4))  
    sim_mixer_1234         = ct.ReactorNet([mixer_1234])
    sim_mixer_1234.advance_to_steady_state()
    
    # ----------------------------------------------------------------
    # ------------------------ Mixing 5-6-7-8 ------------------------
    # ----------------------------------------------------------------     

    outlet_5678            = ct.MassFlowController(mixer_5678, mixer_12345678, mdot = (mass_flow_rate_5 + mass_flow_rate_6 + mass_flow_rate_7 + mass_flow_rate_8))  
    sim_mixer_5678         = ct.ReactorNet([mixer_5678])
    sim_mixer_5678.advance_to_steady_state()
    
    # ----------------------------------------------------------------
    # -------------------- Mixing 1-2-3-4-5-6-7-8 --------------------
    # ----------------------------------------------------------------     

    mixer_air_1            = ct.IdealGasReactor(Fuel_1)
    outlet_12345678        = ct.MassFlowController(mixer_12345678, mixer_air_1, mdot = (mass_flow_rate_1 + mass_flow_rate_2 + mass_flow_rate_3 + mass_flow_rate_4 + mass_flow_rate_5 + mass_flow_rate_6 + mass_flow_rate_7 + mass_flow_rate_8))  
    sim_mixer_12345678     = ct.ReactorNet([mixer_12345678])
    sim_mixer_12345678.advance_to_steady_state()   
    
    # ----------------------------------------------------------------
    # ------------------------- Mixing 1-Air -------------------------
    # ----------------------------------------------------------------     

    Air_1                  = Air
    Air_1.TPX              = T_stag_0, P_stag_0, dict_oxy
    rho_air_1              = Air_1.density    
    res_air_1              = ct.Reservoir(Air_1)
    PFR_1                  = ct.IdealGasConstPressureReactor(Fuel_1)  
    PFR_1.volume           = A_SZ*(L_SZ/3)
    inlet_air_1            = ct.MassFlowController(res_air_1, mixer_air_1, mdot=m_dot_air_SZ)
    m_dot_air_1            = mass_flow_rate_1 + mass_flow_rate_2 + mass_flow_rate_3 + mass_flow_rate_4 + mass_flow_rate_5 + mass_flow_rate_6 + mass_flow_rate_7 + mass_flow_rate_8 + m_dot_air_SZ
    outlet_air_1           = ct.MassFlowController(mixer_air_1, PFR_1, mdot=m_dot_air_SZ)
    sim_mixer_air_1        = ct.ReactorNet([mixer_air_1])
    sim_mixer_air_1.advance(t_mix)
    
    # ----------------------------------------------------------------
    # ---------------------------- PFR #1 ----------------------------
    # ---------------------------------------------------------------- 
    
    mixer_air_2            = ct.IdealGasReactor(Fuel_1)
    outlet_PFR_1           = ct.MassFlowController(PFR_1, mixer_air_2, mdot=m_dot_air_1)                                                          
    t_res_PFR_1            = (Fuel_1.density * PFR_1.volume) / (m_dot_air_1)
    sim_PFR_1              = ct.ReactorNet([PFR_1])
    sim_PFR_1.advance(t_res_PFR_1)
    
    # ----------------------------------------------------------------
    # ------------------------- Mixing 2-Air -------------------------
    # ----------------------------------------------------------------     

    Air_2                  = Air
    Air_2.TPX              = T_stag_0, P_stag_0, dict_oxy
    rho_air_2              = Air_2.density    
    res_air_2              = ct.Reservoir(Air_2)
    PFR_2                  = ct.IdealGasConstPressureReactor(Fuel_1)  
    PFR_2.volume           = A_SZ*(L_SZ/3)
    inlet_air_2            = ct.MassFlowController(res_air_2, mixer_air_2, mdot=m_dot_air_SZ)
    m_dot_air_2            = m_dot_air_1 + m_dot_air_SZ
    outlet_air_2           = ct.MassFlowController(mixer_air_2, PFR_2, mdot=m_dot_air_2)
    sim_mixer_air_2        = ct.ReactorNet([mixer_air_2])
    sim_mixer_air_2.advance(t_mix)
    
    # ----------------------------------------------------------------
    # ---------------------------- PFR #2 ----------------------------
    # ---------------------------------------------------------------- 
    
    mixer_air_3            = ct.IdealGasReactor(Fuel_1)
    outlet_PFR_2           = ct.MassFlowController(PFR_2, mixer_air_3, mdot=m_dot_air_2)                                                          
    t_res_PFR_2            = (Fuel_1.density * PFR_2.volume) / (m_dot_air_2)
    sim_PFR_2              = ct.ReactorNet([PFR_2])
    sim_PFR_2.advance(t_res_PFR_2)    
    
    # ----------------------------------------------------------------
    # ------------------------- Mixing 3-Air -------------------------
    # ----------------------------------------------------------------     

    Air_3                  = Air
    Air_3.TPX              = T_stag_0, P_stag_0, dict_oxy
    rho_air_3              = Air_3.density    
    res_air_3              = ct.Reservoir(Air_3)
    PFR_3                  = ct.IdealGasConstPressureReactor(Fuel_1)  
    PFR_3.volume           = A_SZ*(L_SZ/3)
    inlet_air_3            = ct.MassFlowController(res_air_3, mixer_air_3, mdot=m_dot_air_SZ)
    m_dot_air_3            = m_dot_air_2 + m_dot_air_SZ
    outlet_air_3           = ct.MassFlowController(mixer_air_3, PFR_3, mdot=m_dot_air_3)
    sim_mixer_air_3        = ct.ReactorNet([mixer_air_3])
    sim_mixer_air_3.advance(t_mix)
    
    # ----------------------------------------------------------------
    # ---------------------------- PFR #3 ----------------------------
    # ---------------------------------------------------------------- 
                                                           
    t_res_PFR_3            = (Fuel_1.density * PFR_3.volume) / (m_dot_air_3)
    sim_PFR_3              = ct.ReactorNet([PFR_3])
    sim_PFR_3.advance(t_res_PFR_3)      
       
    # ----------------------------------------------------------------
    # --------------------- Additional computations ------------------
    # ----------------------------------------------------------------    
          
    m_dot_input_combustor  = m_dot_fuel + m_dot_air_id                    # [kg/s] Total mass flow rate entering a single combustor (air + fuel)
    Emission_Index = Fuel_1.Y * (m_dot_input_combustor)/m_dot_fuel 
    
    # Extract properties of combustor flow 
    a_out      = Fuel_1.sound_speed                                       # Speed of sound at PFR outlet
    rho_out    = Fuel_1.density                                           # density of the Fuel_1 in the combustor
    gamma      = Fuel_1.cp_mass / Fuel_1.cv_mass
    h          = Fuel_1.h                                                 # enthalpy
    vel_out    = (m_dot_input_combustor) / (rho_out * A_SZ)               # [m/s] Outlet velocity 
    M_out      = vel_out / a_out                                          # Outlet Mach number
    
    phi        = Fuel_1.equivalence_ratio(fuel = dict_fuel, oxidizer = dict_oxy) 
    
    # Stagnation temperature 
    T_stag_out = Fuel_1.T * (1 + 0.5 * (gamma - 1) * (M_out)**2)
    
    # stagnation pressure 
    P_stag_out = Fuel_1.P * (1 + 0.5 * (gamma - 1) * (M_out)**2)**(gamma / (gamma - 1))
    
    # Stagnation enthalpy 
    h_stag_out = T_stag_out  * Fuel_1.cp_mass
    
    # Fuel_1-to-air ratio (FAR)
    FAR      = m_dot_fuel / (m_dot_air_id)   
    
    return (Fuel_1, Emission_Index, T_stag_out, P_stag_out, h_stag_out, FAR) 

if __name__ == '__main__': 
    main()
    plt.show()