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
    m_dot_air_tot            = 40                                        # [kg/s] Air mass flow going through all combustors
    m_dot_air_TO_tot         = 44                                        # [kg/s] Air mass flow going through all combustors during TO 
    m_dot_fuel_tot           = m_dot_air_tot*FAR                          # [kg/s] Fuel mass flow going through all combustors
    m_dot_fuel_TO_tot        = m_dot_air_TO_tot*FAR_TO                    # [kg/s] Fuel mass flow going through all combustors during TO 
    N_comb                  = 9                                         # [-]    Number of can-annular combustors
    m_dot_air_id             = m_dot_air_tot/N_comb                       # [kg/s] Ideal Air mass flow inside each combustor, scaled inside each PSR to vary the Equivalence Ratio
    m_dot_air_TO_id          = m_dot_air_TO_tot/N_comb                    # [kg/s] Air mass flow inside each combustor during TO
    m_dot_fuel               = m_dot_fuel_tot/N_comb                      # [kg/s] Fuel mass flow inside each combustor
    m_dot_fuel_TO            = m_dot_fuel_TO_tot/N_comb                   # [kg/s] Fuel mass flow inside each combustor during TO    
    
    # Primary Zone Inputs
    N_PZ                    = 4                                         # [-]    Number of PSR (EVEN)
    V_PZ                    = 0.0023                                    # [m**3] Volume of the Primary Zone in a SINGLE combustor, must be split into the different PSRs       
    phi_PZ_des              = 1.77                                      # [-]    Design Equivalence Ratio
    S_PZ                    = 0.39                                      # [-]    Mixing parameter, used to define the Equivalence Ratio standard deviation         
    LHV_input_fuel          = 1                                         # [-]    FIX VALUE
    LHV_model_fuel          = 1                                         # [-]    FIX VALUE
    F_SC                    = LHV_input_fuel/LHV_model_fuel             # [-]    Fuel scaler, used to define the fraction of total air present in the combustor that enters the Primary Zone
    
    #Secondary Zone Inputs
    A_SZ                    = 0.15                                      # [m**2] Secondary Zone cross-sectional area
    L_SZ                    = 0.075                                     # [m]    Secondary Zone length
    l_SA_SM                 = 0.55                                      # [-]          
    l_SA_FM                 = 0.055                                     # [-]    
    l_DA_start              = 0.95                                      # [-]          
    l_DA_end                = 0.1                                       # [-] 
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
        gas, EI, T_stag_out, P_stag_out, h_stag_out, FAR = combustor(Fuel, Air, dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, m_dot_fuel, m_dot_fuel_TO, m_dot_air_id, m_dot_air_TO_id, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des, l_SA_SM, l_SA_FM, F_SC, A_SZ, L_SZ, l_DA_start, l_DA_end, f_SM)
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
 
def combustor(Fuel, Air, dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, m_dot_fuel, m_dot_fuel_TO, m_dot_air_id, m_dot_air_TO, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des, l_SA_SM, l_SA_FM, F_SC, A_SZ, L_SZ, l_DA_start, l_DA_end, f_SM):
    
    # ----------------------------------------------------------------
    # ---------------------- Initial Parameters ----------------------
    # ----------------------------------------------------------------    
      
    f_air_PZ               = (m_dot_fuel_TO*F_SC)/(phi_PZ_des*m_dot_air_TO*FAR_st) # Fraction of total air present in the combustor that enters the Primary Zone
    f_air_SZ               = 1 - f_air_PZ                                          # Fraction of total air present in the combustor that enters the Secondary Zone  
    m_dot_air              = f_air_PZ*m_dot_air_id                                 # Air mass flow going through the PSRs
    phi_sign               = ((m_dot_fuel*F_SC)/m_dot_air)/(FAR_st)                # Mean Equivalence Ratio
    sigma_phi              = S_PZ*phi_sign                                         # Standard deviation of the Equivalence Ratio    
    m_dot_air_PSR          = m_dot_air/N_PZ                                        # Air mass flow going through each PSR
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
    outlet_1               = ct.MassFlowController(PSR_1, mixer_12) 
    Y_fuel_1               = Fuel_1[comp_fuel].Y                                                  
    t_res_PSR_1            = V_PZ_PSR / (rho_1 * mass_flow_rate_1)
    sim_PSR_1              = ct.ReactorNet([PSR_1])
    sim_PSR_1.advance(t_res_PSR_1)
    #sim_psr.advance_to_steady_state()
    Y_fuel_1              = Fuel_1[comp_fuel].Y
    
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
    mixer_12               = ct.IdealGasReactor(Fuel_2)
    PSR_2                  = ct.IdealGasReactor(Fuel_2)
    PSR_2.volume           = V_PZ_PSR
    inlet_2                = ct.MassFlowController(upstream_2, PSR_2)
    inlet_2.mass_flow_rate = mass_flow_rate_2    
    outlet_2               = ct.MassFlowController(PSR_2, mixer_12) 
    Y_fuel_2               = Fuel_2[comp_fuel].Y                                                          
    t_res_PSR_2            = V_PZ_PSR / (rho_2 * mass_flow_rate_2)
    sim_PSR_2              = ct.ReactorNet([PSR_2])
    sim_PSR_2.advance(t_res_PSR_2)
    #sim_psr.advance_to_steady_state()
    Y_fuel_2              = Fuel_2[comp_fuel].Y     

    

    
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
    f_air_SA   = m_dot_fuel_TO/(phi_SZ_des*FAR_st*m_dot_air_TO)
    f_air_DA   = 1 - f_air_PZ - f_air_SA
    
    beta_SA_SM = (f_air_SA*f_SM*m_dot_air)/(l_SA_SM * L_SZ)
    beta_SA_FM = (f_air_SA*f_FM*m_dot_air)/(l_SA_FM * L_SZ)
    beta_DA    = (f_air_DA*m_dot_air)/((l_DA_end - l_DA_start)*L_SZ)
    
    sim = ct.ReactorNet([PFR_1])
    sim.rtol = 1e-4  # Relax the relative tolerance further
    sim.atol = 1e-8  # Adjust absolute tolerance to control smaller values
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
            
        air_mass_flow = max(0, beta_air_in * air.density * dz)
        mfc_air.mass_flow_rate = air_mass_flow
        outlet = ct.MassFlowController(PFR_1, mixer_2)
        outlet.mass_flow_rate = m_dot_air + m_dot_fuel + air_mass_flow
        
        try:
            sim.advance(z + dz / 100)
            # After advancing the simulation step:
            if Fuel.T < 300:  # Set a reasonable minimum temperature (e.g., 300 K)
                print(f"Warning: Unphysical temperature detected: T = {Fuel.T:.2f} K at position {z:.2f} m")
                break  # Stop simulation or try a corrective measure            
        except ct.CanteraError as e:
            print(f"Warning: Error during simulation at position {z:.2f} m: {e}")
        
        # Print results for this segment
        print(f"Position: {z:.2f} m, Temperature: {Fuel.T:.2f} K, Pressure: {Fuel.P:.2f} Pa")    
        
    ## ----------------------------------------------------------------
    ## ---------------------------- PFR #2 ----------------------------
    ## ----------------------------------------------------------------
    
    sim = ct.ReactorNet([PFR_2])
    sim.rtol = 1e-4  # Relax the relative tolerance further
    sim.atol = 1e-8  # Adjust absolute tolerance to control smaller values 
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
            
        air_mass_flow = max(0, beta_air_in * air.density * dz)
        mfc_air.mass_flow_rate = air_mass_flow
        outlet = ct.MassFlowController(PFR_2, mixer_2)
        outlet.mass_flow_rate = m_dot_air + m_dot_fuel + air_mass_flow
        
        try:
            sim.advance(z + dz / 100)
            # After advancing the simulation step:
            if Fuel.T < 300:  # Set a reasonable minimum temperature (e.g., 300 K)
                print(f"Warning: Unphysical temperature detected: T = {Fuel.T:.2f} K at position {z:.2f} m")
                break  # Stop simulation or try a corrective measure            
        except ct.CanteraError as e:
            print(f"Warning: Error during simulation at position {z:.2f} m: {e}")
        
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
    Emission_Index = Fuel.Y * (m_dot_fuel + m_dot_air)/m_dot_fuel 
    
    # Extract properties of combustor flow 
    a_out      = Fuel.sound_speed  # Speed of sound at PFR outlet
    rho_out    = Fuel.density # density of the Fuel in the combustor
    gamma      = Fuel.cp_mass / Fuel.cv_mass
    h          = Fuel.h # enthalpy
    vel_out    = (m_dot_fuel + m_dot_air) / (rho_out * A_SZ)  # Outlet velocity (m/s)  
    M_out      = vel_out / a_out  # Outlet Mach number
    
    phi        = Fuel.equivalence_ratio(fuel = dict_fuel, oxidizer = dict_oxy) 
    
    # Stagnation temperature 
    T_stag_out = Fuel.T * (1 + 0.5 * (gamma - 1) * (M_out)**2)
    
    # stagnation pressure 
    P_stag_out = Fuel.P * (1 + 0.5 * (gamma - 1) * (M_out)**2)**(gamma / (gamma - 1))
    
    # Stagnation enthalpy 
    h_stag_out = T_stag_out  * Fuel.cp_mass
    
    # Fuel-to-air ratio (FAR)
    FAR      = m_dot_fuel / (m_dot_air)   
    
    return (Fuel, Emission_Index, T_stag_out, P_stag_out, h_stag_out, FAR) 

if __name__ == '__main__': 
    main()
    plt.show()