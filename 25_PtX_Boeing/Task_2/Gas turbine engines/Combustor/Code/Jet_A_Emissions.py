import cantera              as ct
import numpy                as np
import pandas               as pd
import matplotlib.pyplot    as plt
import time 

def main():

    ti                      = time.time()                                           # [s]       Define the initial simulation time
                                                                                    
    # ------------------------------------------------------------------------------              
    # ------------------------------ Combustor Inputs ------------------------------              
    # ------------------------------------------------------------------------------              
                                                                                    
    # General Inputs                                                                
    full_kin_mech           = False                                                 # [-]       True (around 300 s): Computes EI_CO2, EI_CO, EI_H2O, EI_NO2, EI_NO, EI_CSOOT; False (around ): Computes EI_CO2, EI_CO, EI_H2O
    T_stag_0                = 800                                                   # [K]       Stagnation Temperature entering all combustors
    P_stag_0                = 2000000                                               # [Pa]      Stagnation Pressure entering all combustors
    FAR                     = 0.02                                                  # [-]       Fuel-to-Air ratio
    FAR_TO                  = 0.0275                                                # [-]       Fuel-to-Air ratio during TO 
    FAR_st                  = 0.068                                                 # [-]       Stoichiometric Fuel-to-Air ratio
    m_dot_air_tot           = 40                                                    # [kg/s]    Air mass flow going through all combustors
    m_dot_air_TO_tot        = 44                                                    # [kg/s]    Air mass flow going through all combustors during TO 
    m_dot_fuel_tot          = m_dot_air_tot*FAR                                     # [kg/s]    Fuel mass flow going through all combustors
    m_dot_fuel_TO_tot       = m_dot_air_TO_tot*FAR_TO                               # [kg/s]    Fuel mass flow going through all combustors during TO 
    N_comb                  = 10                                                    # [-]       Number of can-annular combustors
    m_dot_air_id            = m_dot_air_tot/N_comb                                  # [kg/s]    Ideal Air mass flow inside each combustor, scaled inside each PSR to vary the Equivalence Ratio
    m_dot_air_TO            = m_dot_air_TO_tot/N_comb                               # [kg/s]    Air mass flow inside each combustor during TO
    m_dot_fuel              = m_dot_fuel_tot/N_comb                                 # [kg/s]    Fuel mass flow inside each combustor
    m_dot_fuel_TO           = m_dot_fuel_TO_tot/N_comb                              # [kg/s]    Fuel mass flow inside each combustor during TO    
                                                                                                
    # Primary Zone Inputs                                                                       
    N_PZ                    = 8                                                     # [-]       Number of PSR (EVEN, must match the number of PSR below)
    V_PZ                    = 0.0023                                                # [m**3]    Volume of the Primary Zone in a SINGLE combustor, must be split into the different PSRs       
    phi_PZ_des              = 0.7                                                   # [-]       Primary Zone Design Equivalence Ratio
    S_PZ                    = 0.4                                                   # [-]       Mixing parameter, used to define the Equivalence Ratio standard deviation  
    F_SC                    = 0.8                                                   # [-]       Fuel scaler, used to define the fraction of total air present in the combustor that enters the Primary Zone
                                                                                                
    #Secondary Zone Inputs                                                                      
    A_SZ                    = 0.15                                                  # [m**2]    Secondary Zone cross-sectional area
    L_SZ                    = 0.075                                                 # [m]       Secondary Zone length  
    phi_SZ_des_1            = 0.7                                                   # [-]       Design Equivalence Ratio for PFR #1
    phi_SZ_des_2            = 0.2                                                   # [-]       Design Equivalence Ratio for PFR #2
    phi_SZ_des_3            = 0.2                                                   # [-]       Design Equivalence Ratio for PFR #3
                                                                                    
    if full_kin_mech:                                                               
        dict_fuel           = {'NC10H22':0.16449, 'NC12H26':0.34308,'NC16H34':0.10335, 'IC8H18':0.08630,'NC7H14':0.07945, 'C6H5C2H5': 0.07348,'C6H5C4H9': 0.05812, 'C10H7CH3': 0.10972} # [-]       Fuel species and corresponding mole fractions for full fuel model
    else:                                                                                       
        dict_fuel           = {'N-C12H26':0.6, 'A1CH3':0.2, 'A1':0.2}               # [-]       Fuel species and corresponding mole fractions for surrogate fuel model 
                                                                                                
    dict_oxy                = {'O2':0.2095, 'N2':0.7809, 'AR':0.0096}               # [-]       Air species and corresponding mole fractions     
                                                                                                
    if full_kin_mech:                                                                           
        list_sp             = ['CO', 'CO2', 'H2O', 'NO', 'NO2', 'CSOLID']           # [-]       Fuel species for Emission Index analysis
    else:                                                                                       
        list_sp             = ['CO2', 'CO', 'H2O']                                  # [-]       Fuel species for Emission Index analysis
                                                                                                
    col_names = ['Tout(K)', 'T_stag_out','P_stag_out', 'h_stag_out'] + ['X_' +str(sp) for sp in list_sp] + ['Y_' +str(sp) for sp in list_sp] + ['EI_' +str(sp) for sp in list_sp] # [-]       Define output variables 
    df                      = pd.DataFrame(columns=col_names)                       # [-]       Assign output variables space to df
                                                                                    
    for n in range(1):                                                              
        gas, EI, T_stag_out, P_stag_out, h_stag_out = combustor(full_kin_mech, dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO, FAR_st, m_dot_fuel, m_dot_fuel_TO, m_dot_air_id, m_dot_air_TO, N_PZ, V_PZ, phi_PZ_des, S_PZ, phi_SZ_des_1, phi_SZ_des_2, phi_SZ_des_3, F_SC, A_SZ, L_SZ)                 # [-]       Run combustor function
                                                                                                
        sp_idx              = [gas.species_index(sp) for sp in list_sp]             # [-]       Retrieve the species index
        data_n              = [gas.T, T_stag_out, P_stag_out, h_stag_out] + list(gas.X[sp_idx]) + list(gas.Y[sp_idx]) + list(EI[sp_idx])              # [-]       Assign output variables  
        df.loc[n]           = data_n                                                # [-]       Assign output variables to df 
                                                                                                
    print(df['EI_CO2'])                                                             # [-]       Print the value of EI_CO2
    print(df['EI_CO'])                                                              # [-]       Print the value of EI_CO
    print(df['EI_H2O'])                                                             # [-]       Print the value of EI_H2O
                                                                                                
    if full_kin_mech:                                                                           
        print(df['EI_NO'])                                                          # [-]       Print the value of EI_NO
        print(df['EI_NO2'])                                                         # [-]       Print the value of EI_NO2
        print(df['EI_CSOLID'])                                                      # [-]       Print the value of EI_CSOLID
                                                                                                
    tf                      = time.time()                                           # [s]       Define the final simulation time
    elapsed_time            = round((tf-ti),2)                                      # [s]       Compute the total simulation time
                                                                                    
    print('Simulation Time: ' + str(elapsed_time) + ' seconds per timestep')        # [-]       Print the value of total simulation time
    
    return 
 
def combustor(full_kin_mech, dict_fuel, dict_oxy, T_stag_0, P_stag_0, FAR, FAR_TO,FAR_st, m_dot_fuel, m_dot_fuel_TO, m_dot_air_id, m_dot_air_TO, N_PZ,V_PZ, phi_PZ_des, S_PZ, phi_SZ_des_1, phi_SZ_des_2, phi_SZ_des_3,F_SC, A_SZ, L_SZ):
    
    # ------------------------------------------------------------------------------
    # ----------------------------- Initial Parameters -----------------------------
    # ------------------------------------------------------------------------------    
      
    f_air_PZ                = (m_dot_fuel_TO*F_SC)/(phi_PZ_des*m_dot_air_TO*FAR_st) # [-]       Fraction of total air present in the combustor that enters the Primary Zone
    f_air_SZ                = 1 - f_air_PZ                                          # [-]       Fraction of total air present in the combustor that enters the Secondary Zone  
    m_dot_air_PZ            = f_air_PZ*m_dot_air_id                                 # [kg/s]    Air mass flow going through the Primary Zone
    m_dot_air_SZ            = (f_air_SZ*m_dot_air_id)/3                             # [kg/s]    Air mass flow going through each dilution air inlet (3 inlets)
    phi_sign                = ((m_dot_fuel*F_SC)/m_dot_air_PZ)/(FAR_st)             # [-]       Mean Equivalence Ratio
    sigma_phi               = S_PZ*phi_sign                                         # [-]       Standard deviation of the Equivalence Ratio    
    m_dot_air_PSR           = m_dot_air_PZ                                          # [kg/s]    Air mass flow going through each PSR
    m_dot_fuel_PSR          = m_dot_fuel                                            # [kg/s]    Fuel mass flow going through each PSR
    V_PZ_PSR                = V_PZ/N_PZ                                             # [m**3]    Volume of each PSR
    phi_PSR                 = np.linspace(0.001, 2*phi_sign, N_PZ)                  # [-]       Distribution of Equivalence Ratio through the PSRs
    Delta_phi               = np.abs(phi_PSR[0] - phi_PSR[1])                       # [-]       Difference between two subsequent Equivalence Ratios
    comp_fuel               = list(dict_fuel.keys())                                # [-]       Fuel components
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PSR #1 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    f_PZ_1                  = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[0] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi # [-]       Fraction of mass flow entering the PSR at the PSR equivalence ratio 
    if full_kin_mech:                                                                           
        Fuel_1              = ct.Solution('chem.yaml')                              # [-]       Import full fuel kinematic mechanism
    else:                                                                                       
        Fuel_1              = ct.Solution('JetFuelSurrogate.yaml')                  # [-]       Import surrogate fuel kinematic mechanism
    Fuel_1.TP               = T_stag_0, P_stag_0                                    # [-]       Set the fuel temperature and pressure
    Fuel_1.set_equivalence_ratio(phi_PSR[0], fuel=dict_fuel, oxidizer=dict_oxy)     # [-]       Set the euivalence ratio inside the PSR
    Fuel_1.equilibrate('HP')                                                        # [-]       Fix the specific enthalpy and pressure 
    rho_1                   = Fuel_1.density                                        # [kg/m**3] Fuel density
    m_dot_air_1             = m_dot_air_PSR * f_PZ_1                                # [kg/s]    Air mass flow inside the PSR    
    m_dot_fuel_1            = m_dot_fuel_PSR * f_PZ_1                               # [kg/s]    Fuel mass flow inside the PSR
    mass_flow_rate_1        = m_dot_fuel_1 + m_dot_air_1                            # [kg/s]    Total mass flow inside the PSR  
    upstream_1              = ct.Reservoir(Fuel_1)                                  # [-]       Create a resevoir for the Fuel upstream of the PSR
    mixer_12                = ct.IdealGasReactor(Fuel_1)                            # [-]       Create the reactor for the downstream mixer
    PSR_1                   = ct.IdealGasReactor(Fuel_1)                            # [-]       Create the reactor for the PSR 
    PSR_1.volume            = V_PZ_PSR                                              # [m**3]    Set the PSR volume
    inlet_1                 = ct.MassFlowController(upstream_1, PSR_1)              # [-]       Connect the upstream resevoir with the PSR 
    inlet_1.mass_flow_rate  = mass_flow_rate_1                                      # [kg/s]    Prescribe the inlet mass flow rate 
    outlet_1                = ct.MassFlowController(PSR_1, mixer_12, mdot=mass_flow_rate_1)          # [-]       Connect the PSR with the downstream mixer with the same mass flow rate                                     
    t_res_PSR_1             = (rho_1 * V_PZ_PSR) / (mass_flow_rate_1)               # [s]       Compute the residence time inside the PSR
    sim_PSR_1               = ct.ReactorNet([PSR_1])                                # [-]       Set the PSR simulation
    sim_PSR_1.advance(t_res_PSR_1)                                                  # [-]       Run the simulation until the residence time is reached
    Y_fuel_1                = Fuel_1[comp_fuel].Y                                   # [-]       Store mass fractions inside the PSR
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PSR #2 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    f_PZ_2                  = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[1] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi # [-]       Fraction of mass flow entering the PSR at the PSR equivalence ratio 
    if full_kin_mech:                                                                           
        Fuel_2              = ct.Solution('chem.yaml')                              # [-]       Import full fuel kinematic mechanism
    else:                                                                                       
        Fuel_2              = ct.Solution('JetFuelSurrogate.yaml')                  # [-]       Import surrogate fuel kinematic mechanism
    Fuel_2.TP               = T_stag_0, P_stag_0                                    # [-]       Set the fuel temperature and pressure
    Fuel_2.set_equivalence_ratio(phi_PSR[1], fuel=dict_fuel, oxidizer=dict_oxy)     # [-]       Set the euivalence ratio inside the PSR
    Fuel_2.equilibrate('HP')                                                        # [-]       Fix the specific enthalpy and pressure 
    rho_2                   = Fuel_2.density                                        # [kg/m**3] Fuel density
    m_dot_air_2             = m_dot_air_PSR * f_PZ_2                                # [kg/s]    Air mass flow inside the PSR    
    m_dot_fuel_2            = m_dot_fuel_PSR * f_PZ_2                               # [kg/s]    Fuel mass flow inside the PSR
    mass_flow_rate_2        = m_dot_fuel_2 + m_dot_air_2                            # [kg/s]    Total mass flow inside the PSR  
    upstream_2              = ct.Reservoir(Fuel_2)                                  # [-]       Create a resevoir for the Fuel upstream of the PSR
    PSR_2                   = ct.IdealGasReactor(Fuel_2)                            # [-]       Create the reactor for the PSR 
    PSR_2.volume            = V_PZ_PSR                                              # [m**3]    Set the PSR volume
    inlet_2                 = ct.MassFlowController(upstream_2, PSR_2)              # [-]       Connect the upstream resevoir with the PSR 
    inlet_2.mass_flow_rate  = mass_flow_rate_2                                      # [kg/s]    Prescribe the inlet mass flow rate 
    outlet_2                = ct.MassFlowController(PSR_2, mixer_12, mdot=mass_flow_rate_2)          # [-]       Connect the PSR with the downstream mixer with the same mass flow rate               
    t_res_PSR_2             = (rho_2 * V_PZ_PSR) / (mass_flow_rate_2)               # [s]       Compute the residence time inside the PSR               
    sim_PSR_2               = ct.ReactorNet([PSR_2])                                # [-]       Set the PSR simulation
    sim_PSR_2.advance(t_res_PSR_2)                                                  # [-]       Run the simulation until the residence time is reached
    Y_fuel_2                = Fuel_2[comp_fuel].Y                                   # [-]       Store mass fractions inside the PSR
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PSR #3 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    f_PZ_3                  = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[2] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi # [-]       Fraction of mass flow entering the PSR at the PSR equivalence ratio 
    if full_kin_mech:                                                                           
        Fuel_3              = ct.Solution('chem.yaml')                              # [-]       Import full fuel kinematic mechanism
    else:                                                                                       
        Fuel_3              = ct.Solution('JetFuelSurrogate.yaml')                  # [-]       Import surrogate fuel kinematic mechanism
    Fuel_3.TP               = T_stag_0, P_stag_0                                    # [-]       Set the fuel temperature and pressure
    Fuel_3.set_equivalence_ratio(phi_PSR[2], fuel=dict_fuel, oxidizer=dict_oxy)     # [-]       Set the euivalence ratio inside the PSR
    Fuel_3.equilibrate('HP')                                                        # [-]       Fix the specific enthalpy and pressure 
    rho_3                   = Fuel_3.density                                        # [kg/m**3] Fuel density
    m_dot_air_3             = m_dot_air_PSR * f_PZ_3                                # [kg/s]    Air mass flow inside the PSR    
    m_dot_fuel_3            = m_dot_fuel_PSR * f_PZ_3                               # [kg/s]    Fuel mass flow inside the PSR
    mass_flow_rate_3        = m_dot_fuel_3 + m_dot_air_3                            # [kg/s]    Total mass flow inside the PSR  
    upstream_3              = ct.Reservoir(Fuel_3)                                  # [-]       Create a resevoir for the Fuel upstream of the PSR
    mixer_34                = ct.IdealGasReactor(Fuel_3)                            # [-]       Create the reactor for the downstream mixer
    PSR_3                   = ct.IdealGasReactor(Fuel_3)                            # [-]       Create the reactor for the PSR 
    PSR_3.volume            = V_PZ_PSR                                              # [m**3]    Set the PSR volume
    inlet_3                 = ct.MassFlowController(upstream_3, PSR_3)              # [-]       Connect the upstream resevoir with the PSR 
    inlet_3.mass_flow_rate  = mass_flow_rate_3                                      # [kg/s]    Prescribe the inlet mass flow rate 
    outlet_3                = ct.MassFlowController(PSR_3, mixer_34,mdot=mass_flow_rate_3)          # [-]       Connect the PSR with the downstream mixer with the same mass flow rate                                                 
    t_res_PSR_3             = (rho_3 * V_PZ_PSR) / (mass_flow_rate_3)               # [s]       Compute the residence time inside the PSR
    sim_PSR_3               = ct.ReactorNet([PSR_3])                                # [-]       Set the PSR simulation
    sim_PSR_3.advance(t_res_PSR_3)                                                  # [-]       Run the simulation until the residence time is reached
    Y_fuel_3                = Fuel_3[comp_fuel].Y                                   # [-]       Store mass fractions inside the PSR
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PSR #4 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    f_PZ_4                  = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[3] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi # [-]       Fraction of mass flow entering the PSR at the PSR equivalence ratio 
    if full_kin_mech:                                                                           
        Fuel_4              = ct.Solution('chem.yaml')                              # [-]       Import full fuel kinematic mechanism
    else:                                                                                       
        Fuel_4              = ct.Solution('JetFuelSurrogate.yaml')                  # [-]       Import surrogate fuel kinematic mechanism
    Fuel_4.TP               = T_stag_0, P_stag_0                                    # [-]       Set the fuel temperature and pressure
    Fuel_4.set_equivalence_ratio(phi_PSR[3], fuel=dict_fuel, oxidizer=dict_oxy)     # [-]       Set the euivalence ratio inside the PSR
    Fuel_4.equilibrate('HP')                                                        # [-]       Fix the specific enthalpy and pressure 
    rho_4                   = Fuel_4.density                                        # [kg/m**3] Fuel density
    m_dot_air_4             = m_dot_air_PSR * f_PZ_4                                # [kg/s]    Air mass flow inside the PSR    
    m_dot_fuel_4            = m_dot_fuel_PSR * f_PZ_4                               # [kg/s]    Fuel mass flow inside the PSR
    mass_flow_rate_4        = m_dot_fuel_4 + m_dot_air_4                            # [kg/s]    Total mass flow inside the PSR  
    upstream_4              = ct.Reservoir(Fuel_4)                                  # [-]       Create a resevoir for the Fuel upstream of the PSR
    PSR_4                   = ct.IdealGasReactor(Fuel_4)                            # [-]       Create the reactor for the PSR 
    PSR_4.volume            = V_PZ_PSR                                              # [m**3]    Set the PSR volume
    inlet_4                 = ct.MassFlowController(upstream_4, PSR_4)              # [-]       Connect the upstream resevoir with the PSR 
    inlet_4.mass_flow_rate  = mass_flow_rate_4                                      # [kg/s]    Prescribe the inlet mass flow rate 
    outlet_4                = ct.MassFlowController(PSR_4, mixer_34, mdot=mass_flow_rate_4)          # [-]       Connect the PSR with the downstream mixer with the same mass flow rate                                                   
    t_res_PSR_4             = (rho_4 * V_PZ_PSR) / (mass_flow_rate_4)               # [s]       Compute the residence time inside the PSR               
    sim_PSR_4               = ct.ReactorNet([PSR_4])                                # [-]       Set the PSR simulation
    sim_PSR_4.advance(t_res_PSR_4)                                                  # [-]       Run the simulation until the residence time is reached
    Y_fuel_4                = Fuel_4[comp_fuel].Y                                   # [-]       Store mass fractions inside the PSR
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PSR #5 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    f_PZ_5                  = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[4] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi # [-]       Fraction of mass flow entering the PSR at the PSR equivalence ratio 
    if full_kin_mech:                                                                           
        Fuel_5              = ct.Solution('chem.yaml')                              # [-]       Import full fuel kinematic mechanism
    else:                                                                                       
        Fuel_5              = ct.Solution('JetFuelSurrogate.yaml')                  # [-]       Import surrogate fuel kinematic mechanism
    Fuel_5.TP               = T_stag_0, P_stag_0                                    # [-]       Set the fuel temperature and pressure
    Fuel_5.set_equivalence_ratio(phi_PSR[4], fuel=dict_fuel, oxidizer=dict_oxy)     # [-]       Set the euivalence ratio inside the PSR
    Fuel_5.equilibrate('HP')                                                        # [-]       Fix the specific enthalpy and pressure 
    rho_5                   = Fuel_5.density                                        # [kg/m**3] Fuel density
    m_dot_air_5             = m_dot_air_PSR * f_PZ_5                                # [kg/s]    Air mass flow inside the PSR    
    m_dot_fuel_5            = m_dot_fuel_PSR * f_PZ_5                               # [kg/s]    Fuel mass flow inside the PSR
    mass_flow_rate_5        = m_dot_fuel_5 + m_dot_air_5                            # [kg/s]    Total mass flow inside the PSR  
    upstream_5              = ct.Reservoir(Fuel_5)                                  # [-]       Create a resevoir for the Fuel upstream of the PSR
    mixer_56                = ct.IdealGasReactor(Fuel_5)                            # [-]       Create the reactor for the downstream mixer
    PSR_5                   = ct.IdealGasReactor(Fuel_5)                            # [-]       Create the reactor for the PSR 
    PSR_5.volume            = V_PZ_PSR                                              # [m**3]    Set the PSR volume
    inlet_5                 = ct.MassFlowController(upstream_5, PSR_5)              # [-]       Connect the upstream resevoir with the PSR 
    inlet_5.mass_flow_rate  = mass_flow_rate_5                                      # [kg/s]    Prescribe the inlet mass flow rate 
    outlet_5                = ct.MassFlowController(PSR_5, mixer_56, mdot=mass_flow_rate_5)          # [-]       Connect the PSR with the downstream mixer with the same mass flow rate                                                  
    t_res_PSR_5             = (rho_5 * V_PZ_PSR) / (mass_flow_rate_5)               # [s]       Compute the residence time inside the PSR
    sim_PSR_5               = ct.ReactorNet([PSR_5])                                # [-]       Set the PSR simulation
    sim_PSR_5.advance(t_res_PSR_5)                                                  # [-]       Run the simulation until the residence time is reached
    Y_fuel_5                = Fuel_5[comp_fuel].Y                                   # [-]       Store mass fractions inside the PSR
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PSR #6 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    f_PZ_6                  = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[5] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi # [-]       Fraction of mass flow entering the PSR at the PSR equivalence ratio 
    if full_kin_mech:                                                                           
        Fuel_6              = ct.Solution('chem.yaml')                              # [-]       Import full fuel kinematic mechanism
    else:                                                                                       
        Fuel_6              = ct.Solution('JetFuelSurrogate.yaml')                  # [-]       Import surrogate fuel kinematic mechanism
    Fuel_6.TP               = T_stag_0, P_stag_0                                    # [-]       Set the fuel temperature and pressure
    Fuel_6.set_equivalence_ratio(phi_PSR[5], fuel=dict_fuel, oxidizer=dict_oxy)     # [-]       Set the euivalence ratio inside the PSR
    Fuel_6.equilibrate('HP')                                                        # [-]       Fix the specific enthalpy and pressure 
    rho_6                   = Fuel_6.density                                        # [kg/m**3] Fuel density
    m_dot_air_6             = m_dot_air_PSR * f_PZ_6                                # [kg/s]    Air mass flow inside the PSR    
    m_dot_fuel_6            = m_dot_fuel_PSR * f_PZ_6                               # [kg/s]    Fuel mass flow inside the PSR
    mass_flow_rate_6        = m_dot_fuel_6 + m_dot_air_6                            # [kg/s]    Total mass flow inside the PSR  
    upstream_6              = ct.Reservoir(Fuel_6)                                  # [-]       Create a resevoir for the Fuel upstream of the PSR
    PSR_6                   = ct.IdealGasReactor(Fuel_6)                            # [-]       Create the reactor for the PSR 
    PSR_6.volume            = V_PZ_PSR                                              # [m**3]    Set the PSR volume
    inlet_6                 = ct.MassFlowController(upstream_6, PSR_6)              # [-]       Connect the upstream resevoir with the PSR 
    inlet_6.mass_flow_rate  = mass_flow_rate_6                                      # [kg/s]    Prescribe the inlet mass flow rate 
    outlet_6                = ct.MassFlowController(PSR_6, mixer_56, mdot=mass_flow_rate_6)          # [-]       Connect the PSR with the downstream mixer with the same mass flow rate                                                  
    t_res_PSR_6             = (rho_6 * V_PZ_PSR) / (mass_flow_rate_6)               # [s]       Compute the residence time inside the PSR               
    sim_PSR_6               = ct.ReactorNet([PSR_6])                                # [-]       Set the PSR simulation
    sim_PSR_6.advance(t_res_PSR_6)                                                  # [-]       Run the simulation until the residence time is reached
    Y_fuel_6                = Fuel_6[comp_fuel].Y                                   # [-]       Store mass fractions inside the PSR
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PSR #7 -----------------------------------
    # ------------------------------------------------------------------------------ 

    f_PZ_7                  = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[6] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi # [-]       Fraction of mass flow entering the PSR at the PSR equivalence ratio 
    if full_kin_mech:                                                                           
        Fuel_7              = ct.Solution('chem.yaml')                              # [-]       Import full fuel kinematic mechanism
    else:                                                                                       
        Fuel_7              = ct.Solution('JetFuelSurrogate.yaml')                  # [-]       Import surrogate fuel kinematic mechanism
    Fuel_7.TP               = T_stag_0, P_stag_0                                    # [-]       Set the fuel temperature and pressure
    Fuel_7.set_equivalence_ratio(phi_PSR[6], fuel=dict_fuel, oxidizer=dict_oxy)     # [-]       Set the euivalence ratio inside the PSR
    Fuel_7.equilibrate('HP')                                                        # [-]       Fix the specific enthalpy and pressure 
    rho_7                   = Fuel_7.density                                        # [kg/m**3] Fuel density
    m_dot_air_7             = m_dot_air_PSR * f_PZ_7                                # [kg/s]    Air mass flow inside the PSR    
    m_dot_fuel_7            = m_dot_fuel_PSR * f_PZ_7                               # [kg/s]    Fuel mass flow inside the PSR
    mass_flow_rate_7        = m_dot_fuel_7 + m_dot_air_7                            # [kg/s]    Total mass flow inside the PSR  
    upstream_7              = ct.Reservoir(Fuel_7)                                  # [-]       Create a resevoir for the Fuel upstream of the PSR
    mixer_78                = ct.IdealGasReactor(Fuel_7)                            # [-]       Create the reactor for the downstream mixer
    PSR_7                   = ct.IdealGasReactor(Fuel_7)                            # [-]       Create the reactor for the PSR 
    PSR_7.volume            = V_PZ_PSR                                              # [m**3]    Set the PSR volume
    inlet_7                 = ct.MassFlowController(upstream_7, PSR_7)              # [-]       Connect the upstream resevoir with the PSR 
    inlet_7.mass_flow_rate  = mass_flow_rate_7                                      # [kg/s]    Prescribe the inlet mass flow rate 
    outlet_7                = ct.MassFlowController(PSR_7, mixer_78,mdot=mass_flow_rate_7)          # [-]       Connect the PSR with the downstream mixer with the same mass flow rate                                                
    t_res_PSR_7             = (rho_7 * V_PZ_PSR) / (mass_flow_rate_7)               # [s]       Compute the residence time inside the PSR
    sim_PSR_7               = ct.ReactorNet([PSR_7])                                # [-]       Set the PSR simulation
    sim_PSR_7.advance(t_res_PSR_7)                                                  # [-]       Run the simulation until the residence time is reached
    Y_fuel_7                = Fuel_7[comp_fuel].Y                                   # [-]       Store mass fractions inside the PSR
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PSR #8 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    f_PZ_8                  = (1 / (np.sqrt(2 * np.pi) * sigma_phi)) * np.exp((-(phi_PSR[7] - phi_sign) ** 2) / (2 * sigma_phi ** 2)) * Delta_phi # [-]       Fraction of mass flow entering the PSR at the PSR equivalence ratio 
    if full_kin_mech:                                                                           
        Fuel_8              = ct.Solution('chem.yaml')                              # [-]       Import full fuel kinematic mechanism
    else:                                                                                       
        Fuel_8              = ct.Solution('JetFuelSurrogate.yaml')                  # [-]       Import surrogate fuel kinematic mechanism
    Fuel_8.TP               = T_stag_0, P_stag_0                                    # [-]       Set the fuel temperature and pressure
    Fuel_8.set_equivalence_ratio(phi_PSR[7], fuel=dict_fuel, oxidizer=dict_oxy)     # [-]       Set the euivalence ratio inside the PSR
    Fuel_8.equilibrate('HP')                                                        # [-]       Fix the specific enthalpy and pressure 
    rho_8                   = Fuel_8.density                                        # [kg/m**3] Fuel density
    m_dot_air_8             = m_dot_air_PSR * f_PZ_8                                # [kg/s]    Air mass flow inside the PSR    
    m_dot_fuel_8            = m_dot_fuel_PSR * f_PZ_8                               # [kg/s]    Fuel mass flow inside the PSR
    mass_flow_rate_8        = m_dot_fuel_8 + m_dot_air_8                            # [kg/s]    Total mass flow inside the PSR  
    upstream_8              = ct.Reservoir(Fuel_8)                                  # [-]       Create a resevoir for the Fuel upstream of the PSR
    PSR_8                   = ct.IdealGasReactor(Fuel_8)                            # [-]       Create the reactor for the PSR 
    PSR_8.volume            = V_PZ_PSR                                              # [m**3]    Set the PSR volume
    inlet_8                 = ct.MassFlowController(upstream_8, PSR_8)              # [-]       Connect the upstream resevoir with the PSR 
    inlet_8.mass_flow_rate  = mass_flow_rate_8                                      # [kg/s]    Prescribe the inlet mass flow rate 
    outlet_8                = ct.MassFlowController(PSR_8, mixer_78,mdot=mass_flow_rate_8)          # [-]       Connect the PSR with the downstream mixer with the same mass flow rate                                                
    t_res_PSR_8             = (rho_8 * V_PZ_PSR) / (mass_flow_rate_8)               # [s]       Compute the residence time inside the PSR               
    sim_PSR_8               = ct.ReactorNet([PSR_8])                                # [-]       Set the PSR simulation
    sim_PSR_8.advance(t_res_PSR_8)                                                  # [-]       Run the simulation until the residence time is reached
    Y_fuel_8                = Fuel_8[comp_fuel].Y                                   # [-]       Store mass fractions inside the PSR
    
    # ------------------------------------------------------------------------------
    # --------------------------------- Mixing 1-2 ---------------------------------
    # ------------------------------------------------------------------------------     
    
    mixer_1234              = ct.IdealGasReactor(Fuel_1)                            # [-]       Create the reactor for the downstream mixer
    outlet_12               = ct.MassFlowController(mixer_12, mixer_1234, mdot = (mass_flow_rate_1 + mass_flow_rate_2)) # [-]       Connect the mixer with the downstream mixer     
    sim_mixer_12            = ct.ReactorNet([mixer_12])                             # [-]       Set the mixer simulation  
    sim_mixer_12.advance_to_steady_state()                                          # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # --------------------------------- Mixing 3-4 ---------------------------------
    # ------------------------------------------------------------------------------     
    
    outlet_34               = ct.MassFlowController(mixer_34, mixer_1234, mdot = (mass_flow_rate_3 + mass_flow_rate_4)) # [-]       Connect the mixer with the downstream mixer
    sim_mixer_34            = ct.ReactorNet([mixer_34])                             # [-]       Set the mixer simulation  
    sim_mixer_34.advance_to_steady_state()                                          # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # --------------------------------- Mixing 5-6 ---------------------------------
    # ------------------------------------------------------------------------------     
    
    mixer_5678              = ct.IdealGasReactor(Fuel_5)                            # [-]       Create the reactor for the downstream mixer
    outlet_56               = ct.MassFlowController(mixer_56, mixer_5678, mdot = (mass_flow_rate_5 + mass_flow_rate_6)) # [-]       Connect the mixer with the downstream mixer     
    sim_mixer_56            = ct.ReactorNet([mixer_56])                             # [-]       Set the mixer simulation  
    sim_mixer_56.advance_to_steady_state()                                          # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # --------------------------------- Mixing 7-8 ---------------------------------
    # ------------------------------------------------------------------------------     
    
    outlet_78               = ct.MassFlowController(mixer_78, mixer_5678, mdot = (mass_flow_rate_7 + mass_flow_rate_8)) # [-]       Connect the mixer with the downstream mixer
    sim_mixer_78            = ct.ReactorNet([mixer_78])                             # [-]       Set the mixer simulation  
    sim_mixer_78.advance_to_steady_state()                                          # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # ------------------------------- Mixing 1-2-3-4 -------------------------------
    # ------------------------------------------------------------------------------     
    
    mixer_12345678          = ct.IdealGasReactor(Fuel_1)                            # [-]       Create the reactor for the downstream mixer
    outlet_1234             = ct.MassFlowController(mixer_1234, mixer_12345678, mdot = (mass_flow_rate_1 + mass_flow_rate_2 + mass_flow_rate_3 + mass_flow_rate_4)) # [-]       Connect the mixer with the downstream mixer     
    sim_mixer_1234          = ct.ReactorNet([mixer_1234])                           # [-]       Set the mixer simulation  
    sim_mixer_1234.advance_to_steady_state()                                        # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # ------------------------------- Mixing 5-6-7-8 -------------------------------
    # ------------------------------------------------------------------------------     

    outlet_5678             = ct.MassFlowController(mixer_5678, mixer_12345678, mdot = (mass_flow_rate_5 + mass_flow_rate_6 + mass_flow_rate_7 + mass_flow_rate_8))# [-]       Connect the mixer with the downstream mixer     
    sim_mixer_5678          = ct.ReactorNet([mixer_5678])                           # [-]       Set the mixer simulation  
    sim_mixer_5678.advance_to_steady_state()                                        # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # --------------------------- Mixing 1-2-3-4-5-6-7-8 ---------------------------
    # ------------------------------------------------------------------------------     

    mixer_air_1             = ct.IdealGasReactor(Fuel_1)                            # [-]       Create the reactor for the downstream mixer
    outlet_12345678         = ct.MassFlowController(mixer_12345678, mixer_air_1, mdot = (mass_flow_rate_1 + mass_flow_rate_2 + mass_flow_rate_3 + mass_flow_rate_4 + mass_flow_rate_5 + mass_flow_rate_6 + mass_flow_rate_7 + mass_flow_rate_8)) # [-]       Connect the mixer with the downstream mixer     
    sim_mixer_12345678      = ct.ReactorNet([mixer_12345678])                       # [-]       Set the mixer simulation  
    sim_mixer_12345678.advance_to_steady_state()                                    # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # -------------------------------- Mixing 1-Air --------------------------------
    # ------------------------------------------------------------------------------     

    Air_1                   = ct.Solution('Air.yaml')                               # [-]       Import air kinematic mechanism
    Air_1.TPX               = T_stag_0, P_stag_0, dict_oxy                          # [-]       Set the air temperature, pressure and mole fractions
    rho_air_1               = Air_1.density                                         # [kg/m**3] Fuel density
    res_air_1               = ct.Reservoir(Air_1)                                   # [-]       Create a resevoir for the air upstream of the mixer
    PFR_1                   = ct.IdealGasConstPressureReactor(Fuel_1)               # [-]       Create the reactor for the PFR
    PFR_1.volume            = A_SZ*(L_SZ/3)                                         # [m**3]    Set the PFR volume
    inlet_air_1             = ct.MassFlowController(res_air_1, mixer_air_1, mdot= m_dot_air_SZ) # [-]       Connect the upstream resevoir with the mixer
    m_dot_air_1             = mass_flow_rate_1 + mass_flow_rate_2 + mass_flow_rate_3 + mass_flow_rate_4 + mass_flow_rate_5 + mass_flow_rate_6 + mass_flow_rate_7 + mass_flow_rate_8 + m_dot_air_SZ                            # [kg/s]    Total mass flow inside the mixer  
    outlet_air_1            = ct.MassFlowController(mixer_air_1, PFR_1, mdot=m_dot_air_1) # [-]       Connect the mixer with the downstream PFR with the same mass flow rate           
    sim_mixer_air_1         = ct.ReactorNet([mixer_air_1])                          # [-]       Set the mixer simulation  
    sim_mixer_air_1.advance_to_steady_state()                                       # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PFR #1 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    mixer_air_2             = ct.IdealGasReactor(Fuel_1)                            # [-]       Create the reactor for the downstream mixer
    outlet_PFR_1            = ct.MassFlowController(PFR_1, mixer_air_2, mdot=m_dot_air_1) # [-]       Connect the PFR with the downstream mixer                                 
    t_res_PFR_1             = (Fuel_1.density * PFR_1.volume) / (m_dot_air_1)       # [s]       Compute the residence time inside the PFR           
    Fuel_1.set_equivalence_ratio(phi_SZ_des_1, fuel=dict_fuel, oxidizer=dict_oxy)   # [-]       Set the euivalence ratio inside the PFR
    sim_PFR_1               = ct.ReactorNet([PFR_1])                                # [-]       Set the PFR simulation
    sim_PFR_1.advance(t_res_PFR_1)                                                  # [-]       Run the simulation until the residence time is reached
    
    # ------------------------------------------------------------------------------
    # -------------------------------- Mixing 2-Air --------------------------------
    # ------------------------------------------------------------------------------     

    Air_2                   = ct.Solution('Air.yaml')                               # [-]       Import air kinematic mechanism
    Air_2.TPX               = T_stag_0, P_stag_0, dict_oxy                          # [-]       Set the air temperature, pressure and mole fractions
    rho_air_2               = Air_2.density                                         # [kg/m**3] Fuel density
    res_air_2               = ct.Reservoir(Air_2)                                   # [-]       Create a resevoir for the air upstream of the mixer
    PFR_2                   = ct.IdealGasConstPressureReactor(Fuel_1)               # [-]       Create the reactor for the PFR
    PFR_2.volume            = A_SZ*(L_SZ/3)                                         # [m**3]    Set the PFR volume
    inlet_air_2             = ct.MassFlowController(res_air_2, mixer_air_2, mdot=m_dot_air_SZ) # [-]       Connect the upstream resevoir with the mixer
    m_dot_air_2             = m_dot_air_1 + m_dot_air_SZ                            # [kg/s]    Total mass flow inside the mixer        
    outlet_air_2            = ct.MassFlowController(mixer_air_2, PFR_2, mdot=m_dot_air_2) # [-]       Connect the mixer with the downstream PFR with the same mass flow rate  
    sim_mixer_air_2         = ct.ReactorNet([mixer_air_2])                          # [-]       Set the mixer simulation  
    sim_mixer_air_2.advance_to_steady_state()                                       # [-]       Run the simulation until it reaches steady state
                                                                                    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PFR #2 -----------------------------------
    # ------------------------------------------------------------------------------ 
    
    mixer_air_3             = ct.IdealGasReactor(Fuel_1)                            # [-]       Create the reactor for the downstream mixer
    outlet_PFR_2            = ct.MassFlowController(PFR_2, mixer_air_3, mdot=m_dot_air_2) # [-]       Connect the PFR with the downstream mixer                                                               
    t_res_PFR_2             = (Fuel_1.density * PFR_2.volume) / (m_dot_air_2)       # [s]       Compute the residence time inside the PFR           
    Fuel_1.set_equivalence_ratio(phi_SZ_des_2, fuel=dict_fuel, oxidizer=dict_oxy)   # [-]       Set the euivalence ratio inside the PFR 
    sim_PFR_2               = ct.ReactorNet([PFR_2])                                # [-]       Set the PFR simulation
    sim_PFR_2.advance(t_res_PFR_2)                                                  # [-]       Run the simulation until the residence time is reached
    
    # ------------------------------------------------------------------------------
    # -------------------------------- Mixing 3-Air --------------------------------
    # ------------------------------------------------------------------------------     

    Air_3                   = ct.Solution('Air.yaml')                               # [-]       Import air kinematic mechanism
    Air_3.TPX               = T_stag_0, P_stag_0, dict_oxy                          # [-]       Set the air temperature, pressure and mole fractions
    rho_air_3               = Air_3.density                                         # [kg/m**3] Fuel density
    res_air_3               = ct.Reservoir(Air_3)                                   # [-]       Create a resevoir for the air upstream of the mixer
    PFR_3                   = ct.IdealGasConstPressureReactor(Fuel_1)               # [-]       Create the reactor for the PFR
    PFR_3.volume            = A_SZ*(L_SZ/3)                                         # [m**3]    Set the PFR volume
    inlet_air_3             = ct.MassFlowController(res_air_3, mixer_air_3, mdot=m_dot_air_SZ) # [-]       Connect the upstream resevoir with the mixer
    m_dot_air_3             = m_dot_air_2 + m_dot_air_SZ                            # [kg/s]    Total mass flow inside the mixer        
    outlet_air_3            = ct.MassFlowController(mixer_air_3, PFR_3, mdot=m_dot_air_3) # [-]       Connect the mixer with the downstream PFR with the same mass flow rate  
    sim_mixer_air_3         = ct.ReactorNet([mixer_air_3])                          # [-]       Set the mixer simulation  
    sim_mixer_air_3.advance_to_steady_state()                                       # [-]       Run the simulation until it reaches steady state
    
    # ------------------------------------------------------------------------------
    # ----------------------------------- PFR #3 -----------------------------------
    # ------------------------------------------------------------------------------ 
                                                           
    t_res_PFR_3             = (Fuel_1.density * PFR_3.volume) / (m_dot_air_3)       # [s]       Compute the residence time inside the PFR           
    Fuel_1.set_equivalence_ratio(phi_SZ_des_3, fuel=dict_fuel, oxidizer=dict_oxy)   # [-]       Set the euivalence ratio inside the PFR  
    sim_PFR_3               = ct.ReactorNet([PFR_3])                                # [-]       Set the PFR simulation                 
    sim_PFR_3.advance(t_res_PFR_3)                                                  # [-]       Run the simulation until the residence time is reached
                                                                                    
    # ------------------------------------------------------------------------------
    # ---------------------------- Additional computations -------------------------
    # ------------------------------------------------------------------------------    
          
    m_dot_input_combustor   = m_dot_fuel + m_dot_air_id                             # [kg/s]    Total mass flow rate entering a single combustor (air + fuel)
    Emission_Index          = Fuel_1.Y * (m_dot_input_combustor)/m_dot_fuel         # [-]       Computation of the Emission Index
    a_out                   = Fuel_1.sound_speed                                    # [m/s]     Speed of sound at PFR outlet
    rho_out                 = Fuel_1.density                                        # [kg/m**3] Density of the Fuel_1 
    gamma                   = Fuel_1.cp_mass / Fuel_1.cv_mass                       # [-]       Heat capacity ratio
    h                       = Fuel_1.h                                              # [J/kg]    Enthalpy
    vel_out                 = (m_dot_input_combustor) / (rho_out * A_SZ)            # [m/s]     Outlet velocity 
    M_out                   = vel_out / a_out                                       # [-]       Outlet Mach number   
    T_stag_out              = Fuel_1.T * (1 + 0.5 * (gamma - 1) * (M_out)**2)       # [K]       Stagnation temperature      
    P_stag_out              = Fuel_1.P * (1 + 0.5 * (gamma - 1) * (M_out)**2)**(gamma / (gamma - 1)) # [Pa]      Stagnation pressure  
    h_stag_out              = T_stag_out  * Fuel_1.cp_mass                          # [J/kg]    Stagnation enthalpy
    
    return (Fuel_1, Emission_Index, T_stag_out, P_stag_out, h_stag_out) 

if __name__ == '__main__': 
    main()
    plt.show()