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
    P_stag_0                = 1500000                                   # [Pa] 
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
        Fuel.TP = max(T_stag_0, 300), P_stag_0
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
            
        #air_mass_flow = max(0, beta_air_in * air.density * dz)
        #mfc_air.mass_flow_rate = air_mass_flow
        #outlet = ct.MassFlowController(PFR_1, mixer_2)
        #outlet.mass_flow_rate = mdot_air + mdot_fuel + air_mass_flow
        
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
        outlet.mass_flow_rate = mdot_air + mdot_fuel + air_mass_flow
        
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

if __name__ == '__main__': 
    main()
    plt.show()