# Created:  Jul 2023, M. Clarke
# Modified: Jun 2024, M. Guidotti & D.J. Lee

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 

# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core                                         import Units    
from RCAIDE.Library.Methods.Energy.Propulsors.Turboshaft_Propulsor import design_turboshaft, compute_turboshaft_performance
from RCAIDE.Library.Plots                                          import *     

# python imports 
import numpy                                                       as np  
import pickle
import matplotlib.pyplot                                           as plt    

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(): 

    #first_order_analysis()
    second_order_analysis()
    
    return 

def brayton_cycle_with_regenerator(): 
    
    # Inputs 
    M0                                         = 0                                                                                    #[-]
    T0                                         = 520                                                                                  #[R]
    gamma                                      = 1.4                                                                                  #[-]
    Cp                                         = 0.24                                                                                 #[Btu/(lbm*R)]
    LHV                                        = 18400                                                                                #[Btu/lbm]
    Tt4                                        = 2600                                                                                 #[R]
    pi_c                                       = np.linspace(2,18,100)                                                                #[-]
    x                                          = [1.02, 1.04]                                                                         #[-]
                                                                                                                                      
    # Equations                                                                                                                       
    tau_lambda                                 = Tt4/T0                                                                               #[-]                                                  
    tau_r                                      = 1 + ((gamma - 1)/2)*M0**2                                                            #[-]                                                                     
    tau_c                                      = pi_c**((gamma - 1)/gamma)                                                            #[-]
    C_shaft                                    = np.zeros((len(x), len(pi_c)))                                                        #[-]
    Psp                                        = np.zeros((len(x), len(pi_c)))                                                        #[hp/(lbm/sec)]   
    f                                          = np.zeros((len(x), len(pi_c)))                                                        #[-]  
    PSFC                                       = np.zeros((len(x), len(pi_c)))                                                        #[hp/(lb/hr)/hp]  
    eta_T                                      = np.zeros((len(x), len(pi_c)))                                                        #[-] 
    tau_t                                      = np.zeros((len(x), len(pi_c)))                                                        #[-]  
                                                                                                                                      
    for i in range(len(x)):                                                                                                           
        C_shaft[i,:]                           = tau_lambda*(1 - x[i]/(tau_r*tau_c)) - tau_r*(tau_c - 1)                              #[-]                                 
        Psp[i,:]                               = Cp*T0*C_shaft[i]*1.4148532041                                                        #[hp/(lbm/sec)]
        f[i,:]                                 = (Cp*T0*tau_lambda/LHV)*(1 - x[i]/(tau_r*tau_c))                                      #[-]
        PSFC[i,:]                              = f[i]/Psp[i]*3600/1.4148532041                                                        #[hp/(lb/hr)/hp] Units unclear
        eta_T[i,:]                             = 1 - (tau_r*(tau_c - 1))/(tau_lambda*(1 - x[i]/(tau_r*tau_c)))                        #[-]
        tau_t[i,:]                             = (x[i]/(tau_r*tau_c))                                                                 #[-]
        
    return Psp*np.ones(100), f*np.ones(100), PSFC*np.ones(100), eta_T*np.ones(100), C_shaft*np.ones(100), tau_t*np.ones(100), pi_c

def first_order_analysis():

    Psp, f, PSFC, eta_T, C_shaft, tau_t, pi_c  = brayton_cycle_with_regenerator()

    fig = plt.figure()
    axis = fig.add_subplot(1,1,1)
    axis.set_xlabel('Compressor Pressure Ratio [-]')
    axis.set_ylabel('Specific Power [hp/(lbm/sec)]')
    axis.plot(pi_c, Psp[0,:], label = 'x = 1.02' ) 
    axis.plot(pi_c, Psp[1,:], label = 'x = 1.04' ) 
    axis.legend()
    axis.grid()
    fig.tight_layout()
    
    fig = plt.figure()
    axis = fig.add_subplot(1,1,1)
    axis.set_xlabel('Compressor Pressure Ratio [-]')
    axis.set_ylabel('f [-]')
    axis.plot(pi_c, f[0,:], label = 'x = 1.02' ) 
    axis.plot(pi_c, f[1,:], label = 'x = 1.04' )
    axis.legend()
    axis.grid()
    fig.tight_layout()  
    
    fig = plt.figure()
    axis = fig.add_subplot(1,1,1)
    axis.set_xlabel('Compressor Pressure Ratio [-]')
    axis.set_ylabel('PSFC [hp/(lb/hr)/hp]')
    axis.plot(pi_c, PSFC[0,:], label = 'x = 1.02' ) 
    axis.plot(pi_c, PSFC[1,:], label = 'x = 1.04' )
    axis.legend()
    axis.grid()
    fig.tight_layout()     
    
    fig = plt.figure()
    axis = fig.add_subplot(1,1,1)
    axis.set_xlabel('Compressor Pressure Ratio [-]')
    axis.set_ylabel('eta_T [-]')
    axis.plot(pi_c, eta_T[0,:], label = 'x = 1.02' ) 
    axis.plot(pi_c, eta_T[1,:], label = 'x = 1.04' )
    axis.legend()
    axis.grid()
    fig.tight_layout()     
 
    fig = plt.figure()
    axis = fig.add_subplot(1,1,1)
    axis.set_xlabel('Compressor Pressure Ratio [-]')
    axis.set_ylabel('C_shaft [-]')
    axis.plot(pi_c, C_shaft[0,:], label = 'x = 1.02' )
    axis.plot(pi_c, C_shaft[1,:], label = 'x = 1.04' )
    axis.legend()
    axis.grid()
    fig.tight_layout()  
    
    fig = plt.figure()
    axis = fig.add_subplot(1,1,1)
    axis.set_xlabel('Compressor Pressure Ratio [-]')
    axis.set_ylabel('tau_t [-]')
    axis.plot(pi_c, tau_t[0,:], label = 'x = 1.02' )
    axis.plot(pi_c, tau_t[1,:], label = 'x = 1.04' )
    axis.legend()
    axis.grid()
    fig.tight_layout()   
    
    return 

def second_order_analysis():  

    altitude           = np.linspace(0,36000,10)*Units.feet 
    mach               = np.ones_like(altitude)*0.1 

    power_baseline , thermal_efficiency_baseline, PSFC = turboshaft_engine(altitude,mach) 

    fig = plt.figure('Turboshaft')
    axis = fig.add_subplot(1,3,1)
    axis.set_ylabel('altitude [ft]')
    axis.set_xlabel('Power [W]')
    axis.plot(power_baseline, altitude/Units.feet, color = 'black', label = 'baseline')    
    axis.grid()  
    axis.legend()    

    axis2 = fig.add_subplot(1,3,2)
    axis2.set_ylabel('altitude [ft]')
    axis2.set_xlabel('Thermal efficency [-]')
    axis2.plot(thermal_efficiency_baseline, altitude/Units.feet, color = 'black', label = 'baseline')  
    axis2.grid()
    axis2.legend()  
    
    axis3 = fig.add_subplot(1,3,3)
    axis3.set_ylabel('altitude [ft]')
    axis3.set_xlabel('PSFC [-]')
    axis3.plot(PSFC ,altitude/Units.feet, color = 'black', label = 'baseline')  
    axis3.grid()
    axis3.legend()     
    
    return 

def turboshaft_engine(altitude,mach):   

    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------         
    turboshaft                                     = RCAIDE.Library.Components.Propulsors.Turboshaft() 
    turboshaft.tag                                 = 'Turboshaft_propulsor'
    turboshaft.origin                              = [[13.72, 4.86,-1.1]] 
    turboshaft.engine_length                       = 2.71     
    turboshaft.bypass_ratio                        = 0    
    turboshaft.design_altitude                     = 36000.0*Units.ft
    turboshaft.design_mach_number                  = 0.1   
    turboshaft.design_power                        = 522000.0*Units.W 
    turboshaft.mass_flow_rate_design               = 1.9
        
                                                   
    # working fluid                                
    turboshaft.working_fluid                       = RCAIDE.Library.Attributes.Gases.Air() 
    ram                                            = RCAIDE.Library.Components.Propulsors.Converters.Ram()
    ram.tag                                        = 'ram' 
    turboshaft.ram                                 = ram 
                                                   
    # inlet nozzle                                 
    inlet_nozzle                                   = RCAIDE.Library.Components.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                               = 'inlet nozzle'
    inlet_nozzle.polytropic_efficiency             = 0.98
    inlet_nozzle.pressure_ratio                    = 0.98 
    turboshaft.inlet_nozzle                        = inlet_nozzle 
                                                   
    # compressor                                   
    compressor                                     = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    compressor.tag                                 = 'compressor'
    compressor.polytropic_efficiency               = 0.91
    compressor.pressure_ratio                      = 1.9  
    compressor.mass_flow_rate                      = 1.9 
    turboshaft.compressor                          = compressor

    # low pressure turbine  
    low_pressure_turbine                           = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    low_pressure_turbine.tag                       ='lpt'
    low_pressure_turbine.mechanical_efficiency     = 0.99
    low_pressure_turbine.polytropic_efficiency     = 0.93 
    turboshaft.low_pressure_turbine                = low_pressure_turbine
   
    # high pressure turbine     
    high_pressure_turbine                          = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    high_pressure_turbine.tag                      ='hpt'
    high_pressure_turbine.mechanical_efficiency    = 0.99
    high_pressure_turbine.polytropic_efficiency    = 0.93 
    turboshaft.high_pressure_turbine               = high_pressure_turbine 

    # combustor  
    combustor                                      = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                  = 'Comb'
    combustor.efficiency                           = 0.99 
    combustor.alphac                               = 1.0     
    combustor.turbine_inlet_temperature            = 1500
    combustor.pressure_ratio                       = 0.95
    combustor.fuel_data                            = RCAIDE.Library.Attributes.Propellants.Jet_A()  
    turboshaft.combustor                           = combustor

    # core nozzle
    core_nozzle                                    = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    core_nozzle.tag                                = 'core nozzle'
    core_nozzle.polytropic_efficiency              = 0.95
    core_nozzle.pressure_ratio                     = 0.99  
    turboshaft.core_nozzle                         = core_nozzle

    # design turboshaft
    design_turboshaft(turboshaft)    
     
    # connect turboshaft with network
    fuel_line                                      = RCAIDE.Library.Components.Energy.Distribution.Fuel_Line()    
    fuel_tank                                      = RCAIDE.Library.Components.Energy.Fuel_Tanks.Fuel_Tank()  
    fuel                                           = RCAIDE.Library.Attributes.Propellants.Aviation_Gasoline()    
    fuel_tank.fuel                                 = fuel  
    fuel_line.fuel_tanks.append(fuel_tank)  
    fuel_line.propulsors.append(turboshaft) 
    
    # ------------------------------------------------------------------------------------------------------------------------
    # Run engine 
    # ------------------------------------------------------------------------------------------------------------------------ 

    power                                                                         = np.zeros_like(altitude)
    thermal_efficiency                                                            = np.zeros_like(altitude)
    PSFC                                                                          = np.zeros_like(altitude)
                                                                                  
    for i in range(len(altitude)):                                                
                                                                                  
        # define atmospheric properties                                           
        atmosphere                                                                = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
        atmo_data                                                                 = atmosphere.compute_values(altitude[i],0.0)
        planet                                                                    = RCAIDE.Library.Attributes.Planets.Earth() 
        p                                                                         = atmo_data.pressure          
        T                                                                         = atmo_data.temperature       
        rho                                                                       = atmo_data.density          
        a                                                                         = atmo_data.speed_of_sound    
        mu                                                                        = atmo_data.dynamic_viscosity      
                                                                                  
        # initialize the operating conditions                                     
        state                                                                     = RCAIDE.Framework.Mission.Common.State() 
        state.conditions                                                          = RCAIDE.Framework.Mission.Common.Results()
                                                                                  
        # setting conditions for current simulation                               
        state.conditions.freestream.altitude                                      = np.atleast_2d(altitude[i])
        state.conditions.freestream.mach_number                                   = np.atleast_2d(mach[i])
        state.conditions.freestream.pressure                                      = np.atleast_2d(p)
        state.conditions.freestream.temperature                                   = np.atleast_2d(T)
        state.conditions.freestream.density                                       = np.atleast_2d(rho)
        state.conditions.freestream.dynamic_viscosity                             = np.atleast_2d(mu)
        state.conditions.freestream.gravity                                       = np.atleast_2d(planet.sea_level_gravity)
        state.conditions.freestream.isentropic_expansion_factor                   = np.atleast_2d(turboshaft.working_fluid.compute_gamma(T,p))
        state.conditions.freestream.Cp                                            = np.atleast_2d(turboshaft.working_fluid.compute_cp(T,p))
        state.conditions.freestream.R                                             = np.atleast_2d(turboshaft.working_fluid.gas_specific_constant)
        state.conditions.freestream.speed_of_sound                                = np.atleast_2d(a)
        state.conditions.freestream.velocity                                      = np.atleast_2d(a*mach[i])  
        
        # initialize data structure for turboshaft operating conditions (for energy) 
        state.conditions.energy[fuel_line.tag]                                    = RCAIDE.Framework.Mission.Common.Conditions()  
        state.conditions.energy[fuel_line.tag][fuel_tank.tag]                     = RCAIDE.Framework.Mission.Common.Conditions()  
        state.conditions.energy[fuel_line.tag][fuel_tank.tag].mass_flow_rate      = np.zeros((1,1))     
        state.conditions.energy[fuel_line.tag][fuel_tank.tag].mass                = np.zeros((1,1))   
        state.conditions.energy[fuel_line.tag][turboshaft.tag]                    = RCAIDE.Framework.Mission.Common.Conditions() 
        state.conditions.energy[fuel_line.tag][turboshaft.tag].throttle           = np.array([[1.0]])
        state.conditions.energy[fuel_line.tag][turboshaft.tag].thrust             = np.zeros((1,1))
        state.conditions.energy[fuel_line.tag][turboshaft.tag].power              = np.zeros((1,1)) 

        # initialize data structure for turshaft operating conditions (for noise )       
        state.conditions.noise[fuel_line.tag]                                     = RCAIDE.Framework.Mission.Common.Conditions()              
        state.conditions.noise[fuel_line.tag][turboshaft.tag]                     = RCAIDE.Framework.Mission.Common.Conditions() 
        state.conditions.noise[fuel_line.tag][turboshaft.tag].turboshaft          = RCAIDE.Framework.Mission.Common.Conditions() 
                
        total_power, eta_thermal, power_specific_fuel_consumption                 = compute_turboshaft_performance(fuel_line,state)
        
        power[i]                                                                  = total_power[0][0]
        thermal_efficiency[i]                                                     = eta_thermal[0][0]
        PSFC[i]                                                                   = power_specific_fuel_consumption[0][0]

    return power, thermal_efficiency, PSFC
    
# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def save_results(results,filename): 
    pickle_file  =  filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return   

# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results  

if __name__ == '__main__': 
    main()
    plt.show()    