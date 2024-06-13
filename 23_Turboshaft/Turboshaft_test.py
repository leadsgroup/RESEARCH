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
from copy                                                          import deepcopy
import matplotlib.pyplot                                           as plt  
import os   
import matplotlib.cm                                               as cm

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(): 

    #comparison_of_perfomrance()
    first_order_analysis()
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
    mach               = np.ones_like(altitude)*0.001 

    power_baseline , overall_efficiency_baseline = turboshaft_engine(altitude,mach) 

    fig = plt.figure('Turboshaft')
    axis = fig.add_subplot(1,2,1)
    axis.set_ylabel('altitude [ft]')
    axis.set_xlabel('Power [W]')
    axis.plot(power_baseline/Units.Watts,altitude/Units.feet, color = 'black', label = 'baseline')    
    axis.grid()  
    axis.legend()    

    axis2 = fig.add_subplot(1,2,2)
    axis2.set_ylabel('altitude [ft]')
    axis2.set_xlabel('Overall efficency [-]')
    axis2.plot(overall_efficiency_baseline*100 ,altitude/Units.feet, color = 'black', label = 'baseline')  
    axis2.grid()
    axis2.legend()  
            
    fig.tight_layout()      
    
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
    turboshaft.design_altitude                     = 35000.0*Units.ft
    turboshaft.design_mach_number                  = 0.1  
    turboshaft.design_power                        = 522000.0*Units.W 
                                                   
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
    
    # ------------------------------------------------------------------------------------------------------------------------
    # Run engine 
    # ------------------------------------------------------------------------------------------------------------------------ 

    #thrust             = np.zeros_like(altitude )
    #overall_efficiency = np.zeros_like(altitude )
    #for i in range(len(altitude)):

        #atmosphere_sls = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
        #atmo_data      = atmosphere_sls.compute_values(altitude[i],0.0)
        #planet         = RCAIDE.Attributes.Planets.Earth()

        #p   = atmo_data.pressure          
        #T   = atmo_data.temperature       
        #rho = atmo_data.density          
        #a   = atmo_data.speed_of_sound    
        #mu  = atmo_data.dynamic_viscosity      

        ## setup conditions
        #conditions_sls = RCAIDE.Analyses.Mission.Common.Results()            

        ## freestream conditions    
        #conditions_sls.freestream.altitude                    = np.atleast_2d(altitude[i])
        #conditions_sls.freestream.mach_number                 = np.atleast_2d(mach[i])
        #conditions_sls.freestream.pressure                    = np.atleast_2d(p)
        #conditions_sls.freestream.temperature                 = np.atleast_2d(T)
        #conditions_sls.freestream.density                     = np.atleast_2d(rho)
        #conditions_sls.freestream.dynamic_viscosity           = np.atleast_2d(mu)
        #conditions_sls.freestream.gravity                     = np.atleast_2d(planet.sea_level_gravity)
        #conditions_sls.freestream.isentropic_expansion_factor = np.atleast_2d(turboshaft.working_fluid.compute_gamma(T,p))
        #conditions_sls.freestream.Cp                          = np.atleast_2d(turboshaft.working_fluid.compute_cp(T,p))
        #conditions_sls.freestream.R                           = np.atleast_2d(turboshaft.working_fluid.gas_specific_constant)
        #conditions_sls.freestream.speed_of_sound              = np.atleast_2d(a)
        #conditions_sls.freestream.velocity                    = np.atleast_2d(a*mach[i])  

        #net          = RCAIDE.Energy.Networks.turboshaft_Engine()  
        #fuel_line    = RCAIDE.Energy.Distributors.Fuel_Line()  
        #fuel_line.turboshafts.append(turboshaft)       
        #net.fuel_lines.append(fuel_line)   

        #conditions_sls.energy[fuel_line.tag]                                = RCAIDE.Analyses.Mission.Common.Conditions()         
        #conditions_sls.noise[fuel_line.tag]                                 = RCAIDE.Analyses.Mission.Common.Conditions()      
        #sorted_propulsors                                                   = compute_unique_propulsor_groups(fuel_line)
        #fuel_line_results                                                   = conditions_sls.energy[fuel_line.tag] 
        #fuel_line_results[turboshaft.propulsor_group]                         = RCAIDE.Analyses.Mission.Common.Conditions()
        #fuel_line_results[turboshaft.propulsor_group].turboshaft                = RCAIDE.Analyses.Mission.Common.Conditions() 
        #fuel_line_results[turboshaft.propulsor_group].unique_turboshaft_tags    = sorted_propulsors.unique_turboshaft_tags 
        #fuel_line_results.N_turboshafts                                       = sorted_propulsors.N_turboshafts
        #noise_results                                                       = conditions_sls.noise[fuel_line.tag]
        #noise_results[turboshaft.propulsor_group]                             = RCAIDE.Analyses.Mission.Common.Conditions() 
        #noise_results[turboshaft.propulsor_group].turboshaft                    = RCAIDE.Analyses.Mission.Common.Conditions() 
        #fuel_line_results[turboshaft.propulsor_group].turboshaft.throttle       = np.array([[1.0]])  

        #T , P, mdot,sfc,eta_prop,eta = compute_turboshaft_performance(0,fuel_line,turboshaft.propulsor_group,fuel_line.turboshafts,1,conditions_sls)   

        #thrust[i]             = np.linalg.norm(T)
        #overall_efficiency[i] = eta 

    #return thrust , overall_efficiency
    
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