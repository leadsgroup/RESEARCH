# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units   
from RCAIDE.Library.Methods.Propulsors.Turbofan_Propulsor          import design_turbofan

# python imports 
import numpy as np  
import pickle
from copy import deepcopy
import matplotlib.pyplot as plt  
import os   
import matplotlib.cm as cm

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(): 
    # Define Engine  
    turbofan  =  JT9D_7_turbofan_engine( )
  
    # Run engine
    altitude            = np.linspace(0,36000,1000) *Units.feet
    mach                = np.array([0.5])
    motor_work          = None
    regenerator         = None
    shaft_power_offtake = None
    thrust              = np.zeros_like(altitude )
    overall_efficiency  = np.zeros_like(altitude )
    
    for i in range(len(altitude)):

        atmosphere_sls = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
        atmo_data      = atmosphere_sls.compute_values(altitude[i],0.0)
        planet         = RCAIDE.Attributes.Planets.Earth()

        p   = atmo_data.pressure          
        T   = atmo_data.temperature       
        rho = atmo_data.density          
        a   = atmo_data.speed_of_sound    
        mu  = atmo_data.dynamic_viscosity      

        # setup conditions
        conditions_sls = RCAIDE.Analyses.Mission.Common.Results()            

        # freestream conditions    
        conditions_sls.freestream.altitude                    = np.atleast_2d(altitude[i])
        conditions_sls.freestream.mach_number                 = np.atleast_2d(mach[i])
        conditions_sls.freestream.pressure                    = np.atleast_2d(p)
        conditions_sls.freestream.temperature                 = np.atleast_2d(T)
        conditions_sls.freestream.density                     = np.atleast_2d(rho)
        conditions_sls.freestream.dynamic_viscosity           = np.atleast_2d(mu)
        conditions_sls.freestream.gravity                     = np.atleast_2d(planet.sea_level_gravity)
        conditions_sls.freestream.isentropic_expansion_factor = np.atleast_2d(turbofan.working_fluid.compute_gamma(T,p))
        conditions_sls.freestream.Cp                          = np.atleast_2d(turbofan.working_fluid.compute_cp(T,p))
        conditions_sls.freestream.R                           = np.atleast_2d(turbofan.working_fluid.gas_specific_constant)
        conditions_sls.freestream.speed_of_sound              = np.atleast_2d(a)
        conditions_sls.freestream.velocity                    = np.atleast_2d(a*mach[i])  

        net          = RCAIDE.Energy.Networks.Turbofan_Engine()  
        fuel_line    = RCAIDE.Energy.Distributors.Fuel_Line()  
        fuel_line.turbofans.append(turbofan)       
        net.fuel_lines.append(fuel_line)   

        conditions_sls.energy[fuel_line.tag]                                = RCAIDE.Analyses.Mission.Common.Conditions()         
        conditions_sls.noise[fuel_line.tag]                                 = RCAIDE.Analyses.Mission.Common.Conditions()      
        sorted_propulsors                                                   = compute_unique_propulsor_groups(fuel_line)
        fuel_line_results                                                   = conditions_sls.energy[fuel_line.tag] 
        fuel_line_results[turbofan.propulsor_group]                         = RCAIDE.Analyses.Mission.Common.Conditions()
        fuel_line_results[turbofan.propulsor_group].turbofan                = RCAIDE.Analyses.Mission.Common.Conditions() 
        fuel_line_results[turbofan.propulsor_group].unique_turbofan_tags    = sorted_propulsors.unique_turbofan_tags 
        fuel_line_results.N_turbofans                                       = sorted_propulsors.N_turbofans
        noise_results                                                       = conditions_sls.noise[fuel_line.tag]
        noise_results[turbofan.propulsor_group]                             = RCAIDE.Analyses.Mission.Common.Conditions() 
        noise_results[turbofan.propulsor_group].turbofan                    = RCAIDE.Analyses.Mission.Common.Conditions() 
        fuel_line_results[turbofan.propulsor_group].turbofan.throttle       = np.array([[1.0]])  

        T , P, mdot,sfc,eta_prop,eta = compute_propulsor_performance(0,fuel_line,turbofan.propulsor_group,fuel_line.turbofans,1,conditions_sls)   

        thrust[i]             = np.linalg.norm(T)
        overall_efficiency[i] = eta  
 
    return 
    

def JT9D_7_turbofan_engine(): 

    turbofan                                    = RCAIDE.Library.Components.Propulsors.Turbofan() 
    turbofan.tag                                = 'turbofan'
    turbofan.origin                             = [[13.72, 4.86,-1.1]] 
    turbofan.engine_length                      = 2.71     
    turbofan.bypass_ratio                       = 4.8 # checked  
    turbofan.design_altitude                    = 35000.0*Units.ft # checked  
    turbofan.design_mach_number                 = 0.78    # checked  
    turbofan.design_thrust                      = 10000*Units.lbf # checked   

    # fan                
    fan                                         = RCAIDE.Library.Components.Propulsors.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.56 
    turbofan.fan                                = fan        

    # working fluid                   
    turbofan.working_fluid                      = RCAIDE.Library.Attributes.Gases.Air()
    ram                                         = RCAIDE.Library.Components.Propulsors.Converters.Ram()
    ram.tag                                     = 'ram' 
    turbofan.ram                                = ram 

    # inlet nozzle          
    inlet_nozzle                                = RCAIDE.Library.Components.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                            = 'inlet nozzle'
    inlet_nozzle.polytropic_efficiency          = 0.98
    inlet_nozzle.pressure_ratio                 = 0.98 
    turbofan.inlet_nozzle                       = inlet_nozzle 

    # low pressure compressor    
    low_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    low_pressure_compressor.tag                   = 'lpc'
    low_pressure_compressor.polytropic_efficiency = 0.91
    low_pressure_compressor.pressure_ratio        = 2.3 
    turbofan.low_pressure_compressor              = low_pressure_compressor

    # high pressure compressor  
    high_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    high_pressure_compressor.tag                   = 'hpc'
    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 10.2 

    # combustor  
    combustor                                      = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                  = 'Comb'
    combustor.efficiency                           = 0.99 
    combustor.alphac                               = 1.0     
    combustor.turbine_inlet_temperature            = 1273.889
    combustor.pressure_ratio                       = 0.95
    combustor.fuel_data                            = RCAIDE.Library.Attributes.Propellants.Jet_A1()  
    turbofan.combustor                             = combustor

    # low pressure turbine  
    low_pressure_turbine                           = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    low_pressure_turbine.tag                       ='lpt'
    low_pressure_turbine.mechanical_efficiency     = 0.99
    low_pressure_turbine.polytropic_efficiency     = 0.93 
    turbofan.low_pressure_turbine                  = low_pressure_turbine

    # high pressure turbine     
    high_pressure_turbine                          = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    high_pressure_turbine.tag                      ='hpt'
    high_pressure_turbine.mechanical_efficiency    = 0.99
    high_pressure_turbine.polytropic_efficiency    = 0.93 
    turbofan.high_pressure_turbine                 = high_pressure_turbine 

    # core nozzle
    core_nozzle                                    = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    core_nozzle.tag                                = 'core nozzle'
    core_nozzle.polytropic_efficiency              = 0.95
    core_nozzle.pressure_ratio                     = 0.99  
    turbofan.core_nozzle                           = core_nozzle

    # fan nozzle             
    fan_nozzle                                     = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    fan_nozzle.tag                                 = 'fan nozzle'
    fan_nozzle.polytropic_efficiency               = 0.95
    fan_nozzle.pressure_ratio                      = 0.99 
    turbofan.fan_nozzle                            = fan_nozzle 

    # design turbofan
    design_turbofan(turbofan)
    
    
    return turbofan 

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