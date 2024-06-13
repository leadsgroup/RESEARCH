# Created:  Jul 2023, M. Clarke
# Modified: Jun 2024, M. Guidotti & D.J. Lee

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 

# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core                                         import Units    
#from RCAIDE.Library.Methods.Energy.Propulsors.Turboshaft_Propulsor import design_turboshaft, compute_turboshaft_performance
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
    #second_order_analysis()
    
    return 

def brayton_cycle_with_regenerator(): 
    
    # Inputs 
    M0                                         = 0                                  #[-]
    T0                                         = 520                                #[R]
    gamma                                      = 1.4                                #[-]
    Cp                                         = 0.24                               #[Btu/(lbm*R)]
    LHV                                        = 18400                              #[Btu/lbm]
    Tt4                                        = 2600                               #[R]
    pi_c                                       = np.linspace(2,18,100)              #[-]
    x                                          = [1.02, 1.04]                       #[-]
    
    # Equations   
    tau_lambda                                 = Tt4/T0                             #[-]                                                  
    tau_r                                      = 1 + ((gamma - 1)/2)*M0**2          #[-]                                                                     
    tau_c                                      = pi_c**((gamma - 1)/gamma)          #[-]
    C_shaft                                    = np.zeros((len(x), len(pi_c)))      #[-]
    Psp                                        = np.zeros((len(x), len(pi_c)))      #[hp/(lbm/sec)]   
    f                                          = np.zeros((len(x), len(pi_c)))      #[-]  
    PSFC                                       = np.zeros((len(x), len(pi_c)))      #[hp/(lb/hr)/hp]  
    eta_T                                      = np.zeros((len(x), len(pi_c)))      #[-] 
    tau_t                                      = np.zeros((len(x), len(pi_c)))      #[-]  
 
    for i in range(len(x)): 
        C_shaft[i,:]                           = tau_lambda*(1 - x[i]/(tau_r*tau_c)) - tau_r*(tau_c - 1)       #[-]                                 
        Psp[i,:]                               = Cp*T0*C_shaft[i]*1.4148532041                                 #[hp/(lbm/sec)]
        f[i,:]                                 = (Cp*T0*tau_lambda/LHV)*(1 - x[i]/(tau_r*tau_c))               #[-]
        PSFC[i,:]                              = f[i]/Psp[i]*3600/1.4148532041                                 #[hp/(lb/hr)/hp] Units unclear
        eta_T[i,:]                             = 1 - (tau_r*(tau_c - 1))/(tau_lambda*(1 - x[i]/(tau_r*tau_c))) #[-]
        tau_t[i,:]                             = (x[i]/(tau_r*tau_c))                                          #[-]
        
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
    mach               = np.ones_like(altitude)*0.001 # np.linspace(0.1,0.7,100)

    motor_work = None 
    thrust_baseline , overall_efficiency_baseline  = JT9D_7_turbofan_engine(altitude,mach,motor_work) 

    fig = plt.figure('Turbofan')
    axis = fig.add_subplot(1,2,1)
    axis.set_ylabel('altitude')
    axis.set_xlabel('Thrust (lbf)')
    axis.plot(thrust_baseline/Units.lbf,altitude/Units.feet, color = 'black', label = 'baseline')    
    axis.grid()   

    axis2 = fig.add_subplot(1,2,2)
    axis2.set_ylabel('altitude')
    axis2.set_xlabel(r'overall efficency')
    axis2.plot(overall_efficiency_baseline*100 ,altitude/Units.feet, color = 'black', label = 'baseline')  
    axis2.grid()  


    markers = ['s','x','P','v','o','^','D','*','1']
    motor_work  = np.linspace(0,0.99,5)
    colors  = cm.inferno(np.linspace(0.2,0.8,len(motor_work)))  
    for i in range(len(motor_work)):
        thrust_motor , overall_efficiency_motor  = JT9D_7_turbofan_engine(altitude,mach,motor_work[i])    
        motor_tag = str(round(motor_work[i]*100,0)) +  '% $W_m$' 
        axis.plot(thrust_motor/Units.lbf,altitude/Units.feet,linestyle = '--', color =  colors[i], label = motor_tag  )     
        axis2.plot(overall_efficiency_motor*100,altitude/Units.feet,linestyle = '--', color =  colors[i] )    


        regenerator = RCAIDE.Energy.Thermal_Management.Turbine.Regenerator() 
        regenerator.combustor.temperature_ratio = 0.975
        regenerator.nozzle.temperature_ratio    = 0.975 
        motor_regen_tag =str(round(motor_work[i]*100,0)) +  '% $W_m$' +  r', $\tau_x$ = ' + str(regenerator.combustor.temperature_ratio)
        thrust_motor_regen , overall_efficiency_motor_regen  = JT9D_7_turbofan_engine(altitude,mach,motor_work[i],regenerator)   
        axis.plot(thrust_motor_regen/Units.lbf,altitude/Units.feet,linestyle = '-',  color =  colors[i], marker = markers[i], label = motor_regen_tag  )     
        axis2.plot(overall_efficiency_motor_regen*100,altitude/Units.feet, linestyle = '-',color =  colors[i], marker = markers[i]  )  

        #shaft_power_offtake = RCAIDE.Energy.Converters.Shaft_Power_Off_Take()
        #shaft_power_offtake.power_draw = 1000000 
        #regenerator = RCAIDE.Energy.Thermal_Management.Turbine.Regenerator() 
        #regenerator.combustor.temperature_ratio = 1.0
        #regenerator.nozzle.temperature_ratio    = 1.0 
        #motor_regen_tag ='$W_m$ = ' + str(motor_work[i]) +  r', $\tau_x$ = ' + str(regenerator.combustor.temperature_ratio) + ' 1MW Generator'
        #thrust_motor_regen_gen , overall_efficiency_motor_regen_gen  = JT9D_7_turbofan_engine(altitude,mach,motor_work[i],regenerator,shaft_power_offtake)   
        #axis3.plot(thrust_motor_regen_gen/Units.lbf,altitude/Units.feet,linestyle = '-',  color =  colors[i], marker = markers[i], label = motor_regen_tag  )     
        #axis4.plot(overall_efficiency_motor_regen_gen*100,altitude/Units.feet, linestyle = '-',color =  colors[i], marker = markers[i], label = motor_regen_tag )  


    axis.legend()           
    #axis2.legend()               
    fig.tight_layout()      
    return 

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














#def turboshaft_engine_setup(): 
    
    ## ------------------------------------------------------------------
    ##   Initialize the Vehicle
    ## ------------------------------------------------------------------    
    
    #vehicle = RCAIDE.Vehicle()
    #vehicle.tag = ''        

    ## ################################################# Energy Network #######################################################         
    ## Step 1: Define network
    ## Step 2: Define Distribution Type
    ## Step 3: Define Propulsors 
    ## Step 4: Define Enegy Source 

    ##------------------------------------------------------------------------------------------------------------------------- 
    ##  Turbofan Network
    ##-------------------------------------------------------------------------------------------------------------------------   
    #net                                         = RCAIDE.Framework.Networks.Turbofan_Engine_Network() 
    
    ##------------------------------------------------------------------------------------------------------------------------- 
    ## Fuel Distrubition Line 
    ##------------------------------------------------------------------------------------------------------------------------- 
    #fuel_line                                   = RCAIDE.Library.Components.Energy.Distribution.Fuel_Line()  
    
    ##------------------------------------------------------------------------------------------------------------------------------------  
    ## Propulsor: Starboard Propulsor
    ##------------------------------------------------------------------------------------------------------------------------------------         
    #turbofan                                    = RCAIDE.Library.Components.Propulsors.Turbofan() 
    #turbofan.tag                                = 'starboard_propulsor'
    #turbofan.active_fuel_tanks                  = ['fuel_tank']   
    #turbofan.origin                             = [[13.72, 4.86,-1.1]] 
    #turbofan.engine_length                      = 2.71     
    #turbofan.bypass_ratio                       = 5.4    
    #turbofan.design_altitude                    = 35000.0*Units.ft
    #turbofan.design_mach_number                 = 0.78   
    #turbofan.design_thrust                      = 35000.0*Units.N 
             
    ## fan                
    #fan                                         = RCAIDE.Library.Components.Propulsors.Converters.Fan()   
    #fan.tag                                     = 'fan'
    #fan.polytropic_efficiency                   = 0.93
    #fan.pressure_ratio                          = 1.7   
    #turbofan.fan                                = fan        
                   
    ## working fluid                   
    #turbofan.working_fluid                      = RCAIDE.Library.Attributes.Gases.Air() 
    #ram                                         = RCAIDE.Library.Components.Propulsors.Converters.Ram()
    #ram.tag                                     = 'ram' 
    #turbofan.ram                                = ram 
          
    ## inlet nozzle          
    #inlet_nozzle                                = RCAIDE.Library.Components.Propulsors.Converters.Compression_Nozzle()
    #inlet_nozzle.tag                            = 'inlet nozzle'
    #inlet_nozzle.polytropic_efficiency          = 0.98
    #inlet_nozzle.pressure_ratio                 = 0.98 
    #turbofan.inlet_nozzle                       = inlet_nozzle 

    ## low pressure compressor    
    #low_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    #low_pressure_compressor.tag                   = 'lpc'
    #low_pressure_compressor.polytropic_efficiency = 0.91
    #low_pressure_compressor.pressure_ratio        = 1.9   
    #turbofan.low_pressure_compressor              = low_pressure_compressor

    ## high pressure compressor  
    #high_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    #high_pressure_compressor.tag                   = 'hpc'
    #high_pressure_compressor.polytropic_efficiency = 0.91
    #high_pressure_compressor.pressure_ratio        = 10.0    
    #turbofan.high_pressure_compressor              = high_pressure_compressor

    ## low pressure turbine  
    #low_pressure_turbine                           = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    #low_pressure_turbine.tag                       ='lpt'
    #low_pressure_turbine.mechanical_efficiency     = 0.99
    #low_pressure_turbine.polytropic_efficiency     = 0.93 
    #turbofan.low_pressure_turbine                  = low_pressure_turbine
   
    ## high pressure turbine     
    #high_pressure_turbine                          = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    #high_pressure_turbine.tag                      ='hpt'
    #high_pressure_turbine.mechanical_efficiency    = 0.99
    #high_pressure_turbine.polytropic_efficiency    = 0.93 
    #turbofan.high_pressure_turbine                 = high_pressure_turbine 

    ## combustor  
    #combustor                                      = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    #combustor.tag                                  = 'Comb'
    #combustor.efficiency                           = 0.99 
    #combustor.alphac                               = 1.0     
    #combustor.turbine_inlet_temperature            = 1500
    #combustor.pressure_ratio                       = 0.95
    #combustor.fuel_data                            = RCAIDE.Library.Attributes.Propellants.Jet_A()  
    #turbofan.combustor                             = combustor

    ## core nozzle
    #core_nozzle                                    = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    #core_nozzle.tag                                = 'core nozzle'
    #core_nozzle.polytropic_efficiency              = 0.95
    #core_nozzle.pressure_ratio                     = 0.99  
    #turbofan.core_nozzle                           = core_nozzle
             
    ## fan nozzle             
    #fan_nozzle                                     = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    #fan_nozzle.tag                                 = 'fan nozzle'
    #fan_nozzle.polytropic_efficiency               = 0.95
    #fan_nozzle.pressure_ratio                      = 0.99 
    #turbofan.fan_nozzle                            = fan_nozzle 
    
    ## design turbofan
    #design_turbofan(turbofan)  
    ## append propulsor to distribution line  
    
    #fuel_line.propulsors.append(turbofan)  

    ##------------------------------------------------------------------------------------------------------------------------- 
    ##  Energy Source: Fuel Tank
    ##------------------------------------------------------------------------------------------------------------------------- 
    ## fuel tank
    #fuel_tank                                   = RCAIDE.Library.Components.Energy.Fuel_Tanks.Fuel_Tank()
    
    ## append fuel 
    #fuel                                        = RCAIDE.Library.Attributes.Propellants.Aviation_Gasoline()   
    #fuel.mass_properties.mass                   = vehicle.mass_properties.max_takeoff-vehicle.mass_properties.max_fuel
    #fuel.origin                                 = vehicle.wings.main_wing.mass_properties.center_of_gravity      
    #fuel.mass_properties.center_of_gravity      = vehicle.wings.main_wing.aerodynamic_center
    #fuel.internal_volume                        = fuel.mass_properties.mass/fuel.density  
    #fuel_tank.fuel                              = fuel            
    
    ## apend fuel tank to dataclass of fuel tanks on fuel line 
    #fuel_line.fuel_tanks.append(fuel_tank) 

    ## Append fuel line to Network      
    #net.fuel_lines.append(fuel_line)   

    ## Append energy network to aircraft 
    #vehicle.append_energy_network(net)    
      
    #return vehicle
