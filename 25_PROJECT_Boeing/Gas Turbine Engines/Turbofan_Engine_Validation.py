# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units , Data   
from RCAIDE.Library.Methods.Propulsors.Turbofan_Propulsor          import design_turbofan
from RCAIDE.Framework.Mission.Common      import  Conditions
from RCAIDE.Library.Methods.Emissions.emissions_index_CRN_method import  emissions_index_CRN_method

# python imports 
import numpy as np  
import pickle
from copy import deepcopy
import matplotlib.pyplot as plt  
import os   
import matplotlib.cm as cm
import time 

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():  
    # Run engine
    altitude            = np.array([35000])  *Units.feet # np.linspace(0,35000,8) *Units.feet
    mach_number         = np.array([0.78])  # np.linspace(1E-4,0.8,9)
    thrust              = np.zeros((len(altitude),len(mach_number)))
    overall_efficiency  = np.zeros((len(altitude),len(mach_number)))
    thermal_efficiency  = np.zeros((len(altitude),len(mach_number)))
    
    Tt_3           = np.zeros((len(altitude),len(mach_number)))
    Pt_3           = np.zeros((len(altitude),len(mach_number)))
    Tt_4           = np.zeros((len(altitude),len(mach_number)))
    Pt_4           = np.zeros((len(altitude),len(mach_number))) 
    m_dot_core     = np.zeros((len(altitude),len(mach_number)))
    fuel_flow_rate = np.zeros((len(altitude),len(mach_number)))
    
    turbofan       = JT9D_7_turbofan_engine()            
    for i in range(len(altitude)): 
        for j in range(len(mach_number)):
            planet         = RCAIDE.Library.Attributes.Planets.Earth()
            atmosphere     = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
            atmo_data      = atmosphere.compute_values(altitude[i])
            
            p   = atmo_data.pressure          
            T   = atmo_data.temperature       
            rho = atmo_data.density          
            a   = atmo_data.speed_of_sound    
            mu  = atmo_data.dynamic_viscosity     
                
            conditions = RCAIDE.Framework.Mission.Common.Results() 
            conditions.freestream.altitude                    = np.atleast_1d(0)
            conditions.freestream.mach_number                 = np.atleast_1d(mach_number[j])
            conditions.freestream.pressure                    = np.atleast_1d(p)
            conditions.freestream.temperature                 = np.atleast_1d(T)
            conditions.freestream.density                     = np.atleast_1d(rho)
            conditions.freestream.dynamic_viscosity           = np.atleast_1d(mu)
            conditions.freestream.gravity                     = np.atleast_2d(planet.sea_level_gravity)
            conditions.freestream.isentropic_expansion_factor = np.atleast_1d(turbofan.working_fluid.compute_gamma(T,p))
            conditions.freestream.Cp                          = np.atleast_1d(turbofan.working_fluid.compute_cp(T,p))
            conditions.freestream.R                           = np.atleast_1d(turbofan.working_fluid.gas_specific_constant)
            conditions.freestream.speed_of_sound              = np.atleast_1d(a)
            conditions.freestream.velocity                    = np.atleast_1d(a*mach_number[j])  
        
            ## setup conditions  
            fuel_line                = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line()
            segment                  = RCAIDE.Framework.Mission.Segments.Segment()  
            segment.state.conditions = conditions     
            segment.state.conditions.energy[fuel_line.tag] = Conditions()
            segment.state.conditions.noise[fuel_line.tag]  = Conditions()
            turbofan.append_operating_conditions(segment,fuel_line) 
            for tag, item in  turbofan.items(): 
                if issubclass(type(item), RCAIDE.Library.Components.Component):
                    item.append_operating_conditions(segment,fuel_line,turbofan) 
            
            # set throttle
            segment.state.conditions.energy[fuel_line.tag][turbofan.tag].throttle[:,0] = 1.0  
            Thrust,_,_,_,_ = turbofan.compute_performance(segment.state,fuel_line)
                  
            ram                       = turbofan.ram
            inlet_nozzle              = turbofan.inlet_nozzle
            low_pressure_compressor   = turbofan.low_pressure_compressor
            high_pressure_compressor  = turbofan.high_pressure_compressor
            fan                       = turbofan.fan
            combustor                 = turbofan.combustor
            high_pressure_turbine     = turbofan.high_pressure_turbine
            low_pressure_turbine      = turbofan.low_pressure_turbine
            core_nozzle               = turbofan.core_nozzle
            fan_nozzle                = turbofan.fan_nozzle 
            bypass_ratio              = turbofan.bypass_ratio  
        
            # unpack component conditions
            turbofan_conditions     = conditions.energy[fuel_line.tag][turbofan.tag]
            ram_conditions          = turbofan_conditions[ram.tag]    
            fan_conditions          = turbofan_conditions[fan.tag]    
            inlet_nozzle_conditions = turbofan_conditions[inlet_nozzle.tag]
            core_nozzle_conditions  = turbofan_conditions[core_nozzle.tag]
            fan_nozzle_conditions   = turbofan_conditions[fan_nozzle.tag]
            lpc_conditions          = turbofan_conditions[low_pressure_compressor.tag]
            hpc_conditions          = turbofan_conditions[high_pressure_compressor.tag]
            lpt_conditions          = turbofan_conditions[low_pressure_turbine.tag]
            hpt_conditions          = turbofan_conditions[high_pressure_turbine.tag]
            combustor_conditions    = turbofan_conditions[combustor.tag]
            
            # compute emission indices 
            #emissions_index_CRN_method(combustor,turbofan_conditions,conditions)
            
            # extract properties
            U_e             = core_nozzle_conditions.outputs.velocity 
            U_e1            = fan_nozzle_conditions.outputs.velocity  
            mdot_air_core   = turbofan_conditions.core_mass_flow_rate
            mdot_air_fan    = bypass_ratio *  mdot_air_core  
            fuel_enthalpy   = combustor.fuel_data.specific_energy 
            mdot_fuel       = turbofan_conditions.fuel_flow_rate 
            U_0             = a*mach_number[j] 
            h_e1            = fan_nozzle_conditions.outputs.static_enthalpy
            h_e             = core_nozzle_conditions.outputs.static_enthalpy
            h_0             = turbofan.working_fluid.compute_cp(T,p) * T 
            h_t4            = combustor_conditions.outputs.stagnation_enthalpy
            h_t3            = hpc_conditions.outputs.stagnation_enthalpy     
            
            thrust[i,j]             = np.linalg.norm(Thrust)
            
            # Overall Efficiency;  Aircraft and Rocket Propulsion Eqn 2.22 
            overall_efficiency[i,j] = thrust[i,j] * U_0 / (mdot_fuel * fuel_enthalpy)
            
            # Thermal efficiecny ;  Aircraft and Rocket Propulsion Eqn 5.49  
            thermal_efficiency[i,j] = 1 -  (  (mdot_air_core +  mdot_fuel)*(h_e -  h_0) + mdot_air_fan*(h_e1 - h_0) + mdot_fuel *h_0  ) /( (mdot_air_core +  mdot_fuel)*h_t4 - mdot_air_core *h_t3 )
            
            Tt_3[i,j]           = hpc_conditions.outputs.stagnation_temperature 
            Pt_3[i,j]           = hpc_conditions.outputs.stagnation_pressure
            Tt_4[i,j]           = hpt_conditions.inputs.stagnation_temperature 
            Pt_4[i,j]           = hpt_conditions.inputs.stagnation_pressure 
            m_dot_core[i,j]     = turbofan_conditions.core_mass_flow_rate   
            fuel_flow_rate[i,j] = turbofan_conditions.fuel_flow_rate                    

                
    plot_results(altitude,mach_number,thrust,overall_efficiency,thermal_efficiency,Tt_3,Pt_3,Tt_4,Pt_4,m_dot_core,fuel_flow_rate)
    
    return

def plot_results(altitude,mach_number,thrust,overall_efficiency,thermal_efficiency,Tt_3,Pt_3,Tt_4,Pt_4,m_dot_core,fuel_flow_rate):
    ps =  plot_style(number_of_lines = len(mach_number)) 
    
    fig    =  plt.figure('Thrust')
    fig.set_size_inches(7, 6)
    axis_1 = fig.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_1.plot(thrust[:,i]/Units.lbf,altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_1.set_xlabel('Thrust (lbf)')
    axis_1.set_ylabel('Altitude (ft)')
    axis_1.legend()
    fig.tight_layout()
    

    fig_2    =  plt.figure('Thermal Efficiency')
    fig_2.set_size_inches(7, 6)
    axis_2 = fig_2.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_2.plot(thermal_efficiency[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0],linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_2.set_xlabel('Thermal Efficiency')
    axis_2.set_ylabel('Altitude (ft)')
    axis_2.legend()
    fig_2.tight_layout()
    

    fig_3    =  plt.figure('Overall Efficiency')
    fig_3.set_size_inches(7, 6)
    axis_3 = fig_3.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_3.plot(overall_efficiency[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_3.set_xlabel('Overall Efficiency')
    axis_3.set_ylabel('Altitude (ft)')
    axis_3.legend()
    fig_3.tight_layout()
    
    


    fig_4    =  plt.figure('Stagnation Properties')
    fig_4.set_size_inches(7, 6)
    axis_4_1 = fig_4.add_subplot(2,2,1)
    axis_4_2 = fig_4.add_subplot(2,2,2)
    axis_4_3 = fig_4.add_subplot(2,2,3)
    axis_4_4 = fig_4.add_subplot(2,2,4) 
    for i in  range(len(mach_number)):
        axis_4_1.plot(Tt_3[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0], marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2)))
        axis_4_2.plot(Pt_3[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0], marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2)))
        axis_4_3.plot(Tt_4[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0], marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2)))
        axis_4_4.plot(Pt_4[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0], marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_4_1.set_xlabel(r'Entering $T_t$ (K)')
    axis_4_2.set_xlabel(r'Entering $P_t$ (K)')
    axis_4_3.set_xlabel(r'Exiting $T_t$ (Pa)')
    axis_4_4.set_xlabel(r'Exising $P_t$ (Pa)')
    axis_4_1.set_ylabel('Altitude (ft)')
    axis_4_2.set_ylabel('Altitude (ft)')
    axis_4_3.set_ylabel('Altitude (ft)')
    axis_4_4.set_ylabel('Altitude (ft)')
    axis_4_1.legend()
    axis_4_2.legend()
    axis_4_3.legend()
    axis_4_4.legend()
    fig_4.tight_layout()    
    

    fig_5    =  plt.figure('Core Mass Flow')
    fig_5.set_size_inches(7, 6)
    axis_5 = fig_5.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_5.plot(m_dot_core[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_5.set_xlabel(r'$dot{m}_{core}$ (kg/s)')
    axis_5.set_ylabel('Altitude (ft)')
    axis_5.legend()
    fig_5.tight_layout()    
 

    fig_6    =  plt.figure('Fuel Flow Rate')
    fig_6.set_size_inches(7, 6)
    axis_6 = fig_6.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_6.plot(fuel_flow_rate[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_6.set_xlabel(r'$\dot{m}_{core}$ (kg/s)')
    axis_6.set_ylabel('Altitude (ft)')
    axis_6.legend()
    fig_6.tight_layout()
    
    return

def plot_style(number_of_lines= 10): 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18,
                  #figure.dpi': 1200
                  }


    # Universal Plot Settings  
    plt.rcParams.update(parameters)
    plot_parameters                        = Data()
    plot_parameters.line_width             = 1.5  
    plot_parameters.line_style             = ['-','--']
    plot_parameters.marker_size            = 4
    plot_parameters.legend_fontsize        = '12'
    plot_parameters.legend_title_font_size = 14
    plot_parameters.axis_font_size         = 16
    plot_parameters.title_font_size        = 16   
    plot_parameters.markers                =  ['o','x','o','v','P','p','^','D','*']
    plot_parameters.color                  = cm.inferno(np.linspace(0,0.9,number_of_lines)) 


    return plot_parameters

def GE_90_engine(PSR_PFR_combustor_model_flag ):

    
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propulsor: Starboard Propulsor
    #------------------------------------------------------------------------------------------------------------------------------------         
    turbofan                                    = RCAIDE.Library.Components.Propulsors.Turbofan() 
    turbofan.tag                                = 'starboard_ge90_propulsor'
    turbofan.active_fuel_tanks                  = ['b777_fuel_tank']   
    turbofan.origin                             = [[ 25.72797886 , 9.69802 , -2.04  ]]
    turbofan.mass_properties.mass               = 7893
    turbofan.engine_length                      = 7.29
    turbofan.bypass_ratio                       = 9
    turbofan.design_altitude                    = 35000.0*Units.ft
    turbofan.design_mach_number                 = 0.78   
    turbofan.design_thrust                      = 80000  * Units.N  


    # fan                
    fan                                         = RCAIDE.Library.Components.Propulsors.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.7   
    turbofan.fan                                = fan        

    # working fluid                   
    turbofan.working_fluid                      = RCAIDE.Library.Attributes.Gases.Air() 

    
    # Ram inlet 
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
    low_pressure_compressor.pressure_ratio        = 1.9   
    turbofan.low_pressure_compressor              = low_pressure_compressor

    ## high pressure compressor  
    #medium_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    #medium_pressure_compressor.tag                   = 'hpc'
    #medium_pressure_compressor.polytropic_efficiency = 0.91
    #medium_pressure_compressor.pressure_ratio        = 12.38 
    #turbofan.high_pressure_compressor              = high_pressure_compressor

    # high pressure compressor  
    high_pressure_compressor                       = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    high_pressure_compressor.tag                   = 'hpc'
    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 12.38 
    turbofan.high_pressure_compressor              = high_pressure_compressor
    

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

    # combustor  
    combustor                                      = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                  = 'Comb'
    combustor.efficiency                           = 0.99 
    combustor.alphac                               = 1.0     
    combustor.turbine_inlet_temperature            = 1500
    combustor.pressure_ratio                       = 0.95
    combustor.fuel_data                            = RCAIDE.Library.Attributes.Propellants.Jet_A()  
    turbofan.combustor                             = combustor

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
    # append propulsor to distribution line 
    
    return turbofan

def JT9D_7_turbofan_engine(): 

    turbofan                                    = RCAIDE.Library.Components.Propulsors.Turbofan() 
    turbofan.tag                                = 'turbofan'
    turbofan.origin                             = [[13.72, 4.86,-1.1]] 
    turbofan.engine_length                      = 2.71     
    turbofan.engine_diameter                    = 2.37
    turbofan.bypass_ratio                       = 4.8 # checked  
    turbofan.design_altitude                    = 35000.0*Units.ft # checked  
    turbofan.design_mach_number                 = 0.78    # checked
    turbofan.design_thrust                      = 10000*Units.lbf # checked   

    # fan                
    fan                                         = RCAIDE.Library.Components.Propulsors.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.6 
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
    turbofan.high_pressure_compressor              = high_pressure_compressor

    # combustor  
    combustor                                      = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                  = 'Comb'
    combustor.efficiency                           = 0.99 
    combustor.alphac                               = 1.0     
    combustor.turbine_inlet_temperature            = 1273.889
    combustor.pressure_ratio                       = 0.95 
    combustor.fuel_data                            = RCAIDE.Library.Attributes.Propellants.Jet_A()  
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