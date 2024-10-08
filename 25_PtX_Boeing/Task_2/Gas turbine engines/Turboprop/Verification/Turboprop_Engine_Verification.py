# RCAIDE imports 
import RCAIDE
from   RCAIDE.Framework.Core                                 import Units , Data   
from   RCAIDE.Library.Methods.Propulsors.Turboprop_Propulsor import design_turboprop
from   RCAIDE.Framework.Mission.Common                       import  Conditions

# Python imports 
import numpy                                                as np                                             
import matplotlib.pyplot                                    as plt 
import matplotlib.cm                                        as cm
import pickle    
import os                                                   
import time 

'''This code is created to validate the performance of a Turboprop engine in RCAIDE.'''

# Created:  Sep 2024, M. Clarke, M. Guidotti 

# References:
# [1]: 

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():  
    
    ti                      = time.time()                               # [s]       Define the initial simulation time

    altitude                = np.array([35000])*Units.feet 
    mach_number             = np.array([0.8]) 
    
    thrust                  = np.zeros((len(altitude),len(mach_number)))
    propulsive_efficiency   = np.zeros((len(altitude),len(mach_number)))
    thermal_efficiency      = np.zeros((len(altitude),len(mach_number)))
    Tt_3                    = np.zeros((len(altitude),len(mach_number)))
    Pt_3                    = np.zeros((len(altitude),len(mach_number)))
    Tt_4                    = np.zeros((len(altitude),len(mach_number)))
    Pt_4                    = np.zeros((len(altitude),len(mach_number))) 
    fuel_flow_rate          = np.zeros((len(altitude),len(mach_number)))
    m_dot_air_tot           = np.zeros((len(altitude),len(mach_number)))
    TSFC                    = np.zeros((len(altitude),len(mach_number)))
                            
    Turboprop               = Turboprop_engine()
    
    for i in range(len(altitude)): 
        for j in range(len(mach_number)):
            
            planet                                            = RCAIDE.Library.Attributes.Planets.Earth()
            atmosphere                                        = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
            atmo_data                                         = atmosphere.compute_values(altitude[i])
                                                              
            p                                                 = atmo_data.pressure          
            T                                                 = atmo_data.temperature       
            rho                                               = atmo_data.density          
            a                                                 = atmo_data.speed_of_sound    
            mu                                                = atmo_data.dynamic_viscosity     
                
            conditions                                        = RCAIDE.Framework.Mission.Common.Results() 
            conditions.freestream.altitude                    = np.atleast_1d(altitude[i])
            conditions.freestream.mach_number                 = np.atleast_1d(mach_number[j])
            conditions.freestream.pressure                    = np.atleast_1d(p)
            conditions.freestream.temperature                 = np.atleast_1d(T)
            conditions.freestream.density                     = np.atleast_1d(rho)
            conditions.freestream.dynamic_viscosity           = np.atleast_1d(mu)
            conditions.freestream.gravity                     = np.atleast_2d(planet.sea_level_gravity)
            conditions.freestream.isentropic_expansion_factor = np.atleast_1d(Turboprop.working_fluid.compute_gamma(T,p))
            conditions.freestream.Cp                          = np.atleast_1d(Turboprop.working_fluid.compute_cp(T,p))
            conditions.freestream.R                           = np.atleast_1d(Turboprop.working_fluid.gas_specific_constant)
            conditions.freestream.speed_of_sound              = np.atleast_1d(a)
            conditions.freestream.velocity                    = np.atleast_1d(a*mach_number[j])  
        
            # setup conditions  
            fuel_line                                         = RCAIDE.Library.Components.Energy.Distributors.Fuel_Line()
            segment                                           = RCAIDE.Framework.Mission.Segments.Segment()  
            segment.state.conditions                          = conditions     
            segment.state.conditions.energy[fuel_line.tag]    = Conditions()
            segment.state.conditions.noise[fuel_line.tag]     = Conditions()
            
            Turboprop.append_operating_conditions(segment,fuel_line) 
            
            for tag, item in  Turboprop.items(): 
                if issubclass(type(item), RCAIDE.Library.Components.Component):
                    item.append_operating_conditions(segment,fuel_line,Turboprop) 
            
            # set throttle
            segment.state.conditions.energy[fuel_line.tag][Turboprop.tag].throttle[:,0] = 1  
            
            Thrust,_,_,_,_                                    = Turboprop.compute_performance(segment.state,fuel_line)
                  
            ram                                               = Turboprop.ram
            inlet_nozzle                                      = Turboprop.inlet_nozzle
            compressor                                        = Turboprop.compressor
            combustor                                         = Turboprop.combustor
            high_pressure_turbine                             = Turboprop.high_pressure_turbine
            low_pressure_turbine                              = Turboprop.low_pressure_turbine
            core_nozzle                                       = Turboprop.core_nozzle
        
            # unpack component conditions
            Turboprop_conditions                              = conditions.energy[fuel_line.tag][Turboprop.tag]
            ram_conditions                                    = Turboprop_conditions[ram.tag]    
            inlet_nozzle_conditions                           = Turboprop_conditions[inlet_nozzle.tag]
            compressor_conditions                             = Turboprop_conditions[compressor.tag]
            combustor_conditions                              = Turboprop_conditions[combustor.tag]
            hpt_conditions                                    = Turboprop_conditions[high_pressure_turbine.tag]
            lpt_conditions                                    = Turboprop_conditions[low_pressure_turbine.tag]
            core_nozzle_conditions                            = Turboprop_conditions[core_nozzle.tag]          
            
            # extract properties
            U_e_c                                             = core_nozzle_conditions.outputs.velocity 
            mdot_air_core                                     = Turboprop_conditions.core_mass_flow_rate
            fuel_enthalpy                                     = combustor.fuel_data.specific_energy 
            mdot_fuel                                         = Turboprop_conditions.fuel_flow_rate 
            U_0                                               = a*mach_number[j] 
            h_e_c                                             = core_nozzle_conditions.outputs.static_enthalpy
            h_0                                               = Turboprop.working_fluid.compute_cp(T,p) * T 
            h_t4                                              = combustor_conditions.outputs.stagnation_enthalpy
            h_t3                                              = compressor_conditions.outputs.stagnation_enthalpy     
            
            thrust[i,j]                                       = np.linalg.norm(Thrust)
            propulsive_efficiency[i,j]                        = Turboprop_conditions.eta_P
            thermal_efficiency[i,j]                           = Turboprop_conditions.eta_T
            Tt_3[i,j]                                         = compressor_conditions.outputs.stagnation_temperature 
            Pt_3[i,j]                                         = compressor_conditions.outputs.stagnation_pressure
            Tt_4[i,j]                                         = hpt_conditions.inputs.stagnation_temperature 
            Pt_4[i,j]                                         = hpt_conditions.inputs.stagnation_pressure  
            fuel_flow_rate[i,j]                               = Turboprop_conditions.fuel_flow_rate
            m_dot_air_tot[i,j]                                = Turboprop_conditions.core_mass_flow_rate 
            TSFC[i,j]                                         = Turboprop.TSFC # [N/N-s]
            
            
    print(thrust, propulsive_efficiency, thermal_efficiency, Tt_3, Pt_3, Tt_4, Pt_4, fuel_flow_rate, m_dot_air_tot, TSFC)                  
    #plot_results(altitude,mach_number,thrust,overall_efficiency,thermal_efficiency,Tt_3,Pt_3,Tt_4,Pt_4,m_dot_core,fuel_flow_rate,m_dot_air_tot)
    
    tf                      = time.time()                                           # [s]       Define the final simulation time
    elapsed_time            = round((tf-ti),2)                                      # [s]       Compute the total simulation time

    print('Simulation Time: ' + str(elapsed_time) + ' seconds per timestep')        # [-]       Print the value of total simulation time    
    
    return

def plot_results(altitude,mach_number,thrust,overall_efficiency,thermal_efficiency,Tt_3,Pt_3,Tt_4,Pt_4,m_dot_core,fuel_flow_rate,m_dot_air_tot):
    ps =  plot_style(number_of_lines = len(mach_number)) 
    
    fig    =  plt.figure('Thrust')
    fig.set_size_inches(7, 6)
    axis_1 = fig.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_1.plot(thrust[:,i],altitude, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2)))     
    axis_1.set_xlabel('Thrust [N]')
    axis_1.set_ylabel('Altitude [m]')
    
    axis_1.legend()
    fig.tight_layout()
    
    fig_2    =  plt.figure('Thermal Efficiency')
    fig_2.set_size_inches(7, 6)
    axis_2 = fig_2.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_2.plot(thermal_efficiency[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0],linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_2.set_xlabel('Thermal Efficiency')
    axis_2.set_ylabel('Altitude [ft]')
    axis_2.legend()
    fig_2.tight_layout()
    
    fig_3    =  plt.figure('Overall Efficiency')
    fig_3.set_size_inches(7, 6)
    axis_3 = fig_3.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_3.plot(overall_efficiency[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_3.set_xlabel('Overall Efficiency')
    axis_3.set_ylabel('Altitude [ft]')
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
    axis_4_1.set_xlabel(r'Entering $T_t$ [K]')
    axis_4_2.set_xlabel(r'Entering $P_t$ [K]')
    axis_4_3.set_xlabel(r'Exiting $T_t$ [Pa]')
    axis_4_4.set_xlabel(r'Exiting $P_t$ [Pa]')
    axis_4_1.set_ylabel('Altitude [ft]')
    axis_4_2.set_ylabel('Altitude [ft]')
    axis_4_3.set_ylabel('Altitude [ft]')
    axis_4_4.set_ylabel('Altitude [ft]')
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
    axis_5.set_xlabel(r'$\dot{m}_{core}$ [kg/s]')
    axis_5.set_ylabel('Altitude [ft]')
    axis_5.legend()
    fig_5.tight_layout()    
 

    fig_6    =  plt.figure('Fuel Flow Rate')
    fig_6.set_size_inches(7, 6)
    axis_6 = fig_6.add_subplot(1,1,1) 
    for i in  range(len(mach_number)):
        axis_6.plot(fuel_flow_rate[:,i],altitude/Units.feet, color = ps.color[i], linestyle = ps.line_style[0],
                    marker = ps.markers[0], linewidth = ps.line_width, label = 'Mach =' + str( round(mach_number[i], 2))) 
    axis_6.set_xlabel(r'Flow Rate [kg/s]')
    axis_6.set_ylabel('Altitude [ft]')
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

def Turboprop_engine():
        
    Turboprop                                    = RCAIDE.Library.Components.Propulsors.Turboprop()    
    Turboprop.origin                             = [[ 0.0 , 0.0 , 0.0 ]]
    Turboprop.design_altitude                    = 35006.562*Units.ft       
    Turboprop.design_mach_number                 = 0.8                     
    Turboprop.design_thrust                      = 77850 * Units.N       
        
    # working fluid                   
    Turboprop.working_fluid                      = RCAIDE.Library.Attributes.Gases.Air() 
    
    # Propeller and Gearbox efficiency               
    Turboprop.design_propeller_efficiency        = 0.83    
    Turboprop.design_gearbox_efficiency          = 0.99  
    
    # Ram inlet 
    ram                                          = RCAIDE.Library.Components.Propulsors.Converters.Ram()
    ram.tag                                      = 'ram' 
    Turboprop.ram                                = ram 
          
    # inlet nozzle          
    inlet_nozzle                                 = RCAIDE.Library.Components.Propulsors.Converters.Compression_Nozzle()
    inlet_nozzle.tag                             = 'inlet nozzle'                                       
    inlet_nozzle.pressure_ratio                  = 0.98
    inlet_nozzle.compressibility_effects         = False
    Turboprop.inlet_nozzle                       = inlet_nozzle
                                                 
    # compressor                    
    compressor                                   = RCAIDE.Library.Components.Propulsors.Converters.Compressor()    
    compressor.tag                               = 'lpc'                   
    compressor.pressure_ratio                    = 10                   
    Turboprop.compressor                         = compressor

    # combustor  
    combustor                                    = RCAIDE.Library.Components.Propulsors.Converters.Combustor()   
    combustor.tag                                = 'Comb'
    combustor.efficiency                         = 0.99                   
    combustor.turbine_inlet_temperature          = 1370                    
    combustor.pressure_ratio                     = 0.96                    
    combustor.fuel_data                          = RCAIDE.Library.Attributes.Propellants.Jet_A()  
    Turboprop.combustor                          = combustor
    
    # high pressure turbine     
    high_pressure_turbine                        = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    high_pressure_turbine.tag                    ='hpt'
    high_pressure_turbine.mechanical_efficiency  = 0.99                       
    Turboprop.high_pressure_turbine              = high_pressure_turbine 
    
    # low pressure turbine  
    low_pressure_turbine                         = RCAIDE.Library.Components.Propulsors.Converters.Turbine()   
    low_pressure_turbine.tag                     ='lpt'
    low_pressure_turbine.mechanical_efficiency   = 0.99                      
    Turboprop.low_pressure_turbine               = low_pressure_turbine

    # core nozzle
    core_nozzle                                  = RCAIDE.Library.Components.Propulsors.Converters.Expansion_Nozzle()   
    core_nozzle.tag                              = 'core nozzle'          
    core_nozzle.pressure_ratio                   = 0.99
    Turboprop.core_nozzle                        = core_nozzle
    
    # design Turboprop
    design_turboprop(Turboprop)  
    
    return Turboprop

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