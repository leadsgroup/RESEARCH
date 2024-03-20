## @ingroup Visualization-Performance-Energy-Thermal_Management
# RCAIDE/Visualization/Performance/Energy/Thermal_Management/plot_heat_acquisition_system_conditions.py
# 
# 
# Created:  Jul 2023, M. Clarke

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  

from RCAIDE.Core import Units, Data 
from RCAIDE.Visualization.Common import set_axes, plot_style
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np 
import pickle 

def main():     

    # vehicle data  
    file_name_1 = 'X_57_HEX'   
    res_Hex     = load_results(file_name_1) 
    
    file_name_2 = 'X_57_w_no_HEX'     
    res_No_Hex = load_results(file_name_2)
    
    plot_heat_acquisition_system_conditions(res_Hex,res_No_Hex)
    plot_battery_cell_conditions(res_Hex,res_No_Hex)
    return
# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS
# ----------------------------------------------------------------------------------------------------------------------    
## @ingroup Visualization-Performance-Common
def set_axes(axes):
    """This sets the axis parameters for all plots

    Assumptions:
    None

    Source:
    None

    Inputs
    axes

    Outputs:
    axes

    Properties Used:
    N/A 
    """

    axes.minorticks_on()
    #axes.grid(which='major', linestyle='-', linewidth=0.5, color='grey')
    #axes.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')
    axes.grid(True)
    #axes.get_yaxis().get_major_formatter().set_scientific(False)
    #axes.get_yaxis().get_major_formatter().set_useOffset(False)

    return

 
# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Visualization-Performance-Common
def plot_style(): 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 12,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18}
    plt.rcParams.update(parameters)

    # Universal Plot Settings  
    plot_parameters                  = Data()
    plot_parameters.line_width       = 2 
    plot_parameters.line_style       = ['--','-']
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 12
    plot_parameters.axis_font_size   = 14
    plot_parameters.title_font_size  = 18   
    plot_parameters.markers          = [None,None] # ['s','X','o','v','P','p','^','D','*']
    plot_parameters.color            = 'black'
    
    return plot_parameters

# ----------------------------------------------------------------------------------------------------------------------
#   plot_heat_exchanger_system_conditions
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Visualization-Performance-Energy-Thermal_Management
def plot_heat_acquisition_system_conditions(res_Hex,res_No_Hex,
                                  save_figure = False,
                                  show_legend = True,
                                  save_filename = "Heat_Acquisition_System",
                                  file_type = ".png"): 
    
    # get plotting style 
    ps      = plot_style()    
     
    # get line colors for plots 
    line_colors   = cm.inferno(np.linspace(0,0.9,len(res_Hex.segments)))     

    fig = plt.figure(save_filename)
    fig.set_size_inches(5,4) 
    axis_0 = plt.subplot(1,1,1)
    axis_1 = plt.subplot(3,1,1)
    axis_2 = plt.subplot(3,1,2) 
    axis_3 = plt.subplot(3,1,3)          
    b_i = 0 
    for network in res_Hex.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                #axis_0.plot(np.zeros(2),np.zeros(2), color = line_colors[0], marker = ps.markers[0], linewidth = ps.line_width,label= battery.tag) 
                #axis_0.plot(np.zeros(2),np.zeros(2), color = line_colors[0], marker = ps.markers[0], linewidth = ps.line_width,label= battery.tag) 
                axis_0.grid(False)
                axis_0.axis('off')  
               
                for i in range(len(res_Hex.segments)):  
                    time                       = res_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
                    battery_conditions_HEX     = res_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    
                    battery_conditions_No_HEX  = res_No_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    
                    
    
                    outlet_coolant_temperature_HEX    = battery_conditions_HEX.thermal_management_system.HAS.outlet_coolant_temperature[:,0]
                    coolant_mass_flow_rate_HEX        = battery_conditions_HEX.thermal_management_system.HAS.coolant_mass_flow_rate[:,0]
                    power_HEX                         = battery_conditions_HEX.thermal_management_system.HAS.power[:,0]     
                    outlet_coolant_temperature_No_HEX = battery_conditions_No_HEX.thermal_management_system.HAS.outlet_coolant_temperature[:,0]
                    coolant_mass_flow_rate_No_HEX     = battery_conditions_No_HEX.thermal_management_system.HAS.coolant_mass_flow_rate[:,0]
                    power_No_HEX                      = battery_conditions_No_HEX.thermal_management_system.HAS.power[:,0]                 
      
            
                    segment_tag  = res_Hex.segments[i].tag
                    segment_name = segment_tag.replace('_', ' ') 

                    if b_i == 0:                     
                        axis_1.plot(time, outlet_coolant_temperature_HEX, linestyle = ps.line_style[0], color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width)
                    else:
                        axis_1.plot(time, outlet_coolant_temperature_HEX, linestyle = ps.line_style[0], color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width)
                    #axis_1.plot(time, outlet_coolant_temperature_No_HEX, linestyle = ps.line_style[0], color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_1.set_ylabel(r'T$_c$ (K)')  
                    set_axes(axis_1)     
                     
                    axis_2.plot(time, coolant_mass_flow_rate_HEX, linestyle = ps.line_style[0],  color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width)
                    #axis_2.plot(time, coolant_mass_flow_rate_No_HEX, linestyle = ps.line_style[0],  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_2.set_ylabel(r'$\dot{m}_c$ (kg/s)') 
                    set_axes(axis_2)     
             
                    axis_3.plot(time, power_HEX, color = line_colors[i], linestyle = ps.line_style[0],  marker = ps.markers[0], linewidth = ps.line_width)
                    #axis_3.plot(time, power_No_HEX/1000, color = line_colors[i], linestyle = ps.line_style[0],  marker = ps.markers[1], linewidth = ps.line_width)
                    axis_3.set_ylabel(r'Power$_{HAS}$ (W)')
                    axis_3.set_xlabel(r'Time (mins)') 
                    set_axes(axis_3)     
                                  
                b_i += 1 
    fig.tight_layout()
    if save_figure:
        plt.savefig(save_filename + battery.tag + file_type)    
    return 

# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Visualization-Performance-Energy-Battery
def plot_battery_cell_conditions(res_Hex,res_No_Hex,
                                  save_figure = False,
                                  show_legend = True,
                                  save_filename = "Battery_Cell_Conditions_",
                                  file_type = ".png",
                                  width = 12, height = 7): 
    
    # get plotting style 
    ps      = plot_style()   
     
    # get line colors for plots 
    line_colors   = cm.inferno(np.linspace(0,0.9,len(res_Hex.segments)))     

    fig = plt.figure(save_filename)
    fig.set_size_inches(5,4)  
    axis_1 = plt.subplot(2,2,1) 
    axis_3 = plt.subplot(2,2,2) 
    axis_4 = plt.subplot(2,2,3) 
    axis_6 = plt.subplot(2,2,4)      
    b_i = 0 
    for network in res_Hex.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:    
                
                for i in range(len(res_Hex.segments)):  
                    time    = res_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
                    battery_conditions_Hex  = res_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    
                    cell_power_Hex          = battery_conditions_Hex.cell.power[:,0] 
                    cell_current_Hex        = battery_conditions_Hex.cell.current[:,0]
                    cell_SOC_Hex            = battery_conditions_Hex.cell.state_of_charge[:,0]   
                    cell_temperature_Hex    = battery_conditions_Hex.cell.temperature[:,0]  
                    

                    battery_conditions_No_Hex  = res_No_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    
                    cell_power_No_Hex          = battery_conditions_No_Hex.cell.power[:,0] 
                    cell_current_No_Hex        = battery_conditions_No_Hex.cell.current[:,0]
                    cell_SOC_No_Hex            = battery_conditions_No_Hex.cell.state_of_charge[:,0]   
                    cell_temperature_No_Hex    = battery_conditions_No_Hex.cell.temperature[:,0]                      
             

                    if b_i == 0:                     
                        axis_1.plot(time, cell_SOC_Hex, linestyle = ps.line_style[0],  color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width)
                    else:
                        axis_1.plot(time, cell_SOC_Hex, linestyle = ps.line_style[0],  color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width)
                    axis_1.plot(time, cell_SOC_No_Hex, linestyle = ps.line_style[1],  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_1.set_ylabel(r'SOC')
                    axis_1.set_ylim([0,1.1])
                    set_axes(axis_1)     
                      
             
                    axis_3.plot(time, cell_current_Hex, linestyle = ps.line_style[0],  color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width)
                    axis_3.plot(time, cell_current_No_Hex, linestyle = ps.line_style[1],  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_3.set_ylabel(r'Current (A)')
                    set_axes(axis_3)  
             
                    axis_4.plot(time, cell_power_Hex, linestyle = ps.line_style[0],  color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width)
                    axis_4.plot(time, cell_power_No_Hex, linestyle = ps.line_style[1],  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_4.set_ylabel(r'Power (W)')
                    set_axes(axis_4)      
             
                    axis_6.plot(time, cell_temperature_Hex, linestyle = ps.line_style[0],  color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width)
                    axis_6.plot(time, cell_temperature_No_Hex, linestyle = ps.line_style[1],  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_6.set_ylabel(r'Temperature, $\degree$C')
                    axis_6.set_xlabel(r'Time (mins)')
                    set_axes(axis_6)    
    
    fig.tight_layout()
    if save_figure:
        plt.savefig(save_filename + battery.tag + file_type)    
    return fig 

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
