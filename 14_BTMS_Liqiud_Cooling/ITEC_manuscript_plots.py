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
import pandas as pd 
import pickle 

def main():     

    # vehicle data  
    file_name_1 = 'X_57_HEX'   
    res_Hex     = load_results(file_name_1) 

    file_name_2 = 'X_57_w_no_HEX'     
    res_No_Hex = load_results(file_name_2)

    #plot_flight_profile(res_No_Hex)
    #plot_reservoir_conditions(res_Hex)
    #plot_heat_acquisition_system_conditions(res_Hex)
    #plot_heat_exchanger_system_conditions(res_Hex)
    #plot_battery_cell_conditions_1(res_Hex,res_No_Hex)
    #plot_battery_cell_conditions_2(res_Hex,res_No_Hex)
    #plot_cooling_drag(res_Hex,res_No_Hex)
    #plot_percentage_operation(res_Hex)
    save_current_csv(res_Hex)

    return
# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS COMMON AXIS 
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
#  PLOTS STYLE
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Visualization-Performance-Common
def plot_style(): 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18,
                  'figure.dpi': 600
                  }


    # Universal Plot Settings  
    plt.rcParams.update(parameters)
    plot_parameters                        = Data()
    plot_parameters.line_width             = 1  
    plot_parameters.line_style             = ['-','-']
    plot_parameters.marker_size            = 4
    plot_parameters.legend_fontsize        = '12'
    plot_parameters.legend_title_font_size = 14
    plot_parameters.axis_font_size         = 16
    plot_parameters.title_font_size        = 16   
    plot_parameters.markers                =  ['o','x','o','v','P','p','^','D','*']
    plot_parameters.color                  = 'black'


    return plot_parameters

def save_current_csv(res_No_Hex,
                     save_csv = True,
                        show_legend = True,
                        save_filename = "Ambient_temperature",
                        file_type = ".csv"):
    current = []
    t    = []
    for network in res_No_Hex.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:         
    
                for i in range(len(res_No_Hex.segments)): 
        
                    time     = res_No_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.s
                    battery_conditions_No_Hex  = res_No_Hex.segments[i].conditions.energy[bus.tag][battery.tag]         
                    pack_current        = battery_conditions_No_Hex.cell.current[:,0]
                    reservoir_temperature = battery_conditions_No_Hex.thermal_management_system.RES.coolant_temperature[:,0]
                    t_ambient           = res_No_Hex.segments[i].conditions.freestream.temperature[:,0] 
                    
                    current.append(t_ambient)
                    t.append(time)
          
        current = np.hstack(current)
        t       = np.hstack(t)
        df = pd.DataFrame({'Time': t, 'Current': current})
        
        if save_csv:
            df.to_csv(f'IEEE_ITEC_final_plots/{save_filename}{file_type}', index=False)
        return 
    

# ----------------------------------------------------------------------------------------------------------------------
#   Plot Flight Profile and Airspeed
# ----------------------------------------------------------------------------------------------------------------------  
def plot_flight_profile(res_No_Hex,
                        save_figure = True,
                        show_legend = True,
                        save_filename = "Flight_Profile_aircraft_speed",
                        file_type = ".png"):
    # get plotting style 
    ps      = plot_style()    

    # get line colors for plots 
    line_colors   = cm.inferno(np.linspace(0,0.9,len(res_No_Hex.segments)))     
    fig = plt.figure(save_filename)
    fig.set_size_inches(5,6)  
    axis_0 = plt.subplot(1,1,1)
    axis_1 = plt.subplot(3,1,1)
    axis_2 = plt.subplot(3,1,2) 

    axis_0.grid(False)
    axis_0.axis('off')         
    for i in range(len(res_No_Hex.segments)): 

        time     = res_No_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min
        airspeed = res_No_Hex.segments[i].conditions.freestream.velocity[:,0] /   Units['mph']
        altitude = res_No_Hex.segments[i].conditions.freestream.altitude[:,0]/Units.feet

        segment_tag  =  res_No_Hex.segments[i].tag
        segment_name = segment_tag.replace('_', ' ')


        axis_1.plot(time, altitude, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width, label = segment_name)
        axis_1.set_ylabel(r'Altitude (ft)', fontsize = ps.axis_font_size)
        axis_1.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
        set_axes(axis_1)    


        axis_2.plot(time, airspeed, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width )
        axis_2.set_ylabel(r'Airspeed (mph)',fontsize = ps.axis_font_size)
        axis_2.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
        set_axes(axis_2) 


    if show_legend:        
        leg =  fig.legend(bbox_to_anchor=(0.515, 0.05), loc='lower center', ncol = 3, fontsize = ps.legend_fontsize) 
        leg.set_title('Flight Segment', prop={'size': ps.legend_title_font_size, 'weight': 'heavy'})    

    # Adjusting the sub-plots for legend 
    #fig.subplots_adjust(top=0.8)

    # set title of plot 
    #title_text    = 'Flight Profile for eCTOL'      
    #fig.suptitle(title_text, fontsize=ps.title_font_size, fontweight='bold')    
    fig.tight_layout()

    if save_figure:
        plt.savefig(f'IEEE_ITEC_final_plots/{save_filename}{file_type}') 

    return fig


# ----------------------------------------------------------------------------------------------------------------------
#   Plot Reservoir Conditions
# ----------------------------------------------------------------------------------------------------------------------  
def plot_reservoir_conditions(res_Hex,
                              save_figure = True,
                              save_filename = "Reservoir_temperature_conditions",
                              file_type = ".png"):

    # get plotting style df = pd.DataFrame({'Column1': variable1, 'Column2': variable2})
    ps = plot_style()    

    # get line colors for plots 
    line_colors = cm.inferno(np.linspace(0, 0.9, len(res_Hex.segments)))     
    fig, axis = plt.subplots()
    fig.set_size_inches(5, 3)  

    axis.grid(True)
    axis.set_xlabel('Time (mins)', fontsize = ps.axis_font_size)
    axis.set_ylabel(r'T ($\circ$C)',fontsize = ps.axis_font_size)

    for network in res_Hex.segments[0].analyses.energy.networks: 
        busses = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                axis.grid(True)

                for i in range(len(res_Hex.segments)):
                    time = res_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
                    battery_conditions = res_Hex.segments[i].conditions.energy[bus.tag][battery.tag]  
                    reservoir_temperature = battery_conditions.thermal_management_system.RES.coolant_temperature[:,0]-273

                    segment_tag = res_Hex.segments[i].tag
                    segment_name = segment_tag.replace('_', ' ') 

                    axis.plot(time, reservoir_temperature, linestyle=ps.line_style[0], markersize=ps.marker_size, color=line_colors[i], marker=ps.markers[0],  fillstyle = 'none', linewidth=ps.line_width)

    fig.tight_layout()

    if save_figure:
        plt.savefig(f'IEEE_ITEC_final_plots/{save_filename}{file_type}') 

    return fig


# ----------------------------------------------------------------------------------------------------------------------
# PLot Heat Acquisiition System
# ----------------------------------------------------------------------------------------------------------------------   

def plot_heat_acquisition_system_conditions(res_Hex,
                                            save_figure = True,
                                            show_legend = False,
                                            save_filename = "Heat_Acquisition_System",
                                            file_type = ".png"): 

    # get plotting style 
    ps      = plot_style()    

    # get line colors for plots 
    line_colors   = cm.inferno(np.linspace(0,0.9,len(res_Hex.segments)))     
    fig = plt.figure(save_filename)
    fig.set_size_inches(5,5)  
    axis_0 = plt.subplot(1,1,1)
    axis_1 = plt.subplot(2,1,1)
    axis_2 = plt.subplot(2,1,2) 

    axis_0.grid(False)
    axis_0.axis('off')         

    b_i = 0 
    for network in res_Hex.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                #axis_0.plot(np.zeros(2),np.zeros(2), color = line_colors[0], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width,label= battery.tag) 
                #axis_0.plot(np.zeros(2),np.zeros(2), color = line_colors[0], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width,label= battery.tag) 
                axis_0.grid(False)
                axis_0.axis('off')  

                for i in range(len(res_Hex.segments)):  
                    time                       = res_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
                    battery_conditions_HEX     = res_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    



                    coolant_mass_flow_rate_HEX        = battery_conditions_HEX.thermal_management_system.HAS.coolant_mass_flow_rate[:,0]
                    power_HEX                         = battery_conditions_HEX.thermal_management_system.HAS.power[:,0]     


                    segment_tag  = res_Hex.segments[i].tag
                    segment_name = segment_tag.replace('_', ' ') 

                    if b_i == 0:                     
                        axis_1.plot(time, coolant_mass_flow_rate_HEX, linestyle = ps.line_style[0], markersize= ps.marker_size , color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    else:
                        axis_1.plot(time, coolant_mass_flow_rate_HEX, linestyle = ps.line_style[0], markersize= ps.marker_size, color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_1.set_ylabel(r'$\dot{m}_c$ (kg/s)',fontsize = ps.axis_font_size)  
                    set_axes(axis_1)     


                    axis_2.plot(time, power_HEX, color = line_colors[i], linestyle = ps.line_style[0], markersize= ps.marker_size,  marker = ps.markers[1], linewidth = ps.line_width)
                    axis_2.set_ylabel(r'P$_{HAS}$ (W)',fontsize = ps.axis_font_size)
                    axis_2.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size) 
                    set_axes(axis_2)     

                b_i += 1 
    fig.tight_layout()


    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize)     

    if save_figure:
        plt.savefig(f'IEEE_ITEC_final_plots/{save_filename}{file_type}') 

    return fig

# ----------------------------------------------------------------------------------------------------------------------
# PLot Heat Exchanger System
# ----------------------------------------------------------------------------------------------------------------------   

def plot_heat_exchanger_system_conditions(res_Hex,
                                          save_figure = True,
                                          show_legend = False,
                                          save_filename = "Heat_Exchanger_System",
                                          file_type = ".png"): 

    # get plotting style 
    ps      = plot_style()    

    # get line colors for plots 
    line_colors   = cm.inferno(np.linspace(0,0.9,len(res_Hex.segments)))     
    fig = plt.figure(save_filename)
    fig.set_size_inches(5,6) 
    axis_0 = plt.subplot(1,1,1)
    axis_1 = plt.subplot(3,1,1)
    axis_2 = plt.subplot(3,1,2) 
    axis_3 = plt.subplot(3,1,3) 


    axis_0.grid(False)
    axis_0.axis('off')         

    b_i = 0 
    for network in res_Hex.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                #axis_0.plot(np.zeros(2),np.zeros(2), color = line_colors[0], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width,label= battery.tag) 
                #axis_0.plot(np.zeros(2),np.zeros(2), color = line_colors[0], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width,label= battery.tag) 
                axis_0.grid(False)
                axis_0.axis('off')  

                for i in range(len(res_Hex.segments)):  
                    time                       = res_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
                    battery_conditions_HEX     = res_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    



                    coolant_mass_flow_rate = battery_conditions_HEX.thermal_management_system.HEX.coolant_mass_flow_rate[:,0]        
                    power                  = battery_conditions_HEX.thermal_management_system.HEX.power[:,0]     
                    air_mass_flow_rate     = battery_conditions_HEX.thermal_management_system.HEX.air_mass_flow_rate[:,0]     


                    segment_tag  = res_Hex.segments[i].tag
                    segment_name = segment_tag.replace('_', ' ') 

                    if b_i == 0:                     
                        axis_1.plot(time, coolant_mass_flow_rate, linestyle = ps.line_style[0], markersize= ps.marker_size , color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    else:
                        axis_1.plot(time, coolant_mass_flow_rate, linestyle = ps.line_style[0], markersize= ps.marker_size, color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_1.set_ylabel(r'$\dot{m}_c$ (kg/s)',fontsize = ps.axis_font_size)  
                    set_axes(axis_1)     


                    axis_2.plot(time, air_mass_flow_rate, color = line_colors[i], linestyle = ps.line_style[0], markersize= ps.marker_size,  marker = ps.markers[1], linewidth = ps.line_width)
                    axis_2.set_ylabel(r'$\dot{m}_a$ (kg/s)',fontsize = ps.axis_font_size)  
                    set_axes(axis_2) 

                    axis_3.plot(time, power, color = line_colors[i], linestyle = ps.line_style[0], markersize= ps.marker_size,  marker = ps.markers[1], linewidth = ps.line_width)
                    axis_3.set_ylabel(r'P$_{HEX}$ (W)',fontsize = ps.axis_font_size)
                    axis_3.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size) 
                    set_axes(axis_3)                     


                b_i += 1 
    fig.tight_layout()


    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize)     

    if save_figure:
        plt.savefig(f'IEEE_ITEC_final_plots/{save_filename}{file_type}') 

    return fig

# ----------------------------------------------------------------------------------------------------------------------
#  Plot Battery Pack Conditions
# ---------------------------------------------------------------------------------------------------------------------- 

def plot_battery_cell_conditions_1(res_Hex,res_No_Hex,
                                   save_figure = True,
                                   show_legend = True,
                                   save_filename = "Battery_Cell_Conditions_temp",
                                   file_type = ".png",
                                  width = 12, height = 7): 

    # get plotting style 
    ps      = plot_style()   

    # get line colors for plots 
    line_colors   = cm.inferno(np.linspace(0,0.9,len(res_Hex.segments)))     

    fig = plt.figure(save_filename)
    fig = plt.figure(save_filename)
    fig.set_size_inches(5,6)   
    axis_1 = plt.subplot(2,1,1)
    axis_2 = plt.subplot(2,1,2) 



    b_i = 0 
    for network in res_Hex.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:     

                for i in range(len(res_Hex.segments)):  
                    time    = res_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
                    battery_conditions_Hex  = res_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    
                    cell_SOC_Hex            = battery_conditions_Hex.cell.state_of_charge[:,0]   
                    cell_temperature_Hex    = battery_conditions_Hex.cell.temperature[:,0]-273  

                    battery_conditions_No_Hex  = res_No_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    
                    cell_SOC_No_Hex            = battery_conditions_No_Hex.cell.state_of_charge[:,0]   
                    cell_temperature_No_Hex    = battery_conditions_No_Hex.cell.temperature[:,0]-273                      

                    segment_tag  = res_Hex.segments[i].tag
                    segment_name = segment_tag.replace('_', ' ')                     
   
                    if i == 0: 
                        axis_1.plot(time, cell_SOC_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width,label= 'No TMS')
                    else:
                        axis_1.plot(time, cell_SOC_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    if i == 0: 
                        axis_1.plot(time, cell_SOC_No_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width,label= 'TMS')
                    else:
                        axis_1.plot(time, cell_SOC_No_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                        
                    axis_1.set_ylabel(r'SOC',fontsize = ps.axis_font_size)
                    axis_1.set_ylim([0,1.1])
                    set_axes(axis_1)     

                    axis_2.plot(time, cell_temperature_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    axis_2.plot(time, cell_temperature_No_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    axis_2.set_ylabel(r'Temperature, $\degree$C',fontsize = ps.axis_font_size)
                    axis_2.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    set_axes(axis_2)     

                b_i += 1       

    fig.tight_layout()


    if show_legend:  
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.71, 0.86),loc='lower left')     

    if save_figure:
        plt.savefig(f'IEEE_ITEC_final_plots/{save_filename}{file_type}') 

    return fig

def plot_battery_cell_conditions_2(res_Hex,res_No_Hex,
                                   save_figure = True,
                                   show_legend = True,
                                   save_filename = "Battery_Cell_Conditions_Power",
                                   file_type = ".png",
                                  width = 12, height = 7): 

    # get plotting style 
    ps = plot_style()    

    # get line colors for plots 
    line_colors = cm.inferno(np.linspace(0, 0.9, len(res_Hex.segments)))     
    fig, axis = plt.subplots()
    fig.set_size_inches(5, 3)  

    axis.grid(True)
    b_i = 0 
    for network in res_Hex.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                for i in range(len(res_Hex.segments)):  
                    time    = res_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
                    battery_conditions_Hex  = res_Hex.segments[i].conditions.energy[bus.tag][battery.tag]  
                    cell_power_Hex          = battery_conditions_Hex.cell.power[:,0] 
                   

                    battery_conditions_No_Hex  = res_No_Hex.segments[i].conditions.energy[bus.tag][battery.tag]    
                    cell_power_No_Hex          = battery_conditions_No_Hex.cell.power[:,0] 
                 

                    segment_tag  = res_Hex.segments[i].tag
                    segment_name = segment_tag.replace('_', ' ')                     

                    #if b_i == 0:                     
                        #axis_1.plot(time, cell_current_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    #else:
                        #axis_1.plot(time, cell_current_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    #axis_1.plot(time, cell_current_No_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    #axis_1.set_ylabel(r'Current (A)',fontsize = ps.axis_font_size)
                    #set_axes(axis_1)  

                    if i==0:
                        axis.plot(time, cell_power_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width, label= 'No TMS')
                    else: 
                        axis.plot(time, cell_power_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    if i==0:
                        axis.plot(time, cell_power_No_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width, label= 'TMS')
                    else:
                        axis.plot(time, cell_power_No_Hex, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[1], linewidth = ps.line_width)
                    
                    axis.set_ylabel(r'Power (KW)',fontsize = ps.axis_font_size)
                    axis.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    set_axes(axis)      

                b_i += 1       

    fig.tight_layout()


    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.71, 0.725),loc='lower left')   

    if save_figure:
        plt.savefig(f'IEEE_ITEC_final_plots/{save_filename}{file_type}') 

    return fig


# ----------------------------------------------------------------------------------------------------------------------
#  Plot Cooling Drag
# ---------------------------------------------------------------------------------------------------------------------- 

def plot_cooling_drag(res_Hex,res_No_Hex,
                      save_figure = True,
                      show_legend = True,
                      save_filename = "Cooling_drag_variation",
                      file_type = ".png",
                                  width = 12, height = 7): 

     # get plotting style 
    ps = plot_style()    

    # get line colors for plots 
    line_colors = cm.inferno(np.linspace(0, 0.9, len(res_Hex.segments)))     
    fig, axis = plt.subplots()
    fig.set_size_inches(5, 3)  

    axis.grid(True)
    axis.set_xlabel('Time (mins)', fontsize = ps.axis_font_size)
    axis.set_ylabel(r'Drag (KN)',fontsize = ps.axis_font_size)


    for i in range(len(res_Hex.segments)):  

        
        time = res_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
        Drag            = -res_No_Hex.segments[i].conditions.frames.wind.drag_force_vector[:,0]
        Drag_cooling    = -res_Hex.segments[i].conditions.frames.wind.drag_force_vector[:,0]
        

       
        segment_tag  = res_Hex.segments[i].tag
        segment_name = segment_tag.replace('_', ' ')                     
        
        if i==0:
            axis.plot(time,Drag/1000 , linestyle=ps.line_style[0], markersize=ps.marker_size, color=line_colors[i], marker=ps.markers[0],  fillstyle = 'none',linewidth=ps.line_width,label= 'No TMS')
        else:
            axis.plot(time,Drag/1000 , linestyle=ps.line_style[0], markersize=ps.marker_size, color=line_colors[i], marker=ps.markers[0],  fillstyle = 'none',linewidth=ps.line_width)
        if i==0:  
            axis.plot(time,Drag_cooling/1000 , linestyle=ps.line_style[0], markersize=ps.marker_size, color=line_colors[i], marker=ps.markers[1], linewidth=ps.line_width,label= 'TMS')
        else:
            axis.plot(time,Drag_cooling/1000 , linestyle=ps.line_style[0], markersize=ps.marker_size, color=line_colors[i], marker=ps.markers[1], linewidth=ps.line_width)
            

    fig.tight_layout()


    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.15, 0.725),loc='lower left')   


    if save_figure:
        plt.savefig(f'IEEE_ITEC_final_plots/{save_filename}{file_type}') 

    return fig




# ------------------------------------------------------------------
#   Plot percentage operation 
# ------------------------------------------------------------------   
def plot_percentage_operation(res_Hex,
                              save_figure = True,
                              save_filename = "Percentage_operation_plot",
                              file_type = ".png"):
    
    # get plotting style 
    ps = plot_style()    

    #update ytick label size beacuse it is too big
    plt.rcParams.update({'ytick.labelsize': 12})
    
    # get line colors for plots 
    fig, axis = plt.subplots()
    fig.set_size_inches(10, 6)  
     
    
    for i in range(len(res_Hex.segments)):
        segment_tag = res_Hex.segments[i].tag
        segment_name = segment_tag.replace('_', ' ')     
        
        if res_Hex.segments[i].tag in ['DER', 'ICA','Climb', 'Reserve_Climb','Baseleg','Final_Approach']:
            plt.barh(segment_name,res_Hex.segments[i].percent_operation*100 , color='mediumseagreen', label='Operation of Fan')
        else:
            plt.barh(segment_name,res_Hex.segments[i].percent_operation*100 , color='c', label='Operation of Ram Inlet')
        
    plt.xlabel('Percentage of Design Power (%)',fontsize = ps.axis_font_size)
    
    plt.legend(handles=[plt.Rectangle((0,0),1,1,fc="mediumseagreen", edgecolor = 'none'), 
                        plt.Rectangle((0,0),1,1,fc="c", edgecolor = 'none')],
               labels=['Operation of Fan', 'Operation of Ram Inlet'],
               loc='center right',bbox_to_anchor=(0.95, 0.6), fontsize = ps.legend_fontsize )
    
    plt.gca().invert_yaxis()  
    plt.grid(axis='x')      
    if save_figure:
        plt.savefig(f'IEEE_ITEC_final_plots/{save_filename}{file_type}') 
    
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
