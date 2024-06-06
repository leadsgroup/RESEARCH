# RCAIDE imports 
import RCAIDE  
from RCAIDE.Core                                    import Units 
from RCAIDE.Methods.Energy.Sources.Battery.Common   import initialize_from_circuit_configuration   
import matplotlib.pyplot as plt 
from RCAIDE.Core import Units
from RCAIDE.Visualization.Common import set_axes, plot_style 
import matplotlib.cm as cm
import numpy as np 
 
# ----------------------------------------------------------------------------------------------------------------------
#  Build the Vehicle
# ----------------------------------------------------------------------------------------------------------------------   
def main():   
    vehicle  = vehicle_setup() 
    
    # Set up vehicle configs
    configs  = configs_setup(vehicle)
    
    # create analyses
    analyses = analyses_setup(configs)
    
    # mission analyses
    C_rate   = 2
    Ah       = vehicle.networks.isolated_battery_cell.busses.bus.batteries.lithium_ion_nmc.cell.nominal_capacity 
    mission  = mission_setup(analyses,vehicle,C_rate,Ah ) 
    
    # create mission instances (for multiple types of missions)
    missions = missions_setup(mission) 
     
    # mission analysis 
    results = missions.base_mission.evaluate()    

    plot_battery_cell_conditions(results)    
    
    return  

def analyses_setup(configs):

    analyses = RCAIDE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):    
    #   Initialize the Analyses     
    analyses = RCAIDE.Analyses.Vehicle()  
    
    #  Energy
    energy          = RCAIDE.Analyses.Energy.Energy()
    energy.networks = vehicle.networks 
    analyses.append(energy)
 
    #  Planet Analysis
    planet  = RCAIDE.Analyses.Planets.Planet()
    analyses.append(planet)
 
    #  Atmosphere Analysis
    atmosphere                 = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   
 
    return analyses     

def mission_setup(analyses,vehicle,C_rate,Ah):
 
    #   Initialize the Mission 
    mission            = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag        = 'cell_cycle_test'   
    Segments           = RCAIDE.Analyses.Mission.Segments 
    base_segment       = Segments.Segment()   
    time               = (1/C_rate) * Units.hrs  
    current            = C_rate*Ah
        
    # Discharge Segment 
    segment                                 = Segments.Ground.Battery_Discharge(base_segment) 
    segment.analyses.extend(analyses.base)  
    segment.tag                             = 'Discharge_1' 
    segment.time                            = time 
    segment.current                         = current
    segment.battery_cell_temperature        = 300   
    segment.initial_battery_state_of_charge = 1  
    mission.append_segment(segment)          
    
    ## Charge Segment 
    #segment                                = Segments.Ground.Battery_Recharge(base_segment)      
    #segment.analyses.extend(analyses.base) 
    #segment.tag                            = 'Recharge'
    #segment.current                        = current
    #segment.time                           = time     
    #mission.append_segment(segment)   

    return mission 

def vehicle_setup(): 

    vehicle                                 = RCAIDE.Vehicle() 
    vehicle.tag                             = 'battery'   
    vehicle.reference_area                  = 1
  
    # mass properties
    vehicle.mass_properties.takeoff         = 1 * Units.kg 
    vehicle.mass_properties.max_takeoff     = 1 * Units.kg 
         
    net                                     = RCAIDE.Energy.Networks.Isolated_Battery_Cell_Network() 
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Bus
    #------------------------------------------------------------------------------------------------------------------------------------  
    bus                                                       = RCAIDE.Energy.Networks.Distribution.Electrical_Bus()  
    battery                                                   = RCAIDE.Energy.Sources.Batteries.Lithium_Ion_NMC() 
    HAS                                                       = RCAIDE.Energy.Thermal_Management.Batteries.Heat_Acquisition_Systems.Direct_Air() 
    HAS.convective_heat_transfer_coefficient                  = 7.17
    battery.thermal_management_system.heat_acquisition_system = HAS 
    net.voltage                                               = battery.cell.nominal_voltage 
    initialize_from_circuit_configuration(battery,module_weight_factor = 1.0)  
    bus.voltage                                               = battery.pack.maximum_voltage  
    bus.batteries.append(battery)                                
      
    # append bus   
    net.busses.append(bus) 
    
    # append network 
    vehicle.append_energy_network(net) 
      
    vehicle.mass_properties.takeoff = battery.mass_properties.mass 
    return vehicle


def missions_setup(mission): 
 
    missions         = RCAIDE.Analyses.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions  

def configs_setup(vehicle): 
    configs         = RCAIDE.Components.Configs.Config.Container()  
    base_config     = RCAIDE.Components.Configs.Config(vehicle)
    base_config.tag = 'base' 
    configs.append(base_config)   
    return configs  
# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Visualization-Performance-Energy-Battery
def plot_battery_cell_conditions(results,
                                  save_figure = False,
                                  show_legend = True,
                                  save_filename = "Battery_Cell_Conditions_",
                                  file_type = ".png",
                                  width = 12, height = 7):
    """Plots the cell-level conditions of the battery throughout flight.

    Assumptions:
    None

    Source:
    None

    Inputs:
    results.segments.conditions.
        freestream.altitude
        weights.total_mass
        weights.vehicle_mass_rate
        frames.body.thrust_force_vector

    Outputs: 
    Plots

    Properties Used:
    N/A	
    """ 
    
    # get plotting style 
    ps      = plot_style()  

    parameters = {'axes.labelsize': ps.axis_font_size,
                  'xtick.labelsize': ps.axis_font_size,
                  'ytick.labelsize': ps.axis_font_size,
                  'axes.titlesize': ps.title_font_size}
    plt.rcParams.update(parameters)
     
    # get line colors for plots 
    line_colors   = cm.inferno(np.linspace(0,0.9,len(results.segments)))     

    fig = plt.figure(save_filename)
    fig.set_size_inches(width,height) 
    axis_0 = plt.subplot(1,1,1)
    axis_1 = plt.subplot(3,2,1)
    axis_2 = plt.subplot(3,2,2) 
    axis_3 = plt.subplot(3,2,3) 
    axis_4 = plt.subplot(3,2,4)
    axis_5 = plt.subplot(3,2,5) 
    axis_6 = plt.subplot(3,2,6)      
    b_i = 0 
    for network in results.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                axis_0.plot(np.zeros(2),np.zeros(2), color = line_colors[0], marker = ps.markers[b_i], linewidth = ps.line_width,label= battery.tag) 
                axis_0.grid(False)
                axis_0.axis('off')  
               
                for i in range(len(results.segments)):  
                    time    = results.segments[i].conditions.frames.inertial.time[:,0]    
                    battery_conditions  = results.segments[i].conditions.energy[bus.tag][battery.tag]    
                    cell_power          = battery_conditions.cell.power[:,0]
                    cell_energy         = battery_conditions.cell.energy[:,0]
                    cell_volts          = battery_conditions.cell.voltage_under_load[:,0]
                    cell_volts_oc       = battery_conditions.cell.voltage_open_circuit[:,0]
                    cell_current        = battery_conditions.cell.current[:,0]
                    cell_SOC            = battery_conditions.cell.state_of_charge[:,0]   
                    cell_temperature    = battery_conditions.cell.temperature[:,0]  
            
                    segment_tag  = results.segments[i].tag
                    segment_name = segment_tag.replace('_', ' ') 

                    if b_i == 0:                     
                        axis_1.plot(time, cell_SOC, color = line_colors[i], marker = ps.markers[b_i], linewidth = ps.line_width, label = segment_name)
                    else:
                        axis_1.plot(time, cell_SOC, color = line_colors[i], marker = ps.markers[b_i], linewidth = ps.line_width)
                    axis_1.set_ylabel(r'SOC')
                    axis_1.set_ylim([0,1.1])
                    set_axes(axis_1)     
                     
                    axis_2.plot(time, cell_energy/Units.Wh, color = line_colors[i], marker = ps.markers[b_i], linewidth = ps.line_width)
                    axis_2.set_ylabel(r'Energy (W-hr)')
                    set_axes(axis_2) 
             
                    axis_3.plot(time, cell_current, color = line_colors[i], marker = ps.markers[b_i], linewidth = ps.line_width)
                    axis_3.set_ylabel(r'Current (A)')
                    set_axes(axis_3)  
             
                    axis_4.plot(time, cell_power, color = line_colors[i], marker = ps.markers[b_i], linewidth = ps.line_width)
                    axis_4.set_ylabel(r'Power (W)')
                    set_axes(axis_4)     
                     
                    axis_5.plot(time, cell_volts, color = line_colors[i], marker = ps.markers[b_i], linewidth = ps.line_width) 
                    axis_5.set_ylabel(r'Voltage (V)')
                    axis_5.set_xlabel(r'Time (s)')
                    set_axes(axis_5) 
             
                    axis_6.plot(time, cell_temperature, color = line_colors[i], marker = ps.markers[b_i], linewidth = ps.line_width)
                    axis_6.set_ylabel(r'Temperature, $\degree$C')
                    axis_6.set_xlabel(r'Time (s)')
                    set_axes(axis_6)  
                b_i += 1 
            
    if show_legend:      
        h, l = axis_0.get_legend_handles_labels()
        axis_2.legend(h, l)    
        leg =  fig.legend(bbox_to_anchor=(0.5, 0.95), loc='upper center', ncol = 5) 
        leg.set_title('Flight Segment', prop={'size': ps.legend_font_size, 'weight': 'heavy'})      
    
    # Adjusting the sub-plots for legend 
    fig.subplots_adjust(top=0.8) 
    
    # set title of plot 
    title_text   = 'Battery Cell Conditions'       
    fig.suptitle(title_text) 
    
    if save_figure:
        plt.savefig(save_filename + battery.tag + file_type)    
    return fig 

if __name__ == '__main__':
    main()
    plt.show()