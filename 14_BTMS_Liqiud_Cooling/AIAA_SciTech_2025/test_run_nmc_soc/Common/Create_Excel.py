import pandas as pd
from openpyxl import load_workbook
import os
from RCAIDE.Framework.Core import Units
import numpy as np
def write_data(results, filename, group):
    
    # Create or append to Excel file with group sheets

    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Go one folder back and into Raw_Data
    load_dir = os.path.join(current_dir, '..', 'Raw_Data')
    
    excel_file = os.path.join(load_dir, filename + '.xlsx')
    sheet_name = 'Group_' + str(group)
    
    if not os.path.exists(excel_file):
        # If the file doesn't exist, create a new DataFrame
        df = pd.DataFrame()  # Add any initial data here
        # Save the DataFrame to an Excel file
        df.to_excel(excel_file, index=False, engine='openpyxl')
    
    existing_df = pd.DataFrame()  
  
    # Loop through the networks in the first segment (assuming segments[0] is representative)
    for network in results.segments[0].analyses.energy.vehicle.networks: 
        busses = network.busses
        coolant_lines = network.coolant_lines
        # Zip busses and coolant lines together
        for bus, coolant_line in zip(busses, coolant_lines):
           
            # Loop through all segments and collect the data
            for i in range(len(results.segments)):
                time = results.segments[i].conditions.frames.inertial.time[:, 0] / Units.min    
                # Battery Conditions
                battery_conditions = results.segments[i].conditions.energy[bus.tag].battery_modules.lithium_ion_nmc
                cell_power         = battery_conditions.cell.power[:, 0]
                cell_energy        = battery_conditions.cell.energy[:, 0]
                cell_volts         = battery_conditions.cell.voltage_under_load[:, 0] 
                cell_current       = battery_conditions.cell.current[:, 0]
                cell_SOC           = battery_conditions.cell.state_of_charge[:, 0]   
                cell_temperature   = battery_conditions.cell.temperature[:, 0]
                module_heat_gen    = battery_conditions.heat_energy_generated[:, 0]
                cycle_day          = battery_conditions.cell.cycle_in_day*np.ones_like(time)
                capacity_fade      = battery_conditions.cell.capacity_fade_factor*np.ones_like(time)
                resistance_growth  = battery_conditions.cell.resistance_growth_factor*np.ones_like(time)
                charge_throughput  = battery_conditions.cell.charge_throughput[:, 0]

                # Flight Conditions
                airspeed           = results.segments[i].conditions.freestream.velocity[:,0] /   Units['mph']
                theta              = results.segments[i].conditions.frames.body.inertial_rotations[:,1,None] / Units.deg
                Range              = results.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.nmi
                altitude           = results.segments[i].conditions.freestream.altitude[:,0]/Units.feet

                # Reservoir Conditions
                reservoir_conditions    = results.segments[i].conditions.energy[coolant_line.tag].coolant_reservoir   
                reservoir_temperature =  reservoir_conditions.coolant_temperature[:,0]
                
                # HEX Conditions
                cross_flow_hex_conditions = results.segments[i].conditions.energy[coolant_line.tag].cross_flow_heat_exchanger
                coolant_mass_flow_rate_HEX     = cross_flow_hex_conditions.coolant_mass_flow_rate[:,0]        
                effectiveness_HEX          = cross_flow_hex_conditions.effectiveness_HEX[:,0]   
                power_HEX                  = cross_flow_hex_conditions.power[:,0]                       
                inlet_air_pressure_HEX         = cross_flow_hex_conditions.air_inlet_pressure[:,0]          
                outlet_coolant_temperature_HEX = cross_flow_hex_conditions.outlet_coolant_temperature[:,0]          
                air_mass_flow_rate_HEX         = cross_flow_hex_conditions.air_mass_flow_rate[:,0]     
                
                # Wavy Channel Conditions
                wavy_channel_conditions         = results.segments[i].conditions.energy[coolant_line.tag].wavy_channel_heat_acquisition
                outlet_coolant_temperature_wavy_channel      = wavy_channel_conditions.outlet_coolant_temperature[:,0]
                coolant_mass_flow_rate_wavy_channel          = wavy_channel_conditions.coolant_mass_flow_rate[:,0]
                power_wavy_channel              = wavy_channel_conditions.power[:,0]         
    
                data = {
                    'Time (min)': time,
                    'Cell Power (W)': cell_power,
                    'Cell Energy (J)': cell_energy,
                    'Cell Voltage (V)': cell_volts,
                    'Cell Current (A)': cell_current,
                    'Cell SOC (%)': cell_SOC,
                    'Cell Temperature (C)': cell_temperature,
                    'Cycle Day': cycle_day,
                    'Capacity Fade': capacity_fade,
                    'Module Heat Generation (W)': module_heat_gen,
                    'Resistance Growth': resistance_growth,
                    'Charge Throughput (Ah)': charge_throughput,
                    'Airspeed (mph)': airspeed,
                    'Range (nmi)': Range,
                    'Altitude (ft)': altitude,
                    'Reservoir Temperature (K)': reservoir_temperature,
                    'Coolant Mass Flow Rate HEX (kg/s)': coolant_mass_flow_rate_HEX,
                    'Effectiveness HEX': effectiveness_HEX,
                    'Power HEX (W)': power_HEX,
                    'Inlet Air Pressure HEX (Pa)': inlet_air_pressure_HEX,
                    'Outlet Coolant Temperature HEX (K)': outlet_coolant_temperature_HEX,
                    'Air Mass Flow Rate HEX (kg/s)': air_mass_flow_rate_HEX, 
                    'Outlet Coolant Temperature Wavy Channel (K)': outlet_coolant_temperature_wavy_channel,
                    'Coolant Mass Flow Rate Wavy Channel (kg/s)': coolant_mass_flow_rate_wavy_channel,
                    'Power Wavy Channel (W)': power_wavy_channel,
            
                    }

                df = pd.DataFrame(data)
                existing_df = pd.concat([existing_df, df], ignore_index=True)
               

    with pd.ExcelWriter(excel_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        # This automatically removes the old sheet and replaces it with new data
        existing_df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"Data successfully written to sheet '{sheet_name}'.")
    return
