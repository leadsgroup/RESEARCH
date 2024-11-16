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
    aggregated_df = pd.DataFrame()
  
    # Loop through the networks in the first segment (assuming segments[0] is representative)
    for network in results.segments[0].analyses.energy.vehicle.networks: 
        busses = network.busses
        coolant_lines = network.coolant_lines
        # Zip busses and coolant lines together
        for bus in busses:
            for _,item in enumerate(bus.battery_modules):
                battery_tag = item.tag
                break
          
            # Loop through all segments and collect the data
            for i in range(len(results.segments)):
                time = results.segments[i].conditions.frames.inertial.time[:, 0] / Units.min    

                # Battery Conditions
                battery_conditions = results.segments[i].conditions.energy[bus.tag].battery_modules[battery_tag]
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
                reservoir_conditions    = results.segments[i].conditions.energy.coolant_line.coolant_reservoir   
                reservoir_temperature =  reservoir_conditions.coolant_temperature[:,0]
                
                # HEX Conditions
                cross_flow_hex_conditions = results.segments[i].conditions.energy.coolant_line.cross_flow_heat_exchanger
                coolant_mass_flow_rate_HEX     = cross_flow_hex_conditions.coolant_mass_flow_rate[:,0]        
                effectiveness_HEX          = cross_flow_hex_conditions.effectiveness_HEX[:,0]   
                power_HEX                  = cross_flow_hex_conditions.power[:,0]                       
                inlet_air_pressure_HEX         = cross_flow_hex_conditions.air_inlet_pressure[:,0]          
                outlet_coolant_temperature_HEX = cross_flow_hex_conditions.outlet_coolant_temperature[:,0]          
                air_mass_flow_rate_HEX         = cross_flow_hex_conditions.air_mass_flow_rate[:,0]     
                
                # Wavy Channel Conditions
                wavy_channel_conditions         = results.segments[i].conditions.energy.coolant_line.wavy_channel_heat_acquisition
                outlet_coolant_temperature_wavy_channel      = wavy_channel_conditions.outlet_coolant_temperature[:,0]
                coolant_mass_flow_rate_wavy_channel          = wavy_channel_conditions.coolant_mass_flow_rate[:,0]
                power_wavy_channel              = wavy_channel_conditions.power[:,0]         

                
                data = {
                f'Time (min) [{bus.tag}]': time,
                f'Cell Power (W) [{bus.tag}]': cell_power,
                f'Cell Energy (J) [{bus.tag}]': cell_energy,
                f'Cell Voltage (V) [{bus.tag}]': cell_volts,
                f'Cell Current (A) [{bus.tag}]': cell_current,
                f'Cell SOC (%) [{bus.tag}]': cell_SOC,
                f'Cell Temperature (C) [{bus.tag}]': cell_temperature,
                f'Cycle Day [{bus.tag}]': cycle_day,
                f'Capacity Fade [{bus.tag}]': capacity_fade,
                f'Module Heat Generation (W) [{bus.tag}]': module_heat_gen,
                f'Resistance Growth [{bus.tag}]': resistance_growth,
                f'Charge Throughput (Ah) [{bus.tag}]': charge_throughput,
                f'Airspeed (mph) [{bus.tag}]': airspeed,
                f'Range (nmi) [{bus.tag}]': Range,
                f'Altitude (ft) [{bus.tag}]': altitude,
                f'Reservoir Temperature (K) [{bus.tag}]': reservoir_temperature,
                f'Coolant Mass Flow Rate HEX (kg/s) [{bus.tag}]': coolant_mass_flow_rate_HEX,
                f'Effectiveness HEX [{bus.tag}]': effectiveness_HEX,
                f'Power HEX (W) [{bus.tag}]': power_HEX,
                f'Inlet Air Pressure HEX (Pa) [{bus.tag}]': inlet_air_pressure_HEX,
                f'Outlet Coolant Temperature HEX (K) [{bus.tag}]': outlet_coolant_temperature_HEX,
                f'Air Mass Flow Rate HEX (kg/s) [{bus.tag}]': air_mass_flow_rate_HEX,
                f'Outlet Coolant Temperature Wavy Channel (K) [{bus.tag}]': outlet_coolant_temperature_wavy_channel,
                f'Coolant Mass Flow Rate Wavy Channel (kg/s) [{bus.tag}]': coolant_mass_flow_rate_wavy_channel,
                f'Power Wavy Channel (W) [{bus.tag}]': power_wavy_channel,
                }


            # Convert the data dictionary into a DataFrame
                df = pd.DataFrame(data)
                aggregated_df = pd.concat([aggregated_df, df], ignore_index=True)

            # Append the aggregated_df data as new columns to existing_df
            existing_df = append_columns(existing_df, aggregated_df)
            aggregated_df = pd.DataFrame()  # Reset for the next bus
                

    with pd.ExcelWriter(excel_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        # This automatically removes the old sheet and replaces it with new data
        existing_df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"Data successfully written to sheet '{sheet_name}'.")
    return

def append_columns(existing_df, new_df):
    if existing_df.empty:
        # If the first DataFrame is empty, simply return the new DataFrame
        return new_df
    else:
        # Create an index-aligned DataFrame with NaN values for `new_df` to match `existing_df` index
        new_df_aligned = pd.DataFrame(index=existing_df.index)
        
        # Add the columns from `new_df` to the aligned DataFrame
        for col in new_df.columns:
            new_df_aligned[col] = new_df[col].values
        
        # Concatenate `existing_df` and `new_df_aligned` horizontally
        result_df = pd.concat([existing_df, new_df_aligned], axis=1)
        return result_df

