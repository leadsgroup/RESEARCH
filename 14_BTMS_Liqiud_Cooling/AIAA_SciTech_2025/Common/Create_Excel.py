import pandas as pd
from openpyxl import load_workbook
import os
from RCAIDE.Framework.Core import Units

def write_data(results, filename, group):
    
    # Create or append to Excel file with group sheets
    excel_file = os.path.join("Raw_Data", filename + '.xlsx')
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
        for bus in busses: 
            for b_i, battery in enumerate(bus.battery_modules):
                if b_i == 0 or bus.identical_batteries == False:
                    # Loop through all segments and collect the data
                    for i in range(len(results.segments)):
                        time = results.segments[i].conditions.frames.inertial.time[:, 0] / Units.min    
                        battery_conditions = results.segments[i].conditions.energy[bus.tag].battery_modules[battery.tag]    
                        cell_power = battery_conditions.cell.power[:, 0]
                        cell_energy = battery_conditions.cell.energy[:, 0]
                        cell_volts = battery_conditions.cell.voltage_under_load[:, 0] 
                        cell_current = battery_conditions.cell.current[:, 0]
                        cell_SOC = battery_conditions.cell.state_of_charge[:, 0]   
                        cell_temperature = battery_conditions.cell.temperature[:, 0]
    
                        # Create a DataFrame for the current segment with the collected data
                        data = {
                            'Time (min)': time,
                            'Cell Power (W)': cell_power,
                            'Cell Energy (J)': cell_energy,
                            'Cell Voltage (V)': cell_volts,
                            'Cell Current (A)': cell_current,
                            'Cell SOC (%)': cell_SOC,
                            'Cell Temperature (C)': cell_temperature
                        }
    
                        df = pd.DataFrame(data)
    
                        # Append the new data to the existing data
                        existing_df = pd.concat([existing_df, df], ignore_index=True)
                        
    with pd.ExcelWriter(excel_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        # This automatically removes the old sheet and replaces it with new data
        existing_df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"Data successfully written to sheet '{sheet_name}'.")
    return
