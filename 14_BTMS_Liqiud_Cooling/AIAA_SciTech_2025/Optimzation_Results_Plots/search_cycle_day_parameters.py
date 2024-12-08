import pandas as pd
import os

def main():
    # Directory containing the result files
    directory = os.path.dirname(__file__)
    file_prefix1 = 'Vehicle_Excel_Files/case_'
    file_prefix2  = '/Raw_Data/e_Twin_Otter_nmc_case_'
    file_extension = '.xlsx'
    num_files = 75  # Number of files to process (0 to 74)

    # List to collect the consolidated data
    consolidated_data = []

    for i in range(num_files):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        load_dir = os.path.join(current_dir)
        file_name = f"{file_prefix1}{i}{file_prefix2}{i}.0{file_extension}"
        file_path = os.path.join(load_dir, file_name)
       # file_path = os.path.join(os.path.dirname(__file__), file_name)

        if os.path.exists(file_path):
            print(f"Processing file: {file_name}")
            case_id = f"{i}.0"  # Extract Case ID
            data = extract_exit_condition_from_last_sheet(file_path, case_id)
            if data:
                consolidated_data.append(data)
        else:
            print(f"File not found: {file_name}")

    # Save consolidated data to an Excel file
    consolidated_df = pd.DataFrame(consolidated_data, columns=["Case ID", "Cycle Time", "Cycle Day", "Exit Condition"])
    consolidated_output_path = os.path.join(directory, "consolidated_exit_conditions.xlsx")
    consolidated_df.to_excel(consolidated_output_path, index=False)
    print(f"Consolidated data saved to: {consolidated_output_path}")

def extract_exit_condition_from_last_sheet(file_path, case_id):
    # Open the Excel file
    excel_file = pd.ExcelFile(file_path)
    sheet_names = excel_file.sheet_names
    cumulative_time = 0
    exit_condition_data = None

    for sheet in sheet_names:
        # Read each sheet
        df = pd.read_excel(excel_file, sheet_name=sheet)
        if 'Time (min)' in df.columns:
            # Update cumulative time
            df['Cumulative Time'] = df['Time (min)'] + cumulative_time
            cumulative_time = df['Cumulative Time'].iloc[-1]

    # Process the last sheet for the exit condition
    last_sheet = sheet_names[-1]
    df_last = pd.read_excel(excel_file, sheet_name=last_sheet)
    if 'Cell SOC (%)' in df_last.columns and 'Cell Temperature (C)' in df_last.columns and 'Cycle Day' in df_last.columns:
        # Define exit condition
        soc_condition = df_last['Cell SOC (%)'] < 0.2
        temp_condition = df_last['Cell Temperature (C)'] >= 322.65
        condition = soc_condition | temp_condition

        # Find the first occurrence of the condition
        if condition.any():
            cycle_day = df_last.loc[condition, 'Cycle Day'].min()
            exit_condition = "SOC < 0.2" if df_last.loc[soc_condition, 'Cycle Day'].min() <= df_last.loc[temp_condition, 'Cycle Day'].min() else "Temperature >= 322.65"
            exit_condition_data = [case_id, cumulative_time, cycle_day, exit_condition]

    return exit_condition_data

if __name__ == "__main__":
    main()
 