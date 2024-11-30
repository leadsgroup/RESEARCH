import os
import pandas as pd
import re

def extract_details_from_log(file_path):
    extracted_data = []
    pattern = r"HAS_power:\s([\d\.]+),\sHEX_power:\s([\d\.]+),\sRES_dimensions:\s([\d\.]+)"
    
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                has_power, hex_power, res_dimensions = map(float, match.groups())
                extracted_data.append({
                    "HAS_power": has_power,
                    "HEX_power": hex_power,
                    "RES_dimensions": res_dimensions
                })
    
    return pd.DataFrame(extracted_data)

def extract_weight_details(file_path):
    # Patterns for extracting weight sections and total weight
    section_start_pattern = r"propulsion :"
    total_weight_pattern = r"TMS : ([\d\.]+)$"
    
    extracting = False
    weights = []
    
    with open(file_path, 'r') as file:
        for line in file:
            if section_start_pattern in line:
                extracting = True  # Start extracting
            if extracting:
                match = re.search(total_weight_pattern, line)
                if match:
                    weights.append(float(match.group(1)))
    
    return weights[-1] if weights else None  # Return the last total weight if found

def process_all_logs(directory, file_prefix, file_extension, num_files):
    all_details = []
    for i in range(num_files):
        file_name = f"{file_prefix}{i}{file_extension}"
        file_path = os.path.join(directory, file_name)
        
        if os.path.exists(file_path):
            details_df = extract_details_from_log(file_path)
            total_weight = extract_weight_details(file_path)
            
            if not details_df.empty:
                details_df['File'] = file_name
                details_df['Total Weight'] = total_weight  # Add total weight to each row
                all_details.append(details_df)
        else:
            print(f"File not found: {file_name}")
    
    if all_details:
        return pd.concat(all_details, ignore_index=True)
    return pd.DataFrame()

def main():
    directory = '/Users/sai/Documents/Research/LEADS_Research/RESEARCH/14_BTMS_Liqiud_Cooling/AIAA_SciTech_2025/Optimization_Results_Analysis/outputs'
    file_prefix = 'output_'
    file_extension = '.log'
    num_files = 75
    
    combined_details_df = process_all_logs(directory, file_prefix, file_extension, num_files)
    
    if not combined_details_df.empty:
        # Save the data to a CSV file
        output_csv_path = os.path.join(directory, "consolidated_details.csv")
        combined_details_df.to_csv(output_csv_path, index=False)
        print(f"Data extracted and saved to {output_csv_path}")
    else:
        print("No data extracted from the logs.")

if __name__ == "__main__":
    main()
