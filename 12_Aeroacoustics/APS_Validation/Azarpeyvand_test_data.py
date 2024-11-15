# Azarpeyvand_test_data.py
#
# Created:  Nov 2024, Niranjan Nanjappa

# RCAIDE Imports 

# Imports    
import RCAIDE
from RCAIDE.Framework.Core import Data

# Python Imports  
import pandas as pd
import numpy as np

def excel_column_to_number(column_letter):
    """
    Convert Excel column letter to column number (0-based).
    Example: A -> 0, B -> 1, Z -> 25, AA -> 26, AB -> 27
    
    """
    result = 0
    for char in column_letter.upper():
        result = result * 26 + (ord(char) - ord('A') + 1)
    return result - 1

def number_to_excel_column(n):
    """
    Convert 0-based column number to Excel column letter.
    Example: 0 -> A, 1 -> B, 25 -> Z, 26 -> AA, 27 -> AB
    
    """
    result = ""
    n += 1
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        result = chr(65 + remainder) + result
    return result

def excel_to_arrays(df, file_path, sheet, column_letters=None):
    """
    Read an Excel file and convert specified columns to numpy arrays based on column letters.
    
    """
    try:
        
        # If no specific columns are provided, use all columns
        if column_letters is None:
            max_col = len(df.columns)
            column_letters = [number_to_excel_column(i) for i in range(max_col)]
        
        # Convert column letters to numbers and check validity
        column_numbers = []
        for col_letter in column_letters:
            col_num = excel_column_to_number(col_letter)
            if col_num < 0 or col_num >= len(df.columns):
                raise ValueError(
                    f"Invalid column letter: {col_letter}. "
                    f"Valid range is A to {number_to_excel_column(len(df.columns)-1)}"
                )
            column_numbers.append(col_num)
        
        # Convert specified columns to numpy arrays
        arrays_dict = np.array([])
        for col_num in column_numbers:
            arrays_dict = df.iloc[:, col_num].to_numpy()
        
        return arrays_dict
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None
    except Exception as e:
        print(f"Error: {str(e)}")
        return None

def get_third_octave_bands():
    """
    Returns the nominal 1/3 octave band center frequencies
    """
    return np.array([
        16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500,
        630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,
        10000, 12500, 16000, 20000, 25000, 31500
    ])

def get_band_limits(center_freq):
    """
    Calculate lower and upper frequency limits for each 1/3 octave band
    """
    factor = 2 ** (1/6)  # One-sixth of an octave
    lower = center_freq / factor
    upper = center_freq * factor
    return lower, upper

def convert_to_third_octave(spl, frequencies):
    """
    Convert narrowband SPL data to 1/3 octave bands
    """
    # Remove any data points where frequency is 0 Hz
    mask = frequencies > 0
    spl = spl[mask]
    frequencies = frequencies[mask]
    
    # Get third octave bands
    third_octave_bands = get_third_octave_bands()
    third_octave_levels = np.zeros(len(third_octave_bands))
    lower_freq, upper_freq = get_band_limits(third_octave_bands)
    
    for i in range(len(third_octave_bands)):
        
        # Find indices of frequencies within band limits
        mask = (frequencies >= lower_freq[i]) & (frequencies < upper_freq[i])
        
        if np.any(mask):
            # Handle NaN values in the SPL data
            band_spl = spl[mask]
            valid_data = ~np.isnan(band_spl)
            if np.any(valid_data):
                # Convert valid SPL values to pressure squared
                pressure_squared = 10.0 ** (band_spl[valid_data] / 10.0)
                
                # Calculate mean square pressure in the band
                mean_square_pressure = np.mean(pressure_squared)
                
                # Convert back to SPL
                third_octave_levels[i] = 10.0 * np.log10(mean_square_pressure)
            else:
                third_octave_levels[i] = np.nan
        else:
            third_octave_levels[i] = np.nan
    
    return third_octave_levels

def process_data_structure(data_dict):
    """
    Process the entire nested data structure
    """
    
    frequencies = data_dict.f
    result = Data()
    result.f = get_third_octave_bands()
    
    for velocity_key in ['V__8_7', 'V__26_5']:
        result[velocity_key] = Data()
        
        for alpha_key in ['alpha0', 'alpha30']:
            result[velocity_key][alpha_key] = Data()
            
            for theta_key in ['Theta60', 'Theta90', 'Theta120']:
                spl_data = data_dict[velocity_key][alpha_key][theta_key]
                # Convert to 1/3 octave bands
                third_octave_data = convert_to_third_octave(spl_data, frequencies)
                result[velocity_key][alpha_key][theta_key] = third_octave_data
    
    return result

def Azarpeyvand_test_data():
    Columns = Data()
    Columns.V__8_7 = Data()
    Columns.V__26_5 = Data()
    Columns.V__8_7.alpha0 = Data()
    Columns.V__8_7.alpha30 = Data()
    Columns.V__26_5.alpha0 = Data()
    Columns.V__26_5.alpha30 = Data()
    
    file_path = "Azarpeyvand_SPL_test.xlsx"
    
    column_60_0   = "B"
    column_60_30  = "G"
    column_90_0   = "I"
    column_90_30  = "N"
    column_120_0  = "P"
    column_120_30 = "U"
    
    # Read excel file
    sheet_1 = "8.7ms"
    
    file = pd.read_excel(file_path, sheet_name=sheet_1, header=1)
    
    Columns.V__8_7.alpha0.Theta60    = excel_to_arrays(file, file_path, sheet_1, column_60_0)
    Columns.V__8_7.alpha30.Theta60   = excel_to_arrays(file, file_path, sheet_1, column_60_30)
    Columns.V__8_7.alpha0.Theta90    = excel_to_arrays(file, file_path, sheet_1, column_90_0)
    Columns.V__8_7.alpha30.Theta90   = excel_to_arrays(file, file_path, sheet_1, column_90_30)
    Columns.V__8_7.alpha0.Theta120   = excel_to_arrays(file, file_path, sheet_1, column_120_0)
    Columns.V__8_7.alpha30.Theta120  = excel_to_arrays(file, file_path, sheet_1, column_120_30)
    
    sheet_2 = "26.5ms"
    
    file = pd.read_excel(file_path, sheet_name=sheet_2, header=1)
    
    Columns.V__26_5.alpha0.Theta60    = excel_to_arrays(file, file_path, sheet_2, column_60_0)
    Columns.V__26_5.alpha30.Theta60   = excel_to_arrays(file, file_path, sheet_2, column_60_30)
    Columns.V__26_5.alpha0.Theta90    = excel_to_arrays(file, file_path, sheet_2, column_90_0)
    Columns.V__26_5.alpha30.Theta90   = excel_to_arrays(file, file_path, sheet_2, column_90_30)
    Columns.V__26_5.alpha0.Theta120   = excel_to_arrays(file, file_path, sheet_2, column_120_0)
    Columns.V__26_5.alpha30.Theta120  = excel_to_arrays(file, file_path, sheet_2, column_120_30)
    
    fcolumn = "A"
    Columns.f = excel_to_arrays(file, file_path, sheet_2, fcolumn)
    
    # Process the data
    Azarpeyvand_data = process_data_structure(Columns)

    return Azarpeyvand_data
