import os
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

def main(): 
    file_name = 'e_Twin_Otter_nmc_case_2_.xlsx'
    file_path = os.path.join(os.path.dirname(__file__), file_name)
    df = read_excel_file(file_path)
    variables = [
        ('Time (min)', 'Cell Temperature (C)', 'Cell Temperature vs Time'),
        ('Time (min)', 'Cell SOC (%)', 'Cell SOC vs Time'),
        ('Cycle Day', 'Resistance Growth', 'Resistance Growth vs Cycle Day'),
        ('Cycle Day', 'Capacity Fade', 'Capacity Fade vs Cycle Day'),
        ('Cycle Day', 'Cell Temperature (C)', 'Cell Temperature vs Cycle Day'),
        ('Time (min)', 'Range (nmi)', 'Range vs Time'),
        ('Time (min)', 'Module Heat Generation (W)', 'Module Heat Generation vs Time')
    ]
    plot_data(df, variables)
    return

def read_excel_file(file_path):
    # Read the Excel file
    excel_file = pd.ExcelFile(file_path)
    sheet_names = excel_file.sheet_names
    print("\nAvailable sheets in the Excel file:")
    
    # Dictionary to store concatenated data
    combined_data = {}
    last_time = 0  # Keep track of the last time value
    
    # Iterate through each sheet
    for sheet in sheet_names:
        print(f"- {sheet}")
        # Read the sheet into a DataFrame
        df = pd.read_excel(excel_file, sheet_name=sheet)
        
        # Add the last_time to the current sheet's time values
        if 'Time (min)' in df.columns:
            df['Time (min)'] = df['Time (min)'] + last_time
            # Update last_time for the next sheet
            last_time = df['Time (min)'].iloc[-1]
        
        # For each column in the current sheet
        for column in df.columns:
            # If column name already exists in combined_data, concatenate
            if column in combined_data:
                combined_data[column] = pd.concat([combined_data[column], df[column]], ignore_index=True)
            # If it's a new column name, add it to the dictionary
            else:
                combined_data[column] = df[column]
    
    # Convert the dictionary to a DataFrame
    result_df = pd.DataFrame(combined_data)
    return result_df

def plot_data(df, variables):
    for x_var, y_var, title_suffix in variables:
        fig = go.Figure()

        # Add the line trace
        fig.add_trace(go.Scatter(
            x=df[x_var],
            y=df[y_var],
            mode='lines',
            name=y_var,
            line=dict(color='royalblue', width=2),
            fill='tozeroy'  ,# Adds shading under the line
            fillcolor='rgba(255, 0, 0, 0.3)'
        ))

        # Update layout for styling
        fig.update_layout(
            title=title_suffix,
            xaxis_title=x_var,
            yaxis_title=y_var,
            template="simple_white",
            title_font_size=20,
            legend=dict(title_font_family="Times New Roman", font_size=14),
            xaxis=dict(showgrid=True, gridwidth=0.5),
            yaxis=dict(showgrid=True, gridwidth=0.5),
        )

        # Show plot
        fig.show()

if __name__ == "__main__":
    main()
