#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import pickle
import os
import sys

import RCAIDE
from RCAIDE.Framework.Core import Units, Data   
from RCAIDE.Library.Plots      import plot_weight_breakdown 

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import griddata
from plotly.offline import plot


def main():

    # filename = "LHS_Samples"
    # lhs_samples = load_pickle_results(filename)
    # piecewise_plot(lhs_samples,True)

    # filename = 'consolidated_exit_conditions.xlsx'
    # data = load_excel_data(filename)
    # create_parallel_plot(data,True,False)
    
    #case_numbers = [22,16,26]
    #weight_sunburst(case_numbers,True,True)
    #create_single_flight_profiles(case_numbers,True,False)
    #create_life_cycle_analysis(case_numbers,True,True)
    #create_CF_IR(case_numbers,True,False)

    filename = "Cycle Day"
    sensitivity_analysis = load_pickle_results(filename)
    plot_sensitivity_analysis(sensitivity_analysis,filename,True,True)



    return
def plot_sensitivity_analysis(Si,output_variable,save_figure=True, 
                                  show_figure=True,
                                  show_legend=True,
                                  save_filename="sensitivity_analysis_", 
                                  file_type=".png", 
                                  width=6, 
                                  height=6):
    ps = matplotlibstyle()  
    
    # Extract data from Si
    S1 = np.array(Si["S1"])
    ST = np.array(Si["ST"])
    S1_conf = np.array(Si["S1_conf"])
    ST_conf = np.array(Si["ST_conf"])
    S2 = np.array(Si["S2"])
    S2_conf = np.array(Si["S2_conf"])
    
    # Specify custom parameter names
    parameter_names = [
        "HAS\nCapacity\n", 
        "HEX\nCapacity\n", 
        "RES\nCapacity\n"
    ]
    
    N = len(parameter_names)
    
    # Extract upper-triangular S2 values and their confidences
    param_pairs = []
    S2_vals = []
    S2_vals_conf = []
    for i in range(N):
        for j in range(i+1, N):
            param_pairs.append(f"{parameter_names[i]} - {parameter_names[j]}")
            S2_vals.append(S2[i, j])
            S2_vals_conf.append(S2_conf[i, j])
            
    S2_vals = np.array(S2_vals)
    S2_vals_conf = np.array(S2_vals_conf)
    
    # Set up positions and bar widths
    x_s1_st = np.arange(N)
    x_s2 = np.arange(len(param_pairs))
    width_bar = 0.25

    # --- Figure 1: First and Total Order Sensitivities ---
    #fig_ST, ax_ST = plt.subplots(figsize=(width, height-.41)) # Total Weight
    fig_ST, ax_ST = plt.subplots(figsize=(width, height-.3)) # Cycle Day

    
    # Plot ST
    ax_ST.bar(x_s1_st, ST, width_bar, yerr=ST_conf, label='Total',
              color=ps.color[0], alpha=0.9, capsize=4, edgecolor='none', 
              error_kw=dict(ecolor='black', elinewidth=1))
    ax_ST.set_ylim(0, 1)
    # Overlay S1
    ax_ST.bar(x_s1_st, S1, width_bar, yerr=S1_conf, label='First Order',
              color=ps.color[1], alpha=1.0, capsize=4, edgecolor='none',
              error_kw=dict(ecolor='black', elinewidth=1))
    
    ax_ST.set_xticks(x_s1_st)
    ax_ST.set_xticklabels(parameter_names, rotation=0, ha='center')
    ax_ST.set_ylabel('Sensitivity Index', fontsize=ps.axis_font_size)
    #ax_ST.grid(True)
    if show_legend:
        ax_ST.legend(loc='upper left', fontsize=ps.legend_fontsize)
    
    plt.tight_layout()
    
    # Save figure if requested (First and Total)
    if save_figure:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(current_dir, "Final_Plots") 
        subplot_filename = f"{save_dir}/{save_filename}_first_total_{output_variable}{file_type}"
        fig_ST.savefig(subplot_filename, bbox_inches="tight")


    
    
    # --- Figure 2: Second Order Sensitivities ---
    fig_S2, ax_S2 = plt.subplots(figsize=(width, height))
    
    ax_S2.bar(x_s2, S2_vals, width_bar, yerr=S2_vals_conf, label='Second Order',
              color=ps.color[2], alpha=0.8, capsize=4, edgecolor='none',
              error_kw=dict(ecolor='black', elinewidth=1))
    
    ax_S2.set_xticks(x_s2)
    #ax_S2.grid(True)
    #ax_S2.set_ylim(-.04, 0.06)
    ax_S2.set_xticklabels(param_pairs, rotation=0, ha='center')
    ax_S2.set_ylabel('Sensitivity Index', fontsize=ps.axis_font_size)
    
    if show_legend:
        ax_S2.legend(loc='upper left', fontsize=ps.legend_fontsize)
    
    plt.tight_layout()
    
    # Save figure if requested (Second Order)
    if save_figure:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(current_dir, "Final_Plots") 
        subplot_filename = f"{save_dir}/{save_filename}_second_{output_variable}{file_type}"
        fig_S2.savefig(subplot_filename, bbox_inches="tight")
    
    # Show figure if requested
    if show_figure:
        plt.show()

    return

def create_life_cycle_analysis(case_numbers, 
                               save_figure=True, 
                               show_figure=True,
                               save_filename="Life_Cycle_Analysis_", 
                               file_type=".png",
                               width=8, 
                                height=8):
    base_path = "Vehicle_Excel_Files/"
    save_dir = "Final_Plots"
    fig1 = go.Figure()
    fig2 = go.Figure()
    additional_thickness = 1.1

    # Generate evenly spaced positions for the cases
    y_positions = np.linspace(1, len(case_numbers), len(case_numbers))

    # Get styling parameters
    ps = plotlystyle()

    for idx, case_number in enumerate(case_numbers):
        # Load the data
        file_name = base_path + f'case_{case_number}/Raw_Data/e_Twin_Otter_nmc_case_{case_number}.0.xlsx'
        try:
            data = read_excel_file(file_name)  # Custom read function for your Excel files
        except FileNotFoundError:
            print(f"File not found: {file_name}")
            continue

        # Prepare the data
        x = data['Cycle Day']
        z1 = data['Cell Temperature (C)'] - 273  # Convert to Celsius
        z2 =(1-data['Cell SOC (%)'])*100
        y_pos = y_positions[idx]  # Assign evenly spaced y-value for this case

        # Add the main cell temperature curve for each case
        fig1.add_trace(go.Scatter3d(
            x=x,
            y=[y_pos] * len(x),  # Use the evenly spaced y-value
            z=z1,
            mode='lines',  # Ensure smooth and clean lines with markers
            name=f'Configuration {idx+1}',
            line=dict(width=ps['export_scale']*(ps['line_width']+additional_thickness), color=ps['color'][idx % len(ps['color'])]),
            marker=dict(size=ps['marker_size'], symbol=ps['marker_symbols'][idx % len(ps['marker_symbols'])])
        ))
        fig2.add_trace(go.Scatter3d(
            x=x,
            y=[y_pos] * len(x),  # Use the evenly spaced y-value
            z=z2,
            mode='lines',  # Ensure smooth and clean lines with markers
            name=f'Configuration {idx+1}',
            line=dict(width=ps['export_scale']*(ps['line_width']+additional_thickness), color=ps['color'][idx % len(ps['color'])]),
            marker=dict(size=ps['marker_size'], symbol=ps['marker_symbols'][idx % len(ps['marker_symbols'])])
        ))


    # Add a plane at z=50 for fig1 (Thermal Limit)
    fig1.add_trace(
        go.Surface(
            z=50 * np.ones((10, 10)), 
            x=np.linspace(min(x), max(x), 10), 
            y=np.linspace(min(y_positions), max(y_positions), 10), 
            opacity=0.3,  
            colorscale=[[0, "#ffa500"], [1, "#ffa500"]],
            showscale=False,  
            name="Thermal Limit"
        )
    )


    fig2.add_trace(
        go.Surface(
            z=80 * np.ones((10, 10)),  
            x=np.linspace(min(x), max(x), 10),  
            y=np.linspace(min(y_positions), max(y_positions), 10), 
            opacity=0.3, 
            colorscale=[[0, "#ffa500"], [1, "#ffa500"]], 
            showscale=False, 
            name="Capacity Limit"
        )
    )   
    fig1.update_layout(
    scene=dict(
        camera=dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=-0.17053947243903617, y=-0.03337157001708417, z=-0.2811928384252967),
            eye=dict(x=1.5, y=-1.8, z=1.1),
            projection=dict(type="perspective")
        ),
        xaxis=dict(
            title='Days    ',
            titlefont=dict(size=ps["axis_title_font_size"]),
            tickfont=dict(size=ps["tick_font_size"]),
            range=[0, 200]
        ),
        yaxis=dict(
            title='',  # Remove y-axis title
            showticklabels=False,
            showgrid=False
        ),
        zaxis=dict(
            title='Cell Temperature (°C)',
            titlefont=dict(size=ps["axis_title_font_size"]),
            tickfont=dict(size=ps["tick_font_size"]),
            range=[15, 60]
        )
    ),
    legend=dict(
        font=dict(size=ps["legend_font_size"]),
        title=dict(font=dict(size=ps["legend_title_font_size"]))
    ),
    template='plotly_white',
    plot_bgcolor=ps["plot_bgcolor"],
    paper_bgcolor=ps["paper_bgcolor"],
    showlegend=True  # Add legend for clarity
)

    fig2.update_layout(
    scene=dict(
        camera=dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=-0.17053947243903617, y=-0.03337157001708417, z=-0.2811928384252967),
            eye=dict(x=1.5, y=-1.8, z=1.1),
            projection=dict(type="perspective")
        ),
        xaxis=dict(
            title='Days     ',
            titlefont=dict(size=ps["axis_title_font_size"]),
            tickfont=dict(size=ps["tick_font_size"]),
            range=[0, 200]
        ),
        yaxis=dict(
            title='',  # Remove y-axis title
            showticklabels=False,
            showgrid=False
        ),
        zaxis=dict(
            title='Depth of Discharge (%)',
            titlefont=dict(size=ps["axis_title_font_size"]),
            tickfont=dict(size=ps["tick_font_size"]),
            range=[0.1, 100]
        )
    ),
    legend=dict(
        font=dict(size=ps["legend_font_size"]),
        title=dict(font=dict(size=ps["legend_title_font_size"]))
    ),
    template='plotly_white',
    plot_bgcolor=ps["plot_bgcolor"],
    paper_bgcolor=ps["paper_bgcolor"],
    showlegend=True  # Add legend for clarity
)
    
    fig1.add_annotation(x=0.55, y=0.8,
            text="Thermal Limit",
            showarrow=True,
            arrowhead=2,
            ax=0,  
            ay=-50,
            font=dict(size=ps["axis_title_font_size"], color="black"),)#ADD FONT SIZE
    fig2.add_annotation(x=0.55, y=0.8,
            text="Capacity Limit",
            showarrow=True,
            arrowhead=2,
            ax=0,  
            ay=-50,  
            font=dict(size=ps["axis_title_font_size"], color="black"),)

    
    # Save or show the figure
    if save_figure:
         # Determine the current directory and create a folder for saving plots
        current_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(current_dir, "Final_Plots")
        fig1.write_image(os.path.join(save_dir, f"{save_filename}_Temperature{file_type}"),width=width*100, height=height*100, scale=ps["export_scale"])
        fig2.write_image(os.path.join(save_dir, f"{save_filename}_State_of_charge{file_type}"),width=width*100, height=height*100, scale=ps["export_scale"])
       
    if show_figure:
        fig1.show()
        fig2.show()

    return

def create_single_flight_profiles(case_numbers, 
                                  save_figure=True, 
                                  show_figure=True,
                                  show_legend=True,
                                  save_filename="single_flight_", 
                                  file_type=".png", 
                                  width=6, 
                                  height=6):
    # Get plot style parameters
    ps = matplotlibstyle()

    base_path = "Vehicle_Excel_Files/"
    save_dir = "Final_Plots"

    # Data containers for subplots
    fig_data = [([], [], "Cell SOC (%)", "Time (min)"), 
                ([], [], "Cell Temperature (°C)", "Time (min)"), 
                ([], [], "Heat Generation (W)", "Time (min)")]

    # Populate data for each subplot
    for case_number in case_numbers:
        file_name = base_path + f'case_{case_number}/Raw_Data/e_Twin_Otter_nmc_case_{case_number}.0.xlsx'
        data = load_excel_data(file_name, sheet_name="Group_1")
        transition_index = data[(data["Cell Current (A)"] == 0) & (data["Cell Current (A)"].shift(-1) != 0)].index + 1
        limited_df = data.iloc[:transition_index[0]]

        # Append data for each plot
        fig_data[0][0].append(limited_df["Time (min)"])
        fig_data[0][1].append(limited_df["Cell SOC (%)"] * 100)
        fig_data[1][0].append(limited_df["Time (min)"])
        fig_data[1][1].append(limited_df["Cell Temperature (C)"] - 273)
        fig_data[2][0].append(limited_df["Time (min)"])
        fig_data[2][1].append(limited_df["Module Heat Generation (W)"] / 12)

    # Create individual figures
    for i, (x_data, y_data, ylabel, xlabel) in enumerate(fig_data):
        fig, ax = plt.subplots(figsize=(width, height))

        # Plot all configurations
        for j, (x, y) in enumerate(zip(x_data, y_data)):
            ax.plot(
                x, y,
                label=f"Config {j+1}",
                color=ps.line_color[j],
                linewidth=ps.line_width,
                linestyle=ps.line_style[j % len(ps.line_style)],
                marker=ps.markers[j % len(ps.markers)],
                markersize=ps.marker_size
            )

        ax.set_xlabel(xlabel, fontsize=ps.axis_font_size)
        ax.set_ylabel(ylabel, fontsize=ps.axis_font_size)
        ax.grid(True)

        # Add zoomed-in rectangle and inset for SOC subplot
        if i == 0:  # SOC plot
            rect = patches.Rectangle((30, 40), 10, 10, linewidth=2, edgecolor='black', facecolor='yellow', alpha=0.3)
            ax.add_patch(rect)

            # Add inset
            inset_ax = fig.add_axes([0.6, 0.25, 0.25, 0.25])  # [x, y, width, height]
            for j in range(len(x_data)):
                # Access x and y data for the current case
                x = x_data[j]
                y = y_data[j]
                
                # Filter the x and y values for the zoomed range
                zoomed_x = x[(x >= 30) & (x <= 40)]
                zoomed_y = y[zoomed_x.index[0] : zoomed_x.index[-1]+1]
                
                # Plot the zoomed data
                inset_ax.plot(
                    zoomed_x, 
                    zoomed_y, 
                    label=f"Config {case_numbers[j]}", 
                    color=ps.line_color[j],
                    linewidth=ps.line_width, 
                    linestyle=ps.line_style[j % len(ps.line_style)],
                    marker=ps.markers[j % len(ps.markers)],
                    markersize=ps.marker_size
                )

        # Add legend only to Heat Generation subplot
        if i == 2:  # Heat Generation plot
            ax.legend(fontsize=ps.legend_fontsize)

        # Save figure
        if save_figure:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            save_dir = os.path.join(current_dir, "Final_Plots") 
            subplot_filename = f"{save_dir}/{save_filename}_subplot_{i + 1}{file_type}"
            fig.savefig(subplot_filename, bbox_inches="tight")
            print(f"Saved subplot {i + 1} to: {subplot_filename}")

    return

def weight_sunburst(case_numbers, 
                    save_figure=False, 
                    show_figure=True,
                    title="Sunburst Distribution",
                    save_filename="weight_distribution", 
                    file_type=".png", 
                    width=6, 
                    height=6):
    # Get plot style parameters
    ps = plotlystyle()

    base_path = "Vehicle_Excel_Files/"
    # Look for vehicle file and load the result
    for case_number in case_numbers:
        save_filename_1 = save_filename + '_' + str(case_number)
        file_name = base_path + 'case_' + str(case_number) + "/Raw_Data/e_twin_otter_vehicle"
        vehicle = load_pickle_results(file_name)
        weight_analysis = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
        weight_analysis.vehicle = vehicle
        weight_analysis.settings.miscelleneous_weight_factor = 1.0
        weight = weight_analysis.evaluate()
        print(weight)
        sunburst_plot = plot_weight_breakdown(
            weight_analysis.vehicle,
            show_figure=False,
            show_legend=True,
            SI_Units=True,
            aircraft_name=case_number,
            file_type=".png",
            width=width,
            height=height
        )

        # Apply styling to the plot
        sunburst_plot.update_layout(
            font=dict(family=ps["font_family"]),
        )

        if show_figure:
            sunburst_plot.show()

        # Save the figure if requested
        if save_figure:
            # Determine the current directory and create a folder for saving plots
            current_dir = os.path.dirname(os.path.abspath(__file__))
            save_dir = os.path.join(current_dir, "Final_Plots")
            save_path = os.path.join(save_dir, f"{save_filename_1}{file_type}")
            # Save the plot using Plotly's write_image method
            sunburst_plot.write_image(save_path, width=width*100, height=height*100, scale=ps["export_scale"])
            print(f"Figure saved at: {save_path}")

    return

def create_parallel_plot(dataframe, 
                         save_figure=False, 
                         show_figure = True,
                         title="Multivariable Parallel Plot",
                         save_filename="parallel_plot", 
                         file_type=".png", 
                         width=12, 
                         height=8):
    
    ps = plotlystyle()
    # Columns to include in the plot
    columns = [
        "Case ID",
        "HAS_power",
        "HEX_power",
        "RES_dimensions",
        "TMS  Weight",
        "Total Weight",
        "Cycle Day"
    ]
    # Desired axis labels
    axis_labels = {
        "Case ID": "Configuration",
        "HAS_power": "Heat Acquisition Power",
        "HEX_power": "Heat Exchanger Power",
        "RES_dimensions": "Reservoir Dimension",
        "TMS  Weight": "TMS Weight",
        "Total Weight": "Aircraft Weight",
        "Cycle Day": "Life of Battery Pack"
    }

    # Create the parallel coordinates plot
    fig = px.parallel_coordinates(
        dataframe,
        dimensions=columns,
        color=dataframe["Cycle Day"],  # Color by the last column
        color_continuous_scale=ps["color_scale"],
        labels=axis_labels 
    )


    # Update layout for custom styling
    fig.update_layout(
        font=dict(family=ps["font_family"], size=ps["font_size"]),
        coloraxis_colorbar=dict(
            title=dict(text=columns[-1], font=dict(size=ps["axis_title_font_size"]))
        ),
        height=600,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    if show_figure:
        fig.show()

    # Save the figure if requested
    if save_figure:
        # Determine the current directory and create a folder for saving plots
        current_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(current_dir, "Final_Plots")
        save_path = os.path.join(save_dir, f"{save_filename}{file_type}")
        # Save the plot using Plotly's write_image method
        fig.write_image(save_path, width=width*100, height=height*100)
        print(f"Figure saved at: {save_path}")

    return fig

def piecewise_plot(lhs_samples, save_figure=False, save_filename="LHS_Sample_Pairwise_Plot", 
                   file_type=".png", width=12, height=8):
    lhs_df = pd.DataFrame(lhs_samples, columns=["Heat Acqusition Power (W)", "Heat Exchanger Power (W)", "Reservoir Size (m)"])
    
    # Apply the existing plot style
    ps = matplotlibstyle()
    
    # Set theme and create the pairwise plot with enhancements
    sns.set_theme(style="whitegrid", font_scale=1.2)
    pairplot = sns.pairplot(
        lhs_df, 
        diag_kind="kde", 
        plot_kws={'alpha': 0.8, 's': 70},  # Transparency and marker size
        diag_kws={'shade': True},  # Shaded KDE
    )

    # Adjust title placement and font size
    #pairplot.figure.subplots_adjust(top=0.9)  # Adjust the top margin to fit the title
    #pairplot.figure.suptitle("Matrix of Pairwise Plots for LHS Samples", 
                          #y=0.95, fontsize=ps.title_font_size)

    # Set axis labels style
    for ax in pairplot.axes.flatten():
        if ax is not None:
            ax.set_xlabel(ax.get_xlabel(), fontsize=ps.axis_font_size)
            ax.set_ylabel(ax.get_ylabel(), fontsize=ps.axis_font_size)

    # Customize figure size
    pairplot.figure.set_size_inches(width, height)

    # Save the figure if required
    if save_figure:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        save_dir = os.path.join(current_dir, "Final_Plots")
        save_path = os.path.join(save_dir, f"{save_filename}{file_type}")
        pairplot.savefig(save_path, bbox_inches='tight')
        print(f"Figure saved at: {save_path}")

    return



# ****************************************************************************************************************************************
# ****************************************************************************************************************************************
# ------------------------------------------------------------------
#  Plot Parameters for Matplotlib
# ------------------------------------------------------------------   
def matplotlibstyle():
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 18-4,
                  'ytick.labelsize': 18,
                  'axes.titlesize': 18,
                  'figure.dpi': 128
                  }


    # Universal Plot Settings  
    plt.rcParams.update(parameters)
    plot_parameters                        = Data()
    plot_parameters.line_width             = 1  
    plot_parameters.line_style             = ['-','-']
    plot_parameters.marker_size            = 4
    plot_parameters.legend_fontsize        = '18'
    plot_parameters.legend_title_font_size = 18
    plot_parameters.axis_font_size         = 14+4
    plot_parameters.title_font_size        = 18   
    plot_parameters.markers                =  ['o','x','v','P','p','^','D','*']
    plot_parameters.color                  = [ '#003f5c',
                                            '#ffa600',
                                            '#ff6361',]
    plot_parameters.line_color = [
    '#440154',  # Dark Purple
    '#31688E',  # Blue
    '#1F9E89',  # Green-Blue
   
    '#35B779',  # Green
    '#6DCD59',  # Yellow-Green
    '#FDE725'   # Bright Yellow
    ]

    return plot_parameters
# ------------------------------------------------------------------
#  Plot Parameters for Plotly
# ------------------------------------------------------------------  
def plotlystyle():
    style_params = {
        "font_family": "Times New Roman",
        "font_size": 18,
        "line_width": 1.5,
        "line_dash": ["solid", "dot"],
        "marker_size": 8,
        "marker_symbols": ["circle", "x", "square", "triangle-up", "star"],
        "color_scale": "Inferno",
        "legend_font_size": 14,
        "legend_title_font_size": 16,
        "axis_title_font_size": 18,
        "title_font_size": 20,
        "tick_font_size": 14,  
        "plot_bgcolor": "white",  
        "paper_bgcolor": "white",  
        "grid_color": "black",  
        "grid_width": 0.5  ,
       "color" : ['#440154',  # Dark Purple
    '#31688E',  # Blue
    '#1F9E89',  # Green-Blue
                "#818bf8", 
                "#49dab4", 
                "red", 
                "#FFB7A7",  # Peach
                "#C9ACE6",  # Lavender
                "#FFD7E9"   # Light Pink
            ],
        "export_scale": 5,  # Increase scale for higher resolution
        "line_width": 1.5
    }
    
    return style_params

# ------------------------------------------------------------------
#   Load Picke Results
# ------------------------------------------------------------------   
def load_pickle_results(filename):  
    current_dir = os.path.dirname(os.path.abspath(__file__))
    load_dir = os.path.join(current_dir)
    load_file = os.path.join(load_dir, filename + '.pkl')
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results

# ------------------------------------------------------------------
#   Load All Groups Excel Results
# ------------------------------------------------------------------  
def read_excel_file(filename):
    # Read the Excel file
    file_path = os.path.join(os.path.dirname(__file__), filename)
    excel_file = pd.ExcelFile(file_path)
    sheet_names = excel_file.sheet_names
    print("\nAvailable sheets in the Excel file:",filename)
    
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
# ------------------------------------------------------------------
#   Load Excel Results
# ------------------------------------------------------------------   
def load_excel_data(file_name, sheet_name="Sheet1"):
      file_path = os.path.join(os.path.dirname(__file__), file_name)
      return pd.read_excel(file_path, sheet_name=sheet_name)



if __name__ == '__main__':
    main()
    plt.show()

