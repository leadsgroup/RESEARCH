import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import os




def main():
    file_name = 'Optimization_Results_Analysis/consolidated_exit_conditions.xlsx'
    file_path = os.path.join(os.path.dirname(__file__), file_name)
    sheet_name = 'Sheet1'
    data = load_data(file_path, sheet_name)
    plot1 = create_parallel_plot(data)
    plot1.show()
    return 

# Define plot styling parameters
def plot_style():
    style_params = {
        "font_family": "Times New Roman",
        "font_size": 16,
        "line_width": 5,
        "line_dash": ["solid", "dot"],
        "marker_size": 8,
        "marker_symbols": ["circle", "x", "square", "triangle-up", "star"],
        "color_scale": "Inferno", 
        "legend_font_size": 14,
        "legend_title_font_size": 16,
        "axis_title_font_size": 18,
        "title_font_size": 20,
    }
    return style_params


def create_parallel_plot(dataframe, title="Multivariable Parallel Plot"):
    """
    Creates a parallel coordinates plot using Plotly.
    
    Parameters:
        dataframe (pd.DataFrame): The dataframe containing the data.
        columns (list): The list of columns to include in the plot.
        title (str): The title of the plot.
    
    Returns:
        fig (plotly.graph_objs._figure.Figure): The generated Plotly figure.
    """
    columns = [
    "Case ID",
    "HAS_power",
    "HEX_power",
    "RES_dimensions",
    "TMS Weight",
    "Total Weight",
    "Cycle Day"
    ]

     # Get plot style parameters
    ps = plot_style()

    # Filter the dataframe to include only the specified columns
    filtered_df = dataframe[columns]

    # Create the parallel coordinates plot
    fig = px.parallel_coordinates(
        filtered_df,
        dimensions=columns,
        color=columns[-1],
        color_continuous_scale=ps["color_scale"],
        title=title
    )

    # Update layout for custom styling
    fig.update_layout(
        title=dict(text=title, font=dict(size=ps["title_font_size"], family=ps["font_family"])),
        font=dict(family=ps["font_family"], size=ps["font_size"]),
        coloraxis_colorbar=dict(
            title=dict(text=columns[-1], font=dict(size=ps["axis_title_font_size"]))
        ),
        height=600,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    return fig
    
# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_data(file_path, sheet_name):
      return pd.read_excel(file_path, sheet_name=sheet_name)



if __name__ == '__main__': 
    main()    
    
   
