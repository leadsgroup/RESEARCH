'''





'''

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTS 
# --------------------------------------------------------------------------------------------------------------------------------------------------- 

from RCAIDE.Framework.Core import  Units 
from RCAIDE.load import load as load_results
from RCAIDE.save import save as save_results 

import os 
import numpy as np           
import pandas as pd 
import json
from pyatmos import coesa76

from urllib.request import urlopen

import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import plotly.io as pio 
import plotly.express as px 
from plotly.graph_objs import * 
 
def main(): 
    
    generate_saf_plots()
    generate_electrification_plots()
    generate_electrification_plots()
     
    return



def generate_saf_plots(): 
          
        return  


def generate_electrification_plots(): 
          
        return  

     
def generate_hydrogen_plots(): 
          
        return  

 
    
# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS COMMON AXIS 
# ----------------------------------------------------------------------------------------------------------------------    
## @ingroup Visualization-Performance-Common
def set_axes(axes):
    """This sets the axis parameters for all plots

    Assumptions:
    None

    Source:
    None

    Inputs
    axes

    Outputs:
    axes

    Properties Used:
    N/A 
    """

    axes.minorticks_on()
    axes.grid(which='major', linestyle='-', linewidth=0.5, color='grey')
    axes.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')
    axes.grid(True)
    axes.get_yaxis().get_major_formatter().set_scientific(False)
    axes.get_yaxis().get_major_formatter().set_useOffset(False)

    return
 
def plot_style(): 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18,
                  'figure.dpi': 1200
                  }


    # Universal Plot Settings  
    plt.rcParams.update(parameters)
    plot_parameters                        = Data()
    plot_parameters.line_width             = 1.5  
    plot_parameters.line_style             = ['-','-']
    plot_parameters.marker_size            = 4
    plot_parameters.legend_fontsize        = '12'
    plot_parameters.legend_title_font_size = 14
    plot_parameters.axis_font_size         = 16
    plot_parameters.title_font_size        = 16   
    plot_parameters.markers                =  ['o','x','o','v','P','p','^','D','*']
    plot_parameters.color                  = 'black' 

    return plot_parameters
 
     
if __name__ == '__main__':
    main()
