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
    file_type   =  'png'
    save_figure =  True
    width       =  5
    height      =  5
    generate_saf_plots(file_type,save_figure,width, height)
    generate_electrification_plots(file_type,save_figure,width, height)
    generate_hydrogen_plots(file_type,save_figure,width, height)
    return



def generate_saf_plots(file_type,save_figure,width, height): 
    saf_data =  load_results('saf_data.res')

    return  


def generate_electrification_plots(file_type,save_figure,width, height):  
    electric_data =  load_results('electric_data.res')
    
    
    # get plotting style 
    ps      = plot_style()  

    parameters = {'axes.labelsize': ps.axis_font_size,
                  'xtick.labelsize': ps.axis_font_size,
                  'ytick.labelsize': ps.axis_font_size,
                  'axes.titlesize': ps.title_font_size}
    plt.rcParams.update(parameters)
     
    # get line colors for plots 
    line_colors   = cm.inferno(np.linspace(0,0.9,len(results.segments)))     
     
    fig_1   = plt.figure(save_filename+"_True_Airspeed")
    fig_2   = plt.figure(save_filename+"_Equiv._Airspeed")
    fig_3   = plt.figure(save_filename+"_Calibrated_Airspeed")
    fig_4   = plt.figure(save_filename+"_Mach_Number")
    
    fig_1.set_size_inches(width,height)
    fig_2.set_size_inches(width,height)
    fig_3.set_size_inches(width,height)
    fig_4.set_size_inches(width,height)

    axis_1 = fig_1.add_subplot(1,1,1) 
    axis_2 = fig_2.add_subplot(1,1,1)    
    axis_3 = fig_3.add_subplot(1,1,1) 
      
      
    axis_1.plot(time, velocity, color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = segment_name)
    axis_1.set_ylabel(r'True Airspeed (kts)')
    axis_1.set_xlabel('Time (mins)')        
    set_axes(axis_1)    
    
    axis_2.plot(time, EAS, color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = segment_name) 
    axis_2.set_ylabel(r'Equiv. Airspeed (kts)')
    axis_2.set_xlabel('Time (mins)')
    set_axes(axis_2) 

    axis_3.plot(time, CAS, color = line_colors[i], marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = segment_name)
    axis_3.set_xlabel('Time (mins)')
    axis_3.set_ylabel(r'Calibrated Airspeed (kts)')
    set_axes(axis_3) 
         
    
    if show_legend:    
        leg1 =  fig_1.legend(bbox_to_anchor=(0.5, 1.0), loc='upper center')
        leg2 =  fig_2.legend(bbox_to_anchor=(0.5, 1.0), loc='upper center')
        leg3 =  fig_3.legend(bbox_to_anchor=(0.5, 1.0), loc='upper center') 
        
    # Adjusting the sub-plots for legend 
    fig_1.tight_layout()    
    fig_2.tight_layout()    
    fig_3.tight_layout()    
    fig_4.tight_layout()
     
    
    if save_figure:
        fig_1.savefig(save_filename + file_type)
        fig_2.savefig(save_filename + file_type)
        fig_3.savefig(save_filename + file_type) 
    return    


def generate_hydrogen_plots(file_type,save_figure,width, height): 
    hydrogen_data =  load_results('hydrogen_data.res')    

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
