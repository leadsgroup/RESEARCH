'''
# Simulation_Repeated_Flight_Operations.py
#
# Created: May 2024, S S. Shekar

'''

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units, Data   
import pickle
from RCAIDE.Library.Plots                                           import *  

import time  
import numpy as np
import pylab as plt
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.transforms import Bbox
import pandas as pd
import sys
import os
import matplotlib
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource

def main():
    
    file_name_1 = '../e_Twin_Otter_ORD_January'  # -5.1 Mean temp  
    ord_jan     = load_results(file_name_1) 

    file_name_2 = '../e_Twin_Otter_ORD_December' # -2.25 Mean Temp 
    ord_dec = load_results(file_name_2)
    
    file_name_3 = '../e_Twin_Otter_ORD_March' # 3 Mean Temp 
    ord_mar = load_results(file_name_3)
    
    file_name_4 = '../e_Twin_Otter_ORD_November' # 4 Mean Temp 
    ord_nov = load_results(file_name_4)
    
    file_name_5 = '../e_Twin_Otter_DFW_January' # 7.1 Mean Temp 
    dfw_jan = load_results(file_name_5)            
    
    file_name_6 = '../e_Twin_Otter_DFW_December' # 8.5 Mean Temp 
    dfw_dec = load_results(file_name_6)
    
    file_name_7 = '../e_Twin_Otter_DFW_February' # 10.5 Mean Temp 
    dfw_feb = load_results(file_name_7)
    
    file_name_8 = '../e_Twin_Otter_ORD_October' # 11.5 Mean Temp 
    ord_oct = load_results(file_name_8)    
    
    file_name_9 = '../e_Twin_Otter_DFW_November' # 14 Mean Temp 
    dfw_nov = load_results(file_name_9)        
    
    file_name_10 = '../e_Twin_Otter_ORD_September' # 18 Mean Temp 
    ord_sep = load_results(file_name_10)

    file_name_11 = '../e_Twin_Otter_DFW_October' # 20 Mean Temp 
    dfw_oct = load_results(file_name_11)                  
    
    file_name_12 = '../e_Twin_Otter_ORD_June' # 21 Mean Temp 
    ord_jun = load_results(file_name_12)
    
    file_name_13 = '../e_Twin_Otter_ORD_August' # 22.5 Mean Temp 
    ord_aug = load_results(file_name_13)
    
    file_name_14 = '../e_Twin_Otter_DFW_May' # 23.5 Mean Temp 
    dfw_may = load_results(file_name_14)
    
    file_name_15 = '../e_Twin_Otter_ORD_July' # 24 Mean Temp 
    ord_jul = load_results(file_name_15)        
    
    file_name_16 = '../e_Twin_Otter_DFW_September' # 25.5 Mean Temp 
    dfw_sep = load_results(file_name_16)
    
    file_name_17 = '../e_Twin_Otter_DFW_June' # 28.5 Mean Temp 
    dfw_jun = load_results(file_name_17)            
                
    file_name_18 = '../e_Twin_Otter_DFW_July' # 30.75 Mean Temp 
    dfw_jul = load_results(file_name_18)        
   
   
    no_of_fligths = 1  

    
    #city_plots()
    #multi_city_plots()
    #monthly_weather()
    #plot_flight_profile(dfw_oct, 6, 3)
    #plot_battery_cell_conditions_1(dfw_oct, 6, 3)
    four_dimensional_contour_plot(ord_jan, ord_dec, ord_mar, ord_nov, dfw_jan, dfw_dec,dfw_feb, ord_oct, dfw_nov, ord_sep, dfw_oct, ord_jun, ord_aug, dfw_may, ord_jul, dfw_sep, dfw_jun, dfw_jul, no_of_fligths)
    #batter_temperature_contour_plot(ord_jan, ord_dec, ord_mar, ord_nov, dfw_jan, dfw_dec,dfw_feb, ord_oct, dfw_nov, ord_sep, dfw_oct, ord_jun, ord_aug, dfw_may, ord_jul, dfw_sep, dfw_jun, dfw_jul, no_of_fligths)
    #has_operation_contour_plot(ord_jan, ord_dec, ord_mar, ord_nov, dfw_jan, dfw_dec,dfw_feb, ord_oct, dfw_nov, ord_sep, dfw_oct, ord_jun, ord_aug, dfw_may, ord_jul, dfw_sep, dfw_jun, dfw_jul, no_of_fligths)
    #reservoir_temperature_contour_plot (ord_jan, ord_dec, ord_mar, ord_nov, dfw_jan, dfw_dec,dfw_feb, ord_oct, dfw_nov, ord_sep, dfw_oct, ord_jun, ord_aug, dfw_may, ord_jul, dfw_sep, dfw_jun, dfw_jul, no_of_fligths)    
    
    return

# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS STYLE
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Visualization-Performance-Common
def plot_style(): 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18,
                  'figure.dpi': 200
                  }


    # Universal Plot Settings  
    plt.rcParams.update(parameters)
    plot_parameters                        = Data()
    plot_parameters.line_width             = 1  
    plot_parameters.line_style             = ['-','-']
    plot_parameters.marker_size            = 4
    plot_parameters.legend_fontsize        = '18'
    plot_parameters.legend_title_font_size = 14
    plot_parameters.axis_font_size         = 14# + 2
    plot_parameters.title_font_size        = 16   
    plot_parameters.markers                =  ['o','x','o','v','P','p','^','D','*']
    plot_parameters.color                  = 'red'


    return plot_parameters

# ------------------------------------------------------------------
#   Create Contour Plots for TMS Operation
# ------------------------------------------------------------------
def has_operation_contour_plot (ord_jan, ord_dec, ord_mar, ord_nov, dfw_jan, dfw_dec,dfw_feb, ord_oct,
                                     dfw_nov, ord_sep, dfw_oct, ord_jun, ord_aug, dfw_may, ord_jul, dfw_sep, dfw_jun, dfw_jul, no_of_fligths, 
                   save_figure = True,
                   save_filename = "Thermal_Management_System_Operating_conditions_multiple_flights_plot",
                   show_legend = False, 
                   file_type = ".png",
                   width = 12, height = 7):
    
  
    ord_jan_time = []
    ord_dec_time = []
    ord_mar_time = []
    ord_nov_time = []
    dfw_jan_time = []
    dfw_dec_time = []
    dfw_feb_time = []
    ord_oct_time = []
    dfw_nov_time = []
    ord_sep_time = []
    dfw_oct_time = []
    ord_jun_time = []
    ord_aug_time = []
    dfw_may_time = []
    ord_jul_time = []
    dfw_sep_time = []
    dfw_jun_time = []
    dfw_jul_time = []
    
   
    ord_jan_percent_operation_HAS = []
    ord_dec_percent_operation_HAS = []
    ord_mar_percent_operation_HAS = []
    ord_nov_percent_operation_HAS = []
    dfw_jan_percent_operation_HAS = []
    dfw_dec_percent_operation_HAS = []
    dfw_feb_percent_operation_HAS = []
    ord_oct_percent_operation_HAS = []
    dfw_nov_percent_operation_HAS = []
    ord_sep_percent_operation_HAS = []
    dfw_oct_percent_operation_HAS = []
    ord_jun_percent_operation_HAS = []
    ord_aug_percent_operation_HAS = []
    dfw_may_percent_operation_HAS = []
    ord_jul_percent_operation_HAS = []
    dfw_sep_percent_operation_HAS = []
    dfw_jun_percent_operation_HAS = []
    dfw_jul_percent_operation_HAS = []
    
    ord_jan_percent_operation_HEX = []
    ord_dec_percent_operation_HEX = []
    ord_mar_percent_operation_HEX = []
    ord_nov_percent_operation_HEX = []
    dfw_jan_percent_operation_HEX = []
    dfw_dec_percent_operation_HEX = []
    dfw_feb_percent_operation_HEX = []
    ord_oct_percent_operation_HEX = []
    dfw_nov_percent_operation_HEX = []
    ord_sep_percent_operation_HEX = []
    dfw_oct_percent_operation_HEX = []
    ord_jun_percent_operation_HEX = []
    ord_aug_percent_operation_HEX = []
    dfw_may_percent_operation_HEX = []
    ord_jul_percent_operation_HEX = []
    dfw_sep_percent_operation_HEX = []
    dfw_jun_percent_operation_HEX = []
    dfw_jul_percent_operation_HEX = []
        
    
    mean_temperature_ord_jan = ord_jan.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_dec = ord_dec.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_mar = ord_mar.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_nov = ord_nov.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jan = dfw_jan.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_dec = dfw_dec.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_feb = dfw_feb.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_oct = ord_oct.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_nov = dfw_nov.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_sep = ord_sep.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_oct = dfw_oct.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_jun = ord_jun.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_aug = ord_aug.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_may = dfw_may.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_jul = ord_jul.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_sep = dfw_sep.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jun = dfw_jun.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jul = dfw_jul.segments[0].conditions.freestream.temperature[0] - 273
    
    mean_temperature = np.array((mean_temperature_ord_jan, mean_temperature_ord_dec, mean_temperature_ord_mar,
                                 mean_temperature_ord_nov,  mean_temperature_dfw_jan,mean_temperature_dfw_dec, mean_temperature_dfw_feb ,
                                 mean_temperature_ord_oct, mean_temperature_dfw_nov, mean_temperature_ord_sep , mean_temperature_dfw_oct,
                                 mean_temperature_ord_jun, mean_temperature_ord_aug  , mean_temperature_dfw_may, mean_temperature_ord_jul ,
                                 mean_temperature_dfw_sep, mean_temperature_dfw_jun , mean_temperature_dfw_jul))      

    for network in ord_jan.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                for i in range(14*no_of_fligths): 
                    time_ord_jan = ord_jan.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_dec = ord_dec.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_mar = ord_mar.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_nov = ord_nov.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jan = dfw_jan.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_dec = dfw_dec.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_feb = dfw_feb.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_oct = ord_oct.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_nov = dfw_nov.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_sep = ord_sep.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_oct = dfw_oct.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_jun = ord_jun.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_aug = ord_aug.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_may = dfw_may.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_jul = ord_jul.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_sep = dfw_sep.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jun = dfw_jun.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jul = dfw_jul.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    
                    ord_jan_time.append(time_ord_jan)                               
                    ord_dec_time.append(time_ord_dec)                  
                    ord_mar_time.append(time_ord_mar)                    
                    ord_nov_time.append(time_ord_nov)
                    dfw_jan_time.append(time_dfw_jan)                                       
                    dfw_dec_time.append(time_dfw_dec)                                       
                    dfw_feb_time.append(time_dfw_feb)                                        
                    ord_oct_time.append(time_ord_oct)                                       
                    dfw_nov_time.append(time_dfw_nov)                                       
                    ord_sep_time.append(time_ord_sep)                                      
                    dfw_oct_time.append(time_dfw_oct)                 
                    ord_jun_time.append(time_ord_jun)
                    ord_aug_time.append(time_ord_aug)
                    dfw_may_time.append(time_dfw_may) 
                    ord_jul_time.append(time_ord_jul)                                       
                    dfw_sep_time.append(time_dfw_sep)                                                                                                                                                                                                                                                                     
                    dfw_jun_time.append(time_dfw_jun)
                    dfw_jul_time.append(time_dfw_jul)
                    
                    percent_operation_HAS_ord_jan = ord_jan.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_ord_dec = ord_dec.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_ord_mar = ord_mar.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_ord_nov = ord_nov.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_jan = dfw_jan.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_dec = dfw_dec.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_feb = dfw_feb.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_ord_oct = ord_oct.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_nov = dfw_nov.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_ord_sep = ord_sep.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_oct = dfw_oct.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_ord_jun = ord_jun.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_ord_aug = ord_aug.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_may = dfw_may.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_ord_jul = ord_jul.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_sep = dfw_sep.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_jun = dfw_jun.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    percent_operation_HAS_dfw_jul = dfw_jul.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HAS.percent_operation[:,0] * 100
                    
                    ord_jan_percent_operation_HAS.append(percent_operation_HAS_ord_jan)                               
                    ord_dec_percent_operation_HAS.append(percent_operation_HAS_ord_dec)                  
                    ord_mar_percent_operation_HAS.append(percent_operation_HAS_ord_mar)                    
                    ord_nov_percent_operation_HAS.append(percent_operation_HAS_ord_nov)
                    dfw_jan_percent_operation_HAS.append(percent_operation_HAS_dfw_jan)                                       
                    dfw_dec_percent_operation_HAS.append(percent_operation_HAS_dfw_dec)                                       
                    dfw_feb_percent_operation_HAS.append(percent_operation_HAS_dfw_feb)                                        
                    ord_oct_percent_operation_HAS.append(percent_operation_HAS_ord_oct)                                       
                    dfw_nov_percent_operation_HAS.append(percent_operation_HAS_dfw_nov)                                       
                    ord_sep_percent_operation_HAS.append(percent_operation_HAS_ord_sep)                                      
                    dfw_oct_percent_operation_HAS.append(percent_operation_HAS_dfw_oct)                 
                    ord_jun_percent_operation_HAS.append(percent_operation_HAS_ord_jun)
                    ord_aug_percent_operation_HAS.append(percent_operation_HAS_ord_aug)
                    dfw_may_percent_operation_HAS.append(percent_operation_HAS_dfw_may) 
                    ord_jul_percent_operation_HAS.append(percent_operation_HAS_ord_jul)                                       
                    dfw_sep_percent_operation_HAS.append(percent_operation_HAS_dfw_sep)                                                                                                                                                                                                                                                                     
                    dfw_jun_percent_operation_HAS.append(percent_operation_HAS_dfw_jun)
                    dfw_jul_percent_operation_HAS.append(percent_operation_HAS_dfw_jul)
                    
                    percent_operation_HEX_ord_jan = ord_jan.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_ord_dec = ord_dec.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_ord_mar = ord_mar.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_ord_nov = ord_nov.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_jan = dfw_jan.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_dec = dfw_dec.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_feb = dfw_feb.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_ord_oct = ord_oct.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_nov = dfw_nov.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_ord_sep = ord_sep.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_oct = dfw_oct.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_ord_jun = ord_jun.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_ord_aug = ord_aug.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_may = dfw_may.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_ord_jul = ord_jul.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_sep = dfw_sep.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_jun = dfw_jun.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    percent_operation_HEX_dfw_jul = dfw_jul.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.HEX.percent_operation[:,0] * 100
                    
                    ord_jan_percent_operation_HEX.append(percent_operation_HEX_ord_jan)                               
                    ord_dec_percent_operation_HEX.append(percent_operation_HEX_ord_dec)                  
                    ord_mar_percent_operation_HEX.append(percent_operation_HEX_ord_mar)                    
                    ord_nov_percent_operation_HEX.append(percent_operation_HEX_ord_nov)
                    dfw_jan_percent_operation_HEX.append(percent_operation_HEX_dfw_jan)                                       
                    dfw_dec_percent_operation_HEX.append(percent_operation_HEX_dfw_dec)                                       
                    dfw_feb_percent_operation_HEX.append(percent_operation_HEX_dfw_feb)                                        
                    ord_oct_percent_operation_HEX.append(percent_operation_HEX_ord_oct)                                       
                    dfw_nov_percent_operation_HEX.append(percent_operation_HEX_dfw_nov)                                       
                    ord_sep_percent_operation_HEX.append(percent_operation_HEX_ord_sep)                                      
                    dfw_oct_percent_operation_HEX.append(percent_operation_HEX_dfw_oct)                 
                    ord_jun_percent_operation_HEX.append(percent_operation_HEX_ord_jun)
                    ord_aug_percent_operation_HEX.append(percent_operation_HEX_ord_aug)
                    dfw_may_percent_operation_HEX.append(percent_operation_HEX_dfw_may) 
                    ord_jul_percent_operation_HEX.append(percent_operation_HEX_ord_jul)                                       
                    dfw_sep_percent_operation_HEX.append(percent_operation_HEX_dfw_sep)                                                                                                                                                                                                                                                                     
                    dfw_jun_percent_operation_HEX.append(percent_operation_HEX_dfw_jun)
                    dfw_jul_percent_operation_HEX.append(percent_operation_HEX_dfw_jul)                       
                    
             

    ord_jan_time = np.hstack(ord_jan_time)
    ord_dec_time = np.hstack(ord_dec_time)
    ord_mar_time = np.hstack(ord_mar_time)
    ord_nov_time = np.hstack(ord_nov_time)
    dfw_jan_time = np.hstack(dfw_jan_time)
    dfw_dec_time = np.hstack(dfw_dec_time)
    dfw_feb_time = np.hstack(dfw_feb_time)
    ord_oct_time = np.hstack(ord_oct_time)
    dfw_nov_time = np.hstack(dfw_nov_time)
    ord_sep_time = np.hstack(ord_sep_time)
    dfw_oct_time = np.hstack(dfw_oct_time)
    ord_jun_time = np.hstack(ord_jun_time)
    ord_aug_time = np.hstack(ord_aug_time)
    dfw_may_time = np.hstack(dfw_may_time)
    ord_jul_time = np.hstack(ord_jul_time)
    dfw_sep_time = np.hstack(dfw_sep_time)
    dfw_jun_time = np.hstack(dfw_jun_time)  
    dfw_jul_time = np.hstack(dfw_jul_time)
    
    ord_jan_percent_operation_HAS = np.hstack(ord_jan_percent_operation_HAS)
    ord_dec_percent_operation_HAS = np.hstack(ord_dec_percent_operation_HAS)
    ord_mar_percent_operation_HAS = np.hstack(ord_mar_percent_operation_HAS)
    ord_nov_percent_operation_HAS = np.hstack(ord_nov_percent_operation_HAS)
    dfw_jan_percent_operation_HAS = np.hstack(dfw_jan_percent_operation_HAS)
    dfw_dec_percent_operation_HAS = np.hstack(dfw_dec_percent_operation_HAS)
    dfw_feb_percent_operation_HAS = np.hstack(dfw_feb_percent_operation_HAS)
    ord_oct_percent_operation_HAS = np.hstack(ord_oct_percent_operation_HAS)
    dfw_nov_percent_operation_HAS = np.hstack(dfw_nov_percent_operation_HAS)
    ord_sep_percent_operation_HAS = np.hstack(ord_sep_percent_operation_HAS)
    dfw_oct_percent_operation_HAS = np.hstack(dfw_oct_percent_operation_HAS)
    ord_jun_percent_operation_HAS = np.hstack(ord_jun_percent_operation_HAS)
    ord_aug_percent_operation_HAS = np.hstack(ord_aug_percent_operation_HAS)
    dfw_may_percent_operation_HAS = np.hstack(dfw_may_percent_operation_HAS)
    ord_jul_percent_operation_HAS = np.hstack(ord_jul_percent_operation_HAS)
    dfw_sep_percent_operation_HAS = np.hstack(dfw_sep_percent_operation_HAS)
    dfw_jun_percent_operation_HAS = np.hstack(dfw_jun_percent_operation_HAS)  
    dfw_jul_percent_operation_HAS = np.hstack(dfw_jul_percent_operation_HAS)
    
    ord_jan_percent_operation_HEX = np.hstack(ord_jan_percent_operation_HEX)
    ord_dec_percent_operation_HEX = np.hstack(ord_dec_percent_operation_HEX)
    ord_mar_percent_operation_HEX = np.hstack(ord_mar_percent_operation_HEX)
    ord_nov_percent_operation_HEX = np.hstack(ord_nov_percent_operation_HEX)
    dfw_jan_percent_operation_HEX = np.hstack(dfw_jan_percent_operation_HEX)
    dfw_dec_percent_operation_HEX = np.hstack(dfw_dec_percent_operation_HEX)
    dfw_feb_percent_operation_HEX = np.hstack(dfw_feb_percent_operation_HEX)
    ord_oct_percent_operation_HEX = np.hstack(ord_oct_percent_operation_HEX)
    dfw_nov_percent_operation_HEX = np.hstack(dfw_nov_percent_operation_HEX)
    ord_sep_percent_operation_HEX = np.hstack(ord_sep_percent_operation_HEX)
    dfw_oct_percent_operation_HEX = np.hstack(dfw_oct_percent_operation_HEX)
    ord_jun_percent_operation_HEX = np.hstack(ord_jun_percent_operation_HEX)
    ord_aug_percent_operation_HEX = np.hstack(ord_aug_percent_operation_HEX)
    dfw_may_percent_operation_HEX = np.hstack(dfw_may_percent_operation_HEX)
    ord_jul_percent_operation_HEX = np.hstack(ord_jul_percent_operation_HEX)
    dfw_sep_percent_operation_HEX = np.hstack(dfw_sep_percent_operation_HEX)
    dfw_jun_percent_operation_HEX = np.hstack(dfw_jun_percent_operation_HEX)  
    dfw_jul_percent_operation_HEX = np.hstack(dfw_jul_percent_operation_HEX)
    
    
    # get plotting style 
    ps = plot_style()        
    
    fig = plt.figure(save_filename)
    fig.set_size_inches(5, 8)    
    axis_0 = plt.subplot(1,1,1)
    axis_1 = plt.subplot(2,1,1)
    axis_2 = plt.subplot(2,1,2) 
    
    
    X = np.array((ord_jan_time, ord_dec_time, ord_mar_time,ord_nov_time, dfw_jan_time, dfw_dec_time, dfw_feb_time,   ord_oct_time, dfw_nov_time, ord_sep_time, dfw_oct_time, ord_jun_time, ord_aug_time, dfw_may_time, ord_jul_time, dfw_sep_time, dfw_jun_time, dfw_jul_time))  
    Y = np.tile(mean_temperature,(1, 272*no_of_fligths))  
    Z1 =  [ord_jan_percent_operation_HAS, ord_dec_percent_operation_HAS, ord_mar_percent_operation_HAS, ord_nov_percent_operation_HAS, dfw_jan_percent_operation_HAS, dfw_dec_percent_operation_HAS, dfw_feb_percent_operation_HAS, ord_oct_percent_operation_HAS, dfw_nov_percent_operation_HAS, ord_sep_percent_operation_HAS, dfw_oct_percent_operation_HAS, ord_jun_percent_operation_HAS, ord_aug_percent_operation_HAS, dfw_may_percent_operation_HAS, ord_jul_percent_operation_HAS, dfw_sep_percent_operation_HAS, dfw_jun_percent_operation_HAS, dfw_jul_percent_operation_HAS]
    Z2 =  [ord_jan_percent_operation_HEX, ord_dec_percent_operation_HEX, ord_mar_percent_operation_HEX, ord_nov_percent_operation_HEX, dfw_jan_percent_operation_HEX, dfw_dec_percent_operation_HEX, dfw_feb_percent_operation_HEX, ord_oct_percent_operation_HEX, dfw_nov_percent_operation_HEX, ord_sep_percent_operation_HEX, dfw_oct_percent_operation_HEX, ord_jun_percent_operation_HEX, ord_aug_percent_operation_HEX, dfw_may_percent_operation_HEX, ord_jul_percent_operation_HEX, dfw_sep_percent_operation_HEX, dfw_jun_percent_operation_HEX, dfw_jul_percent_operation_HEX]                                                                                                                                        
          
    axis_0.grid(False)
    axis_0.axis('off')     
    contour = axis_1.contourf(X, Y, Z1, cmap="summer", levels=256)
    cbar1 = fig.colorbar(contour,ticks=np.linspace(0, 100, 11))
    cbar1.set_label('Operating Power of HAS (%)', fontsize = ps.title_font_size-2)
    axis_1.set_xlabel('Flight Time (mins)', fontsize = ps.axis_font_size)
    axis_1.set_ylabel('Ambient Temperature ($\degree$C)', fontsize = ps.axis_font_size)
    #axis_1.set_title('Variation of Operating Power', fontsize = ps.title_font_size)
    
    contour = axis_2.contourf(X, Y, Z2, cmap="summer", levels=256)
    cbar2 = fig.colorbar(contour,ticks=np.linspace(0, 100, 11))
    cbar2.set_label('Operating Power of HEX (%)', fontsize = ps.title_font_size-2)
    axis_2.set_xlabel('Flight Time (mins)', fontsize = ps.axis_font_size)
    axis_2.set_ylabel('Ambient Temperature ($\degree$C)', fontsize = ps.axis_font_size)
    #axis_2.set_title('Variation of Operating Power', fontsize = ps.title_font_size)    
    
    
    fig.tight_layout()


    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.71, 0.725),loc='lower left')   

    if save_figure:
        plt.savefig(f'../Final_Plots/{save_filename}{file_type}') 

    return fig

# ------------------------------------------------------------------
#   Create Contour Plots for Battery Temperature
# ------------------------------------------------------------------
def batter_temperature_contour_plot (ord_jan, ord_dec, ord_mar, ord_nov, dfw_jan, dfw_dec,dfw_feb, ord_oct,
                                     dfw_nov, ord_sep, dfw_oct, ord_jun, ord_aug, dfw_may, ord_jul, dfw_sep, dfw_jun, dfw_jul, no_of_fligths, 
                   save_figure = True,
                   save_filename = "Battery_Temperature_plot_multiple_plots",
                   show_legend = False, 
                   file_type = ".png",
                   width = 12, height = 7):
    
     # get plotting style 
    ps = plot_style()        
    
    fig, axis = plt.subplots()
    fig.set_size_inches(5, 4)  

    
   
    ord_jan_time = []
    ord_dec_time = []
    ord_mar_time = []
    ord_nov_time = []
    dfw_jan_time = []
    dfw_dec_time = []
    dfw_feb_time = []
    ord_oct_time = []
    dfw_nov_time = []
    ord_sep_time = []
    dfw_oct_time = []
    ord_jun_time = []
    ord_aug_time = []
    dfw_may_time = []
    ord_jul_time = []
    dfw_sep_time = []
    dfw_jun_time = []
    dfw_jul_time = []
    
    ord_jan_battery_temp = []
    ord_dec_battery_temp = []
    ord_mar_battery_temp = []
    ord_nov_battery_temp = []
    dfw_jan_battery_temp = []
    dfw_dec_battery_temp = []
    dfw_feb_battery_temp = []
    ord_oct_battery_temp = []
    dfw_nov_battery_temp = []
    ord_sep_battery_temp = []
    dfw_oct_battery_temp = []
    ord_jun_battery_temp = []
    ord_aug_battery_temp = []
    dfw_may_battery_temp = []
    ord_jul_battery_temp = []
    dfw_sep_battery_temp = []
    dfw_jun_battery_temp = []
    dfw_jul_battery_temp = []       
 

    mean_temperature_ord_jan = ord_jan.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_dec = ord_dec.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_mar = ord_mar.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_nov = ord_nov.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jan = dfw_jan.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_dec = dfw_dec.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_feb = dfw_feb.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_oct = ord_oct.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_nov = dfw_nov.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_sep = ord_sep.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_oct = dfw_oct.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_jun = ord_jun.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_aug = ord_aug.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_may = dfw_may.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_jul = ord_jul.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_sep = dfw_sep.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jun = dfw_jun.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jul = dfw_jul.segments[0].conditions.freestream.temperature[0] - 273
    
    mean_temperature = np.array((mean_temperature_ord_jan, mean_temperature_ord_dec, mean_temperature_ord_mar,
                                 mean_temperature_ord_nov,  mean_temperature_dfw_jan,mean_temperature_dfw_dec, mean_temperature_dfw_feb ,
                                 mean_temperature_ord_oct, mean_temperature_dfw_nov, mean_temperature_ord_sep , mean_temperature_dfw_oct,
                                 mean_temperature_ord_jun, mean_temperature_ord_aug  , mean_temperature_dfw_may, mean_temperature_ord_jul ,
                                 mean_temperature_dfw_sep, mean_temperature_dfw_jun , mean_temperature_dfw_jul))      

    for network in ord_jan.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                for i in range(14*no_of_fligths):
                    
                    
                   
                    time_ord_jan = ord_jan.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_dec = ord_dec.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_mar = ord_mar.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_nov = ord_nov.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jan = dfw_jan.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_dec = dfw_dec.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_feb = dfw_feb.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_oct = ord_oct.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_nov = dfw_nov.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_sep = ord_sep.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_oct = dfw_oct.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_jun = ord_jun.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_aug = ord_aug.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_may = dfw_may.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_jul = ord_jul.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_sep = dfw_sep.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jun = dfw_jun.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jul = dfw_jul.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins                    
                    
                    
                    ord_jan_time.append(time_ord_jan)                               
                    ord_dec_time.append(time_ord_dec)                  
                    ord_mar_time.append(time_ord_mar)                    
                    ord_nov_time.append(time_ord_nov)
                    dfw_jan_time.append(time_dfw_jan)                                       
                    dfw_dec_time.append(time_dfw_dec)                                       
                    dfw_feb_time.append(time_dfw_feb)                                        
                    ord_oct_time.append(time_ord_oct)                                       
                    dfw_nov_time.append(time_dfw_nov)                                       
                    ord_sep_time.append(time_ord_sep)                                      
                    dfw_oct_time.append(time_dfw_oct)                 
                    ord_jun_time.append(time_ord_jun)
                    ord_aug_time.append(time_ord_aug)
                    dfw_may_time.append(time_dfw_may) 
                    ord_jul_time.append(time_ord_jul)                                       
                    dfw_sep_time.append(time_dfw_sep)                                                                                                                                                                                                                                                                     
                    dfw_jun_time.append(time_dfw_jun)
                    dfw_jul_time.append(time_dfw_jul)                       
                    
                    battery_temp_ord_jan = ord_jan.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_jan_battery_temp.append(battery_temp_ord_jan)
                    
                    battery_temp_ord_dec = ord_dec.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_dec_battery_temp.append(battery_temp_ord_dec)
                    
                    battery_temp_ord_mar = ord_mar.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_mar_battery_temp.append(battery_temp_ord_mar)
                    
                    battery_temp_ord_nov = ord_nov.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_nov_battery_temp.append(battery_temp_ord_nov)
                    
                    battery_temp_dfw_jan = dfw_jan.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_jan_battery_temp.append(battery_temp_dfw_jan)
                              
                    battery_temp_dfw_dec = dfw_dec.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_dec_battery_temp.append(battery_temp_dfw_dec)
                    
                    battery_temp_dfw_feb = dfw_feb.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_feb_battery_temp.append(battery_temp_dfw_feb)
                    
                    battery_temp_ord_oct = ord_oct.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_oct_battery_temp.append(battery_temp_ord_oct)
                    
                    battery_temp_dfw_nov = dfw_nov.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_nov_battery_temp.append(battery_temp_dfw_nov)
                    
                    battery_temp_ord_sep = ord_sep.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_sep_battery_temp.append(battery_temp_ord_sep)
                    
                    battery_temp_dfw_oct = dfw_oct.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_oct_battery_temp.append(battery_temp_dfw_oct)                  
                    
                    battery_temp_ord_jun = ord_jun.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_jun_battery_temp.append(battery_temp_ord_jun)
                
                    battery_temp_ord_aug = ord_aug.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_aug_battery_temp.append(battery_temp_ord_aug)
                    
                    battery_temp_dfw_may = dfw_may.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_may_battery_temp.append(battery_temp_dfw_may)
                    
                    battery_temp_ord_jul = ord_jul.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_jul_battery_temp.append(battery_temp_ord_jul)
                    
                    battery_temp_dfw_sep = dfw_sep.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_sep_battery_temp.append(battery_temp_dfw_sep)                                                                                                                                                                                                                                                                      
                    
                    battery_temp_dfw_jun = dfw_jun.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_jun_battery_temp.append(battery_temp_dfw_jun)
                  
                    battery_temp_dfw_jul = dfw_jul.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_jul_battery_temp.append(battery_temp_dfw_jul)
                    
                   
                    
             

   
   
    ord_jan_time = np.hstack(ord_jan_time)
    ord_dec_time = np.hstack(ord_dec_time)
    ord_mar_time = np.hstack(ord_mar_time)
    ord_nov_time = np.hstack(ord_nov_time)
    dfw_jan_time = np.hstack(dfw_jan_time)
    dfw_dec_time = np.hstack(dfw_dec_time)
    dfw_feb_time = np.hstack(dfw_feb_time)
    ord_oct_time = np.hstack(ord_oct_time)
    dfw_nov_time = np.hstack(dfw_nov_time)
    ord_sep_time = np.hstack(ord_sep_time)
    dfw_oct_time = np.hstack(dfw_oct_time)
    ord_jun_time = np.hstack(ord_jun_time)
    ord_aug_time = np.hstack(ord_aug_time)
    dfw_may_time = np.hstack(dfw_may_time)
    ord_jul_time = np.hstack(ord_jul_time)
    dfw_sep_time = np.hstack(dfw_sep_time)
    dfw_jun_time = np.hstack(dfw_jun_time)  
    dfw_jul_time = np.hstack(dfw_jul_time)    
    
            
    ord_jan_battery_temp = np.hstack(ord_jan_battery_temp)
    ord_dec_battery_temp = np.hstack(ord_dec_battery_temp)   
    ord_mar_battery_temp = np.hstack(ord_mar_battery_temp)
    ord_nov_battery_temp = np.hstack(ord_nov_battery_temp)
    dfw_jan_battery_temp = np.hstack(dfw_jan_battery_temp)
    dfw_dec_battery_temp = np.hstack(dfw_dec_battery_temp)
    dfw_feb_battery_temp = np.hstack(dfw_feb_battery_temp)
    ord_oct_battery_temp = np.hstack(ord_oct_battery_temp)
    dfw_nov_battery_temp = np.hstack(dfw_nov_battery_temp)
    ord_sep_battery_temp = np.hstack(ord_sep_battery_temp)
    dfw_oct_battery_temp = np.hstack(dfw_oct_battery_temp)
    ord_jun_battery_temp = np.hstack(ord_jun_battery_temp)
    ord_aug_battery_temp = np.hstack(ord_aug_battery_temp)
    dfw_may_battery_temp = np.hstack(dfw_may_battery_temp)
    ord_jul_battery_temp = np.hstack(ord_jul_battery_temp)
    dfw_sep_battery_temp = np.hstack(dfw_sep_battery_temp)
    dfw_jun_battery_temp = np.hstack(dfw_jun_battery_temp)  
    dfw_jul_battery_temp = np.hstack(dfw_jul_battery_temp)
    
    
    for network in ord_jan.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                for i in range(1):
                    ord_jan_idle_time  = ord_jan.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    ord_dec_idle_time  = ord_dec.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    ord_mar_idle_time  = ord_mar.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    ord_nov_idle_time  = ord_nov.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_jan_idle_time  = dfw_jan.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_dec_idle_time  = dfw_dec.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_feb_idle_time  = dfw_feb.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    ord_oct_idle_time  = ord_oct.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_nov_idle_time  = dfw_nov.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    ord_sep_idle_time  = ord_sep.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_oct_idle_time  = dfw_oct.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    ord_jun_idle_time  = ord_jun.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    ord_aug_idle_time  = ord_aug.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_may_idle_time  = dfw_may.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    ord_jul_idle_time  = ord_jul.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_sep_idle_time  = dfw_sep.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_jun_idle_time  = dfw_jun.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins
                    dfw_jul_idle_time  = dfw_jul.segments[i].conditions.frames.inertial.time[-1,0]/ Units.mins                    
                
                       
                    
        
    
      
    
    X = np.array((ord_jan_time, ord_dec_time, ord_mar_time,ord_nov_time, dfw_jan_time, dfw_dec_time, dfw_feb_time,   ord_oct_time, dfw_nov_time, ord_sep_time, dfw_oct_time, ord_jun_time, ord_aug_time, dfw_may_time, ord_jul_time, dfw_sep_time, dfw_jun_time, dfw_jul_time))  
    Y = np.tile(mean_temperature,(1, 272*no_of_fligths))    
    Z =  np.array((ord_jan_battery_temp, ord_dec_battery_temp, ord_mar_battery_temp,ord_nov_battery_temp, dfw_jan_battery_temp, dfw_dec_battery_temp, dfw_feb_battery_temp,   ord_oct_battery_temp, dfw_nov_battery_temp, ord_sep_battery_temp, dfw_oct_battery_temp, ord_jun_battery_temp, ord_aug_battery_temp, dfw_may_battery_temp, ord_jul_battery_temp, dfw_sep_battery_temp, dfw_jun_battery_temp, dfw_jul_battery_temp))
   
  
    contour = axis.contourf(X, Y, Z, cmap="magma", levels=128)
    cb = fig.colorbar(contour)
    cb.set_label('Battery Temperature ($\degree$C) ', fontsize = ps.axis_font_size)
    axis.set_xlabel('Flight Time (mins)', fontsize = ps.axis_font_size)
    axis.set_ylabel('Ambient Temperature ($\degree$C)', fontsize = ps.axis_font_size)
    #axis.set_title('Variation of Battery Temperature', fontsize = ps.title_font_size)
    
    #idle_time = np.array((ord_jan_idle_time, ord_dec_idle_time, ord_mar_idle_time, ord_nov_idle_time, dfw_jan_idle_time, dfw_dec_idle_time, dfw_feb_idle_time, ord_oct_idle_time, dfw_nov_idle_time, ord_sep_idle_time, dfw_oct_idle_time, ord_jun_idle_time,  ord_aug_idle_time, dfw_may_idle_time, ord_jul_idle_time, dfw_sep_idle_time, dfw_jun_idle_time, dfw_jul_idle_time))   
    #axis.plot(idle_time, mean_temperature, color='white')
    #axis.annotate('Idle\n Time', xy=(2, 1), xytext=(3, 0.5), color = 'white', fontsize = ps.axis_font_size)
    #axis.annotate('Flight\n Time', xy=(2, 1), xytext=(50, 2.5), color = 'white', fontsize = ps.axis_font_size)    
    
    fig.tight_layout()


    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.71, 0.725),loc='lower left')   

    if save_figure:
        plt.savefig(f'../Final_Plots/{save_filename}{file_type}') 

    return fig

def reservoir_temperature_contour_plot (ord_jan, ord_dec, ord_mar, ord_nov, dfw_jan, dfw_dec,dfw_feb, ord_oct,
                                     dfw_nov, ord_sep, dfw_oct, ord_jun, ord_aug, dfw_may, ord_jul, dfw_sep, dfw_jun, dfw_jul, no_of_fligths, 
                   save_figure = True,
                   save_filename = "Reservoir_Temperature_plot_flights",
                   show_legend = False, 
                   file_type = ".png",
                   width = 12, height = 7):
    
     # get plotting style 
    ps = plot_style()        
    
    fig, axis = plt.subplots()
    fig.set_size_inches(5, 4)  

    
   
    ord_jan_time = []
    ord_dec_time = []
    ord_mar_time = []
    ord_nov_time = []
    dfw_jan_time = []
    dfw_dec_time = []
    dfw_feb_time = []
    ord_oct_time = []
    dfw_nov_time = []
    ord_sep_time = []
    dfw_oct_time = []
    ord_jun_time = []
    ord_aug_time = []
    dfw_may_time = []
    ord_jul_time = []
    dfw_sep_time = []
    dfw_jun_time = []
    dfw_jul_time = []   
 
    ord_jan_battery_temp = []
    ord_dec_battery_temp = []
    ord_mar_battery_temp = []
    ord_nov_battery_temp = []
    dfw_jan_battery_temp = []
    dfw_dec_battery_temp = []
    dfw_feb_battery_temp = []
    ord_oct_battery_temp = []
    dfw_nov_battery_temp = []
    ord_sep_battery_temp = []
    dfw_oct_battery_temp = []
    ord_jun_battery_temp = []
    ord_aug_battery_temp = []
    dfw_may_battery_temp = []
    ord_jul_battery_temp = []
    dfw_sep_battery_temp = []
    dfw_jun_battery_temp = []
    dfw_jul_battery_temp = []
    
    mean_temperature_ord_jan = ord_jan.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_dec = ord_dec.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_mar = ord_mar.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_nov = ord_nov.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jan = dfw_jan.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_dec = dfw_dec.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_feb = dfw_feb.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_oct = ord_oct.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_nov = dfw_nov.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_sep = ord_sep.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_oct = dfw_oct.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_jun = ord_jun.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_aug = ord_aug.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_may = dfw_may.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_jul = ord_jul.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_sep = dfw_sep.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jun = dfw_jun.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jul = dfw_jul.segments[0].conditions.freestream.temperature[0] - 273
    
    mean_temperature = np.array((mean_temperature_ord_jan, mean_temperature_ord_dec, mean_temperature_ord_mar,
                                 mean_temperature_ord_nov,  mean_temperature_dfw_jan,mean_temperature_dfw_dec, mean_temperature_dfw_feb ,
                                 mean_temperature_ord_oct, mean_temperature_dfw_nov, mean_temperature_ord_sep , mean_temperature_dfw_oct,
                                 mean_temperature_ord_jun, mean_temperature_ord_aug  , mean_temperature_dfw_may, mean_temperature_ord_jul ,
                                 mean_temperature_dfw_sep, mean_temperature_dfw_jun , mean_temperature_dfw_jul))      

    for network in ord_jan.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                for i in range(14*no_of_fligths):
                    
                    
                   
                    time_ord_jan = ord_jan.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_dec = ord_dec.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_mar = ord_mar.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_nov = ord_nov.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jan = dfw_jan.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_dec = dfw_dec.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_feb = dfw_feb.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_oct = ord_oct.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_nov = dfw_nov.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_sep = ord_sep.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_oct = dfw_oct.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_jun = ord_jun.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_aug = ord_aug.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_may = dfw_may.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_ord_jul = ord_jul.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_sep = dfw_sep.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jun = dfw_jun.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins
                    time_dfw_jul = dfw_jul.segments[i].conditions.frames.inertial.time[:,0]/ Units.mins                    
                    
                    
                    ord_jan_time.append(time_ord_jan)                               
                    ord_dec_time.append(time_ord_dec)                  
                    ord_mar_time.append(time_ord_mar)                    
                    ord_nov_time.append(time_ord_nov)
                    dfw_jan_time.append(time_dfw_jan)                                       
                    dfw_dec_time.append(time_dfw_dec)                                       
                    dfw_feb_time.append(time_dfw_feb)                                        
                    ord_oct_time.append(time_ord_oct)                                       
                    dfw_nov_time.append(time_dfw_nov)                                       
                    ord_sep_time.append(time_ord_sep)                                      
                    dfw_oct_time.append(time_dfw_oct)                 
                    ord_jun_time.append(time_ord_jun)
                    ord_aug_time.append(time_ord_aug)
                    dfw_may_time.append(time_dfw_may) 
                    ord_jul_time.append(time_ord_jul)                                       
                    dfw_sep_time.append(time_dfw_sep)                                                                                                                                                                                                                                                                     
                    dfw_jun_time.append(time_dfw_jun)
                    dfw_jul_time.append(time_dfw_jul)                       
                    
                    battery_temp_ord_jan = ord_jan.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_jan_battery_temp.append(battery_temp_ord_jan)
                    
                    battery_temp_ord_dec = ord_dec.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_dec_battery_temp.append(battery_temp_ord_dec)
                    
                    battery_temp_ord_mar = ord_mar.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_mar_battery_temp.append(battery_temp_ord_mar)
                    
                    battery_temp_ord_nov = ord_nov.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_nov_battery_temp.append(battery_temp_ord_nov)
                    
                    battery_temp_dfw_jan = dfw_jan.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_jan_battery_temp.append(battery_temp_dfw_jan)
                              
                    battery_temp_dfw_dec = dfw_dec.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_dec_battery_temp.append(battery_temp_dfw_dec)
                    
                    battery_temp_dfw_feb = dfw_feb.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_feb_battery_temp.append(battery_temp_dfw_feb)
                    
                    battery_temp_ord_oct = ord_oct.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_oct_battery_temp.append(battery_temp_ord_oct)
                    
                    battery_temp_dfw_nov = dfw_nov.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_nov_battery_temp.append(battery_temp_dfw_nov)
                    
                    battery_temp_ord_sep = ord_sep.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_sep_battery_temp.append(battery_temp_ord_sep)
                    
                    battery_temp_dfw_oct = dfw_oct.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_oct_battery_temp.append(battery_temp_dfw_oct)                  
                    
                    battery_temp_ord_jun = ord_jun.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_jun_battery_temp.append(battery_temp_ord_jun)
                
                    battery_temp_ord_aug = ord_aug.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_aug_battery_temp.append(battery_temp_ord_aug)
                    
                    battery_temp_dfw_may = dfw_may.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_may_battery_temp.append(battery_temp_dfw_may)
                    
                    battery_temp_ord_jul = ord_jul.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    ord_jul_battery_temp.append(battery_temp_ord_jul)
                    
                    battery_temp_dfw_sep = dfw_sep.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_sep_battery_temp.append(battery_temp_dfw_sep)                                                                                                                                                                                                                                                                      
                    
                    battery_temp_dfw_jun = dfw_jun.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_jun_battery_temp.append(battery_temp_dfw_jun)
                  
                    battery_temp_dfw_jul = dfw_jul.segments[i].conditions.energy[bus.tag][battery.tag].thermal_management_system.RES.coolant_temperature[:,0] - 273
                    dfw_jul_battery_temp.append(battery_temp_dfw_jul)
                    
                   
                    
             

   
   
    ord_jan_time = np.hstack(ord_jan_time)
    ord_dec_time = np.hstack(ord_dec_time)
    ord_mar_time = np.hstack(ord_mar_time)
    ord_nov_time = np.hstack(ord_nov_time)
    dfw_jan_time = np.hstack(dfw_jan_time)
    dfw_dec_time = np.hstack(dfw_dec_time)
    dfw_feb_time = np.hstack(dfw_feb_time)
    ord_oct_time = np.hstack(ord_oct_time)
    dfw_nov_time = np.hstack(dfw_nov_time)
    ord_sep_time = np.hstack(ord_sep_time)
    dfw_oct_time = np.hstack(dfw_oct_time)
    ord_jun_time = np.hstack(ord_jun_time)
    ord_aug_time = np.hstack(ord_aug_time)
    dfw_may_time = np.hstack(dfw_may_time)
    ord_jul_time = np.hstack(ord_jul_time)
    dfw_sep_time = np.hstack(dfw_sep_time)
    dfw_jun_time = np.hstack(dfw_jun_time)  
    dfw_jul_time = np.hstack(dfw_jul_time)    
    
            
    ord_jan_battery_temp = np.hstack(ord_jan_battery_temp)
    ord_dec_battery_temp = np.hstack(ord_dec_battery_temp)   
    ord_mar_battery_temp = np.hstack(ord_mar_battery_temp)
    ord_nov_battery_temp = np.hstack(ord_nov_battery_temp)
    dfw_jan_battery_temp = np.hstack(dfw_jan_battery_temp)
    dfw_dec_battery_temp = np.hstack(dfw_dec_battery_temp)
    dfw_feb_battery_temp = np.hstack(dfw_feb_battery_temp)
    ord_oct_battery_temp = np.hstack(ord_oct_battery_temp)
    dfw_nov_battery_temp = np.hstack(dfw_nov_battery_temp)
    ord_sep_battery_temp = np.hstack(ord_sep_battery_temp)
    dfw_oct_battery_temp = np.hstack(dfw_oct_battery_temp)
    ord_jun_battery_temp = np.hstack(ord_jun_battery_temp)
    ord_aug_battery_temp = np.hstack(ord_aug_battery_temp)
    dfw_may_battery_temp = np.hstack(dfw_may_battery_temp)
    ord_jul_battery_temp = np.hstack(ord_jul_battery_temp)
    dfw_sep_battery_temp = np.hstack(dfw_sep_battery_temp)
    dfw_jun_battery_temp = np.hstack(dfw_jun_battery_temp)  
    dfw_jul_battery_temp = np.hstack(dfw_jul_battery_temp)
    
      
    
    X = np.array((ord_jan_time, ord_dec_time, ord_mar_time,ord_nov_time, dfw_jan_time, dfw_dec_time, dfw_feb_time,   ord_oct_time, dfw_nov_time, ord_sep_time, dfw_oct_time, ord_jun_time, ord_aug_time, dfw_may_time, ord_jul_time, dfw_sep_time, dfw_jun_time, dfw_jul_time))  
    Y = np.tile(mean_temperature,(1, 272*no_of_fligths))    
    Z =  np.array((ord_jan_battery_temp, ord_dec_battery_temp, ord_mar_battery_temp,ord_nov_battery_temp, dfw_jan_battery_temp, dfw_dec_battery_temp, dfw_feb_battery_temp,   ord_oct_battery_temp, dfw_nov_battery_temp, ord_sep_battery_temp, dfw_oct_battery_temp, ord_jun_battery_temp, ord_aug_battery_temp, dfw_may_battery_temp, ord_jul_battery_temp, dfw_sep_battery_temp, dfw_jun_battery_temp, dfw_jul_battery_temp))
   
  
    contour = axis.contourf(X, Y, Z, cmap="magma", levels=128)
    cb = fig.colorbar(contour)
    cb.set_label('Reservoir Coolant Temperature ($\degree$C) ', fontsize = ps.axis_font_size)
    axis.set_xlabel('Flight Time (mins)', fontsize = ps.axis_font_size)
    axis.set_ylabel('Ambient Temperature ($\degree$C)', fontsize = ps.axis_font_size)
    #axis.set_title('Variation of Battery Temperature', fontsize = ps.title_font_size)
    fig.tight_layout()


    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.71, 0.725),loc='lower left')   

    if save_figure:
        plt.savefig(f'../Final_Plots/{save_filename}{file_type}') 

    return fig
# ------------------------------------------------------------------
#   Create Contour Plots for Battery Temperature and Range Variation
# ------------------------------------------------------------------
def four_dimensional_contour_plot (ord_jan, ord_dec, ord_mar, ord_nov, dfw_jan, dfw_dec,dfw_feb, ord_oct,
                                     dfw_nov, ord_sep, dfw_oct, ord_jun, ord_aug, dfw_may, ord_jul, dfw_sep, dfw_jun, dfw_jul, no_flights, 
                   save_figure = True,
                   save_filename = "Battery_Range_Temperature_plot_",
                   show_legend = False, 
                   file_type = ".png"):
    
       
    ord_jan_time = []
    ord_dec_time = []
    ord_mar_time = []
    ord_nov_time = []
    dfw_jan_time = []
    dfw_dec_time = []
    dfw_feb_time = []
    ord_oct_time = []
    dfw_nov_time = []
    ord_sep_time = []
    dfw_oct_time = []
    ord_jun_time = []
    ord_aug_time = []
    dfw_may_time = []
    ord_jul_time = []
    dfw_sep_time = []
    dfw_jun_time = []
    dfw_jul_time = []
    
    ord_jan_battery_temp = []
    ord_dec_battery_temp = []
    ord_mar_battery_temp = []
    ord_nov_battery_temp = []
    dfw_jan_battery_temp = []
    dfw_dec_battery_temp = []
    dfw_feb_battery_temp = []
    ord_oct_battery_temp = []
    dfw_nov_battery_temp = []
    ord_sep_battery_temp = []
    dfw_oct_battery_temp = []
    ord_jun_battery_temp = []
    ord_aug_battery_temp = []
    dfw_may_battery_temp = []
    ord_jul_battery_temp = []
    dfw_sep_battery_temp = []
    dfw_jun_battery_temp = []
    dfw_jul_battery_temp = []
    
    ord_jan_range = []
    ord_dec_range = []
    ord_mar_range = []
    ord_nov_range = []
    dfw_jan_range = []
    dfw_dec_range = []
    dfw_feb_range = []
    ord_oct_range = []
    dfw_nov_range = []
    ord_sep_range = []
    dfw_oct_range = []
    ord_jun_range = []
    ord_aug_range = []
    dfw_may_range = []
    ord_jul_range = []
    dfw_sep_range = []
    dfw_jun_range = []
    dfw_jul_range = []    
    
    mean_temperature_ord_jan = ord_jan.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_dec = ord_dec.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_mar = ord_mar.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_nov = ord_nov.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jan = dfw_jan.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_dec = dfw_dec.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_feb = dfw_feb.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_oct = ord_oct.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_nov = dfw_nov.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_sep = ord_sep.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_oct = dfw_oct.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_jun = ord_jun.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_aug = ord_aug.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_may = dfw_may.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_ord_jul = ord_jul.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_sep = dfw_sep.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jun = dfw_jun.segments[0].conditions.freestream.temperature[0] - 273
    mean_temperature_dfw_jul = dfw_jul.segments[0].conditions.freestream.temperature[0] - 273
    
    mean_temperature = np.array((mean_temperature_ord_jan, mean_temperature_ord_dec, mean_temperature_ord_mar,
                                 mean_temperature_ord_nov,  mean_temperature_dfw_jan,mean_temperature_dfw_dec, mean_temperature_dfw_feb ,
                                 mean_temperature_ord_oct, mean_temperature_dfw_nov, mean_temperature_ord_sep , mean_temperature_dfw_oct,
                                 mean_temperature_ord_jun, mean_temperature_ord_aug  , mean_temperature_dfw_may, mean_temperature_ord_jul ,
                                 mean_temperature_dfw_sep, mean_temperature_dfw_jun , mean_temperature_dfw_jul))      

    for network in ord_jan.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:   
                for i in range(14*no_flights): 
                    
                    time_ord_jan = ord_jan.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_ord_dec = ord_dec.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_ord_mar = ord_mar.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_ord_nov = ord_nov.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_jan = dfw_jan.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_dec = dfw_dec.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_feb = dfw_feb.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_ord_oct = ord_oct.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_nov = dfw_nov.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_ord_sep = ord_sep.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_oct = dfw_oct.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_ord_jun = ord_jun.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_ord_aug = ord_aug.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_may = dfw_may.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_ord_jul = ord_jul.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_sep = dfw_sep.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_jun = dfw_jun.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs
                    time_dfw_jul = dfw_jul.segments[i].conditions.frames.inertial.time[:,0]/ Units.hrs                    
                    
                    
                    ord_jan_time.append(time_ord_jan)                               
                    ord_dec_time.append(time_ord_dec)                  
                    ord_mar_time.append(time_ord_mar)                    
                    ord_nov_time.append(time_ord_nov)
                    dfw_jan_time.append(time_dfw_jan)                                       
                    dfw_dec_time.append(time_dfw_dec)                                       
                    dfw_feb_time.append(time_dfw_feb)                                        
                    ord_oct_time.append(time_ord_oct)                                       
                    dfw_nov_time.append(time_dfw_nov)                                       
                    ord_sep_time.append(time_ord_sep)                                      
                    dfw_oct_time.append(time_dfw_oct)                 
                    ord_jun_time.append(time_ord_jun)
                    ord_aug_time.append(time_ord_aug)
                    dfw_may_time.append(time_dfw_may) 
                    ord_jul_time.append(time_ord_jul)                                       
                    dfw_sep_time.append(time_dfw_sep)                                                                                                                                                                                                                                                                     
                    dfw_jun_time.append(time_dfw_jun)
                    dfw_jul_time.append(time_dfw_jul)
                    
                    battery_temp_ord_jan = ord_jan.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_jan_battery_temp.append(battery_temp_ord_jan)
                    
                    battery_temp_ord_dec = ord_dec.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_dec_battery_temp.append(battery_temp_ord_dec)
                    
                    battery_temp_ord_mar = ord_mar.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_mar_battery_temp.append(battery_temp_ord_mar)
                    
                    battery_temp_ord_nov = ord_nov.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_nov_battery_temp.append(battery_temp_ord_nov)
                    
                    battery_temp_dfw_jan = dfw_jan.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_jan_battery_temp.append(battery_temp_dfw_jan)
                              
                    battery_temp_dfw_dec = dfw_dec.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_dec_battery_temp.append(battery_temp_dfw_dec)
                    
                    battery_temp_dfw_feb = dfw_feb.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_feb_battery_temp.append(battery_temp_dfw_feb)
                    
                    battery_temp_ord_oct = ord_oct.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_oct_battery_temp.append(battery_temp_ord_oct)
                    
                    battery_temp_dfw_nov = dfw_nov.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_nov_battery_temp.append(battery_temp_dfw_nov)
                    
                    battery_temp_ord_sep = ord_sep.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_sep_battery_temp.append(battery_temp_ord_sep)
                    
                    battery_temp_dfw_oct = dfw_oct.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_oct_battery_temp.append(battery_temp_dfw_oct)                  
                    
                    battery_temp_ord_jun = ord_jun.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_jun_battery_temp.append(battery_temp_ord_jun)
                
                    battery_temp_ord_aug = ord_aug.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_aug_battery_temp.append(battery_temp_ord_aug)
                    
                    battery_temp_dfw_may = dfw_may.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_may_battery_temp.append(battery_temp_dfw_may)
                    
                    battery_temp_ord_jul = ord_jul.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    ord_jul_battery_temp.append(battery_temp_ord_jul)
                    
                    battery_temp_dfw_sep = dfw_sep.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_sep_battery_temp.append(battery_temp_dfw_sep)                                                                                                                                                                                                                                                                      
                    
                    battery_temp_dfw_jun = dfw_jun.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_jun_battery_temp.append(battery_temp_dfw_jun)
                  
                    battery_temp_dfw_jul = dfw_jul.segments[i].conditions.energy[bus.tag][battery.tag].pack.temperature[:,0] - 273
                    dfw_jul_battery_temp.append(battery_temp_dfw_jul)
                    
                    # RANGE Calculation 
                                        
                    range_ord_jan = ord_jan.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_ord_dec = ord_dec.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_ord_mar = ord_mar.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_ord_nov = ord_nov.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_jan = dfw_jan.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_dec = dfw_dec.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_feb = dfw_feb.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_ord_oct = ord_oct.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_nov = dfw_nov.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_ord_sep = ord_sep.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_oct = dfw_oct.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_ord_jun = ord_jun.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_ord_aug = ord_aug.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_may = dfw_may.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_ord_jul = ord_jul.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_sep = dfw_sep.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_jun = dfw_jun.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    range_dfw_jul = dfw_jul.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi
                    
                    ord_jan_range.append(range_ord_jan)                               
                    ord_dec_range.append(range_ord_dec)                  
                    ord_mar_range.append(range_ord_mar)                    
                    ord_nov_range.append(range_ord_nov)
                    dfw_jan_range.append(range_dfw_jan)                                       
                    dfw_dec_range.append(range_dfw_dec)                                       
                    dfw_feb_range.append(range_dfw_feb)                                        
                    ord_oct_range.append(range_ord_oct)                                       
                    dfw_nov_range.append(range_dfw_nov)                                       
                    ord_sep_range.append(range_ord_sep)                                      
                    dfw_oct_range.append(range_dfw_oct)                 
                    ord_jun_range.append(range_ord_jun)
                    ord_aug_range.append(range_ord_aug)
                    dfw_may_range.append(range_dfw_may) 
                    ord_jul_range.append(range_ord_jul)                                       
                    dfw_sep_range.append(range_dfw_sep)                                                                                                                                                                                                                                                                     
                    dfw_jun_range.append(range_dfw_jun)
                    dfw_jul_range.append(range_dfw_jul)                        
                                        
                    
             

    ord_jan_time = np.hstack(ord_jan_time)
    ord_dec_time = np.hstack(ord_dec_time)
    ord_mar_time = np.hstack(ord_mar_time)
    ord_nov_time = np.hstack(ord_nov_time)
    dfw_jan_time = np.hstack(dfw_jan_time)
    dfw_dec_time = np.hstack(dfw_dec_time)
    dfw_feb_time = np.hstack(dfw_feb_time)
    ord_oct_time = np.hstack(ord_oct_time)
    dfw_nov_time = np.hstack(dfw_nov_time)
    ord_sep_time = np.hstack(ord_sep_time)
    dfw_oct_time = np.hstack(dfw_oct_time)
    ord_jun_time = np.hstack(ord_jun_time)
    ord_aug_time = np.hstack(ord_aug_time)
    dfw_may_time = np.hstack(dfw_may_time)
    ord_jul_time = np.hstack(ord_jul_time)
    dfw_sep_time = np.hstack(dfw_sep_time)
    dfw_jun_time = np.hstack(dfw_jun_time)  
    dfw_jul_time = np.hstack(dfw_jul_time)
    
    ord_jan_battery_temp = np.hstack(ord_jan_battery_temp)
    ord_dec_battery_temp = np.hstack(ord_dec_battery_temp)
    ord_mar_battery_temp = np.hstack(ord_mar_battery_temp)
    ord_nov_battery_temp = np.hstack(ord_nov_battery_temp)
    dfw_jan_battery_temp = np.hstack(dfw_jan_battery_temp)
    dfw_dec_battery_temp = np.hstack(dfw_dec_battery_temp)
    dfw_feb_battery_temp = np.hstack(dfw_feb_battery_temp)
    ord_oct_battery_temp = np.hstack(ord_oct_battery_temp)
    dfw_nov_battery_temp = np.hstack(dfw_nov_battery_temp)
    ord_sep_battery_temp = np.hstack(ord_sep_battery_temp)
    dfw_oct_battery_temp = np.hstack(dfw_oct_battery_temp)
    ord_jun_battery_temp = np.hstack(ord_jun_battery_temp)
    ord_aug_battery_temp = np.hstack(ord_aug_battery_temp)
    dfw_may_battery_temp = np.hstack(dfw_may_battery_temp)
    ord_jul_battery_temp = np.hstack(ord_jul_battery_temp)
    dfw_sep_battery_temp = np.hstack(dfw_sep_battery_temp)
    dfw_jun_battery_temp = np.hstack(dfw_jun_battery_temp)  
    dfw_jul_battery_temp = np.hstack(dfw_jul_battery_temp)
    
    ord_jan_range = np.hstack(ord_jan_range)
    ord_dec_range = np.hstack(ord_dec_range)
    ord_mar_range = np.hstack(ord_mar_range)
    ord_nov_range = np.hstack(ord_nov_range)
    dfw_jan_range = np.hstack(dfw_jan_range)
    dfw_dec_range = np.hstack(dfw_dec_range)
    dfw_feb_range = np.hstack(dfw_feb_range)
    ord_oct_range = np.hstack(ord_oct_range)
    dfw_nov_range = np.hstack(dfw_nov_range)
    ord_sep_range = np.hstack(ord_sep_range)
    dfw_oct_range = np.hstack(dfw_oct_range)
    ord_jun_range = np.hstack(ord_jun_range)
    ord_aug_range = np.hstack(ord_aug_range)
    dfw_may_range = np.hstack(dfw_may_range)
    ord_jul_range = np.hstack(ord_jul_range)
    dfw_sep_range = np.hstack(dfw_sep_range)
    dfw_jun_range = np.hstack(dfw_jun_range)  
    dfw_jul_range = np.hstack(dfw_jul_range)    
    
    
     
    
    
  
    
    X = np.array((ord_jan_time, ord_dec_time, ord_mar_time,ord_nov_time, dfw_jan_time, dfw_dec_time, dfw_feb_time,   ord_oct_time, dfw_nov_time, ord_sep_time, dfw_oct_time, ord_jun_time, ord_aug_time, dfw_may_time, ord_jul_time, dfw_sep_time, dfw_jun_time, dfw_jul_time))  
    Y = np.tile(mean_temperature,(1, 272*no_flights))    
    Z1 =  np.array((ord_jan_battery_temp, ord_dec_battery_temp, ord_mar_battery_temp,ord_nov_battery_temp, dfw_jan_battery_temp, dfw_dec_battery_temp, dfw_feb_battery_temp,   ord_oct_battery_temp, dfw_nov_battery_temp, ord_sep_battery_temp, dfw_oct_battery_temp, ord_jun_battery_temp, ord_aug_battery_temp, dfw_may_battery_temp, ord_jul_battery_temp, dfw_sep_battery_temp, dfw_jun_battery_temp, dfw_jul_battery_temp))
    Z2 = np.array((ord_jan_range, ord_dec_range, ord_mar_range,ord_nov_range, dfw_jan_range, dfw_dec_range, dfw_feb_range, ord_oct_range, dfw_nov_range, ord_sep_range, dfw_oct_range, ord_jun_range, ord_aug_range, dfw_may_range, ord_jul_range, dfw_sep_range, dfw_jun_range, dfw_jul_range))

    # get plotting style 
    ps = plot_style()          
                                                                        
    # Set up plot
    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
    fig.set_size_inches(9, 7)     
    ax.view_init(elev=40, azim=-120, roll=0)
    ls = LightSource(270, 45)
    # To use a custom hillshading mode, override the built-in shading and pass
    # in the rgb colors of the shaded surface calculated from "shade".
    rgb = ls.shade(Z2, cmap=cm.turbo, vert_exag=0.25, blend_mode='overlay')
    surf = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, facecolors=rgb,
                           linewidth=0, antialiased=False, shade=False)
    
    mappable = cm.ScalarMappable(cmap=cm.turbo)
    mappable.set_array(Z2)
    cbar = fig.colorbar(mappable, ax=ax,  shrink=0.65)
    cbar.set_label('Distance Traveled (miles)', fontsize = ps.axis_font_size+2)    
   
   
       
    ax.set_zlim(0, 40)
    ax.set_ylim(-10, 40) 
    ax.set_xlabel('Operational Flight Time (Hours)', fontsize = ps.axis_font_size +2)
    ax.set_ylabel('Ambient Temperature ($\degree$C)', fontsize = ps.axis_font_size+2)
    ax.set_zlabel('Battery Temperature ($\degree$C)', fontsize = ps.axis_font_size+2)
    
    fig.tight_layout()


    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.71, 0.725),loc='lower left')   

    if save_figure:
        plt.savefig(f'../Final_Plots/{save_filename}{no_flights}{file_type}') 

    return fig


# ------------------------------------------------------------------
#   Create Contour Plots for city
# ------------------------------------------------------------------
def city_plots (save_figure = True,
                              save_filename = "LAX_range_contour",
                              file_type = ".png"):
    
    
    # Load the "Temperature Ordered" sheet
    range_data_df = pd.read_excel('../Simulation_plan.xlsx', sheet_name='Temperature_Ordered')

    flight_1_data_temp_ordered = range_data_df[range_data_df['Flight_No'] == 1]
    mean_temperatures_flight_1 = flight_1_data_temp_ordered['Mean_Temperature'].values
    ranges_flight_1 = (flight_1_data_temp_ordered['Range_miles'].values)
   
    azimuths = np.radians(np.linspace(0, 360, 20))
    zeniths = ranges_flight_1 
    
    r, theta = np.meshgrid(zeniths, azimuths)
    values = []
    for i in range (20):
        values.append(np.transpose(mean_temperatures_flight_1))
    
    
    
    
    #----------Plot Range-----------------    
    ps = plot_style()
    
    fig = plt.figure(figsize=(8, 6))
    # Set up the polar plot
    polar_ax = fig.add_axes([0.167, 0.08, 0.668, 0.84], projection='polar')  
    polar_ax.patch.set_alpha(0.3)
    contour = polar_ax.contourf(theta, r, values, cmap="seismic", levels=256, alpha=0.2)
    polar_ax.set_xticks([])  # Remove azimuthal ticks
    polar_ax.set_xlabel('Range (miles) ')

    #------------ Add colorbar-------------
    cbar_ax = fig.add_axes([0.1, 0.1, 0.03, 0.5])  
    cbar = fig.colorbar(contour, cax=cbar_ax, label='Ambient Temperature ($\degree$C)', location='left')
    
    
    #-------------Overlay the background image---------------
    bg_image = plt.imread("folium_maps/LAX_map.png")
    background_ax = plt.axes([0, 0, 1, 1]) 
    background_ax.set_zorder(-1) 
    background_ax.imshow(bg_image, aspect='auto')   
    
    
    if save_figure:
        plt.savefig(f'../Final_Plots/{save_filename}{file_type}') 

    return fig    
    
def multi_city_plots (save_figure = True,
                              save_filename = "midwest_range_contour",
                              file_type = ".png"):
    
    
    # Load the "Temperature Ordered" sheet
    range_data_df = pd.read_excel('../Simulation_plan.xlsx', sheet_name='Temperature_Ordered')

    flight_1_data_temp_ordered = range_data_df[range_data_df['Flight_No'] == 1]
    mean_temperatures_flight_1 = flight_1_data_temp_ordered['Mean_Temperature'].values
    ranges_flight_1 = (flight_1_data_temp_ordered['Range_miles'].values)
   
    azimuths = np.radians(np.linspace(0, 360, 20))
    zeniths = ranges_flight_1 
    
    r, theta = np.meshgrid(zeniths, azimuths)
    values = []
    for i in range (20):
        values.append(np.transpose(mean_temperatures_flight_1))
    
    
    
    
    #----------Plot Range-----------------    
    ps = plot_style()
    
    fig = plt.figure(figsize=(8, 6))
    # Set up the polar plot
    polar_ax1 = fig.add_axes([0.3175, 0.1274, 0.145, 0.37], projection='polar')  
    polar_ax1.patch.set_alpha(0.1)
    contour = polar_ax1.contourf(theta, r, values, cmap="seismic", levels=256, alpha=0.2)
    polar_ax1.set_xticks([])  # Remove azimuthal ticks
    polar_ax1.set_yticks([])  # Remove azimuthal ticks
    #polar_ax1.set_xlabel('Range (miles) ')
    
    polar_ax2 = fig.add_axes([0.435, 0.195, 0.145, 0.37], projection='polar')  
    polar_ax2.patch.set_alpha(0.1)
    contour = polar_ax2.contourf(theta, r, values, cmap="seismic", levels=256, alpha=0.2)
    polar_ax2.set_xticks([])  # Remove azimuthal ticks
    polar_ax2.set_yticks([])  # Remove azimuthal ticks
    
    polar_ax3 = fig.add_axes([0.475, 0.32, 0.145, 0.37], projection='polar')  
    polar_ax3.patch.set_alpha(0.3)
    contour = polar_ax3.contourf(theta, r, values, cmap="seismic", levels=256, alpha=0.2)
    polar_ax3.set_xticks([])  # Remove azimuthal ticks
    polar_ax3.set_yticks([])  # Remove azimuthal ticks
    
    polar_ax4 = fig.add_axes([0.515, 0.37, 0.145, 0.37], projection='polar')  
    polar_ax4.patch.set_alpha(0.1)
    contour = polar_ax4.contourf(theta, r, values, cmap="seismic", levels=256, alpha=0.2)
    polar_ax4.set_xticks([])  # Remove azimuthal ticks
    polar_ax4.set_yticks([])  # Remove azimuthal ticks
    polar_ax5 = fig.add_axes([0.2975, 0.3557, 0.145, 0.37], projection='polar')  
    polar_ax5.patch.set_alpha(0.1)
    contour = polar_ax5.contourf(theta, r, values, cmap="seismic", levels=256, alpha=0.2)
    polar_ax5.set_xticks([])  # Remove azimuthal ticks
    polar_ax5.set_yticks([])  # Remove azimuthal ticks
    polar_ax6 = fig.add_axes([0.375, 0.44, 0.145, 0.37], projection='polar')  
    polar_ax6.patch.set_alpha(0.1)
    contour = polar_ax6.contourf(theta, r, values, cmap="seismic", levels=256, alpha=0.2)
    polar_ax6.set_xticks([])  # Remove azimuthal ticks
    polar_ax6.set_yticks([])  # Remove azimuthal ticks
    polar_ax7 = fig.add_axes([0.5775, 0.435, 0.145, 0.37], projection='polar')  
    polar_ax7.patch.set_alpha(0.1)
    contour = polar_ax7.contourf(theta, r, values, cmap="seismic", levels=256, alpha=0.2)
    polar_ax7.set_xticks([])  # Remove azimuthal ticks     
    polar_ax7.set_yticks([])  # Remove azimuthal ticks
    #------------ Add colorbar-------------
    cbar_ax = fig.add_axes([0.94, 0.1, 0.03, 0.5])  
    fig.colorbar(contour, cax=cbar_ax, label='Ambient Temperature ($\degree$C)', location='left', alpha = 1)
    
    
    #-------------Overlay the background image---------------
    bg_image = plt.imread("folium_maps/midwest_map.png")
    background_ax = plt.axes([0, 0, 1, 1]) 
    background_ax.set_zorder(-1) 
    background_ax.imshow(bg_image, aspect='auto')   
    
    
    if save_figure:
        plt.savefig(f'../Final_Plots/{save_filename}{file_type}') 

    return fig    
        
   
# ------------------------------------------------------------------
#   Monthly Temperature Variation 2023
# ------------------------------------------------------------------
# https://www.weather.gov/fwd/dmotemp

def monthly_weather(save_figure = True, show_legend = True,
                    save_filename = "monthly_weather_plot", file_type = ".png"):
    
    # get plotting style 
    ps = plot_style()    

    
    # get line colors for plots 
    
     
    line_colors   = cm.inferno(np.linspace(0,0.9,5))
    parameters = {'axes.labelsize': ps.axis_font_size+12,
                      'xtick.labelsize': ps.axis_font_size + 12,
                      'ytick.labelsize': ps.axis_font_size + 12,
                      'axes.titlesize': ps.title_font_size + 12}
    plt.rcParams.update(parameters)
    fig, ax = plt.subplots()  
    fig.set_size_inches(10, 6)         
    # Months and temperatures in Celsius
    cities = ["DFW", "ORD", "SFO", "JFK", "LAX"]
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    dfw_temps =  np.array((11.2, 11.3, 16.0, 18.5, 24.2, 28.9, 31.8, 33.2, 29.2, 20.7, 14.3, 11.6))
    ord_temps =  np.array((-1.6, 0.8, 5.8, 11.1, 17.2, 22.0, 24.3, 23.7, 19.0, 11.7, 4.4, -0.8))
    sfo_temps =  np.array((11.1, 12.2, 12.8, 13.9, 15.0, 16.7, 17.8, 18.3, 18.0, 16.7, 13.9, 11.7))
    jfk_temps =  np.array((0.6, 1.7, 6.1, 11.1, 16.7, 21.7, 25.0, 24.4, 20.0, 14.4, 8.9, 3.3))
    lax_temps =  np.array((14.4, 15.0, 16.1, 17.2, 18.3, 20.0, 21.1, 21.7, 21.1, 19.4, 17.2, 14.4))
    temperatures = np.array((dfw_temps, ord_temps, sfo_temps, jfk_temps, lax_temps))
    

    ax.plot(months, dfw_temps, label="DFW", linewidth= ps.line_width+1, linestyle = ps.line_style[0], color = line_colors[0] ) 
    ax.plot(months, ord_temps, label="ORD", linewidth= ps.line_width+1, linestyle = ps.line_style[0], color = line_colors[1])
    ax.plot(months, sfo_temps, label="SFO", linewidth= ps.line_width+1, linestyle = ps.line_style[0], color = line_colors[2])
    ax.plot(months, jfk_temps, label="JFK", linewidth= ps.line_width+1, linestyle = ps.line_style[0], color = line_colors[3])
    ax.plot(months, lax_temps, label="LAX", linewidth= ps.line_width+1, linestyle = ps.line_style[0], color = line_colors[4])
    
    #ax.set_title("2023 Monthly Average Temperatures in Celsius", fontsize = ps.title_font_size)
    ax.set_xlabel("Month", fontsize=ps.axis_font_size+12)
    ax.set_ylabel("Temperature (°C)", fontsize=ps.axis_font_size+12)

    fig.tight_layout()

    if show_legend:        
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.15, 0.625),loc='lower left')   
    if save_figure:
        plt.savefig(f'../Final_Plots/{save_filename}{file_type}') 
    
        return fig
    
# ----------------------------------------------------------------------------------------------------------------------
#   Plot Flight Profile and Airspeed
# ----------------------------------------------------------------------------------------------------------------------  
def plot_flight_profile(dfw_oct, height, width, 
                        save_figure = True,
                        show_legend = False,
                        save_filename = "Flight_Profile_aircraft_speed",
                        file_type = ".png"):

    
    line_colors   = cm.inferno(np.linspace(0,0.9,14))
    # get plotting style 
    ps      = plot_style()   

    fig1, axis_1 = plt.subplots()
    fig1.set_size_inches(height, width)


    fig2, axis_2 = plt.subplots()
    fig2.set_size_inches(height, width)
    
    fig3, axis_3 = plt.subplots()
    fig3.set_size_inches(height, width)
      
    
    for i in range(14): 

        time     = dfw_oct.segments[i].conditions.frames.inertial.time[:,0] / Units.min
        airspeed = dfw_oct.segments[i].conditions.freestream.velocity[:,0] /   Units['mph']
        altitude = dfw_oct.segments[i].conditions.freestream.altitude[:,0]/Units.feet
        Range    = dfw_oct.segments[i].conditions.frames.inertial.aircraft_range[:,0]/ Units.mi

        segment_tag  =  dfw_oct.segments[i].tag
        segment_name = segment_tag.replace('_F_1_D0', ' ').replace('Charge_Day', 'Charging Segment').replace('_', ' ').replace('CLimb', 'Climb')


        axis_1.plot(time, altitude, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width, label = segment_name)
        axis_1.set_ylabel(r'Altitude (ft)', fontsize = ps.axis_font_size)
        axis_1.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
        set_axes(axis_1)    


        axis_2.plot(time, airspeed, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width )
        axis_2.set_ylabel(r'Airspeed (mph)',fontsize  = ps.axis_font_size)
        axis_2.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
        set_axes(axis_2)
        
        axis_3.plot(time, Range, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width )
        axis_3.set_xlabel('Time (mins)', fontsize = ps.axis_font_size)
        axis_3.set_ylabel(r'Range (mi)', fontsize = ps.axis_font_size)
        set_axes(axis_3)         


    if show_legend:        
        leg =  fig.legend(bbox_to_anchor=(0.5, 0.735), loc='lower center', ncol = 2, fontsize = ps.legend_fontsize) 
        leg.set_title('Flight Segment', prop={'size': ps.legend_title_font_size, 'weight': 'heavy'})    

    # Adjusting the sub-plots for legend 
    #fig.subplots_adjust(top=0.8)

    # set title of plot 
    #title_text    = 'Flight Profile for eCTOL'      
    #fig.suptitle(title_text, fontsize=ps.title_font_size, fontweight='bold')    
    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
   

    if save_figure:
        fig1.savefig(f'../Final_Plots/{save_filename}_Altitude{file_type}') 
        fig2.savefig(f'../Final_Plots/{save_filename}_Airspeed{file_type}') 
        fig3.savefig(f'../Final_Plots/{save_filename}_Range{file_type}') 

   
# ----------------------------------------------------------------------------------------------------------------------
#  Plot Battery Pack Conditions
# ---------------------------------------------------------------------------------------------------------------------- 
def plot_battery_cell_conditions_1(dfw_oct, height, width, 
                        save_figure = True,
                        show_legend = False,
                        save_filename = "Battery_Cell_Conditions",
                        file_type = ".png"): 
    # get plotting style 
    ps      = plot_style()   

    line_colors   = cm.inferno(np.linspace(0,0.9,14))
    
    fig1, axis_1 = plt.subplots()
    fig1.set_size_inches(height, width)

    fig2, axis_2 = plt.subplots()
    fig2.set_size_inches(height, width)
    
    fig3, axis_3 = plt.subplots()
    fig3.set_size_inches(height, width)
    
    fig4, axis_4 = plt.subplots()
    fig4.set_size_inches(height, width)
    
    fig5, axis_5 = plt.subplots()
    fig5.set_size_inches(height, width)
    
    fig6, axis_6 = plt.subplots()
    fig6.set_size_inches(height, width)
    
    fig7, axis_7 = plt.subplots()
    fig7.set_size_inches(height, width)
    
    

  
    b_i = 0 
    for network in dfw_oct.segments[0].analyses.energy.networks: 
        busses  = network.busses
        for bus in busses: 
            for battery in bus.batteries:     

                for i in range(14):  
                    time    = dfw_oct.segments[i].conditions.frames.inertial.time[:,0] / Units.min    
                    battery_conditions  = dfw_oct.segments[i].conditions.energy[bus.tag][battery.tag]    
                    cell_power          = battery_conditions.cell.power[:,0]
                    cell_energy         = battery_conditions.cell.energy[:,0]
                    cell_volts          = battery_conditions.cell.voltage_under_load[:,0]
                    cell_current        = battery_conditions.cell.current[:,0]
                    cell_SOC            = battery_conditions.cell.state_of_charge[:,0]   
                    cell_temperature    = battery_conditions.cell.temperature[:,0] -273
                    cell_heat_generated =battery_conditions.pack.heat_energy_generated[:, 0] / 1000

                            

                    #segment_tag  = dfw_oct.segments[i].tag
                    #segment_name = segment_tag.replace('_', ' ')                     
   
                    axis_1.plot(time, cell_SOC, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    axis_1.set_ylabel(r'SOC',fontsize = ps.axis_font_size)
                    axis_1.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    axis_1.set_ylim([0,1.1])
                    set_axes(axis_1)     

                    axis_2.plot(time, cell_energy/Units.Wh,  linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    axis_2.set_ylabel(r'Energy (W-hr)',fontsize = ps.axis_font_size)
                    axis_2.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    set_axes(axis_2) 
             
                    axis_3.plot(time, cell_current,  linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    axis_3.set_ylabel(r'Current (A)',fontsize = ps.axis_font_size)
                    axis_3.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    set_axes(axis_3)  
             
                    axis_4.plot(time, cell_power, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    axis_4.set_ylabel(r'Power (W)',fontsize = ps.axis_font_size)
                    axis_4.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    set_axes(axis_4)     
                     
                    axis_5.plot(time, cell_volts, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    axis_5.set_ylabel(r'Voltage (V)',fontsize = ps.axis_font_size)
                    axis_5.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    set_axes(axis_5) 
             
                    axis_6.plot(time, cell_temperature, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    axis_6.set_ylabel(r'Temperature, ($\degree$C)',fontsize = ps.axis_font_size)
                    axis_6.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    set_axes(axis_6)
                    
                    axis_7.plot(time, cell_heat_generated, linestyle = ps.line_style[0], markersize= ps.marker_size,  color = line_colors[i], marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
                    axis_7.set_ylabel(r'Heat Generated (KW) ',fontsize = ps.axis_font_size)
                    axis_7.set_xlabel(r'Time (mins)',fontsize = ps.axis_font_size)
                    set_axes(axis_7)
                    
                      

                b_i += 1       

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()
    fig6.tight_layout()
    fig7.tight_layout()
    

    if show_legend:  
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.71, 0.86),loc='lower left')     

    if save_figure:
        fig1.savefig(f'../Final_Plots/{save_filename}_SOC{file_type}') 
        fig2.savefig(f'../Final_Plots/{save_filename}_Energy{file_type}') 
        fig3.savefig(f'../Final_Plots/{save_filename}_Current{file_type}') 
        fig4.savefig(f'../Final_Plots/{save_filename}_Power{file_type}') 
        fig5.savefig(f'../Final_Plots/{save_filename}_Voltage{file_type}') 
        fig6.savefig(f'../Final_Plots/{save_filename}_Temperature{file_type}') 
        fig7.savefig(f'../Final_Plots/{save_filename}_HeatGenerated{file_type}')
    
    

# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results  


if __name__ == '__main__': 
    main()    
    plt.show()
