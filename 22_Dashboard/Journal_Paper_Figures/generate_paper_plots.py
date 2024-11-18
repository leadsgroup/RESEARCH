'''





'''

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTS 
# --------------------------------------------------------------------------------------------------------------------------------------------------- 

from   RCAIDE.Framework.Core import Units 
from   RCAIDE.load           import load as load_results
from   RCAIDE.save           import save as save_results 
import os 
import numpy                 as np           
import pandas                as pd 
import json                  
from   pyatmos               import coesa76
from   urllib.request        import urlopen
import matplotlib as mpl
import matplotlib.pyplot     as plt 
import matplotlib.cm         as cm
import matplotlib.colors     as mcolors
import plotly.io             as pio 
import plotly.express        as px 
from   plotly.graph_objs     import * 
from   mpl_toolkits.mplot3d  import Axes3D
from   scipy.interpolate     import CubicSpline
from matplotlib.gridspec     import GridSpec

def main():
    file_type   =  'png'
    save_figure =  True
    width       =  5
    height      =  4   
    generate_saf_plots(file_type,save_figure,width, height)
    generate_electrification_plots(file_type,save_figure,width, height)
    generate_hydrogen_plots(file_type,save_figure,width, height)
    return

def generate_saf_plots(file_type,save_figure,width, height): 
    saf_data    =  load_results('saf_data.res')

    show_legend_S = 'Yes'
    
    saf_plot_1_selected_fuels_list       = saf_data.saf_plot_1_selected_fuels_list
    saf_plot_1_selected_feedstock_list   = saf_data.saf_plot_1_selected_feedstock_list
    saf_plot_1_selected_airpots_list     = saf_data.saf_plot_1_selected_airpots_list
    saf_plot_1_percent_adoption_list     = saf_data.saf_plot_1_percent_adoption_list
    saf_plot_1_land_area                 = saf_data.saf_plot_1_land_area
    saf_plot_1_CASM_wo_SAF_Aircraft      = saf_data.saf_plot_1_CASM_wo_SAF_Aircraft
    saf_plot_1_CASM_w_SAF_Aircraft       = saf_data.saf_plot_1_CASM_w_SAF_Aircraft
    saf_plot_2_selected_feedstock_list   = saf_data.saf_plot_2_selected_feedstock_list
    saf_plot_2_selected_airpots_list     = saf_data.saf_plot_2_selected_airpots_list
    saf_plot_2_percent_adoption_list     = saf_data.saf_plot_2_percent_adoption_list
    saf_plot_2_percent_fuel_use_list     = saf_data.saf_plot_2_percent_fuel_use_list
    saf_plot_2_Emissions_w_SAF_Aircraft  = saf_data.saf_plot_2_Emissions_w_SAF_Aircraft
    saf_plot_2_Emissions_wo_SAF_Aircraft = saf_data.saf_plot_2_Emissions_wo_SAF_Aircraft
        
    save_filename_15  = "15_saf_plot_Land_HEFA"
    save_filename_16  = "16_saf_plot_Land_ATJ"
    save_filename_17  = "17_saf_plot_CASM_HEFA"
    save_filename_18  = "18_saf_plot_CASM_ATJ"
    save_filename_19  = "19_saf_plot_Emissions"
    
    # get plotting style 
    ps                = plot_style()  
                      
    parameters        = {'axes.labelsize' : ps.axis_font_size,
                         'xtick.labelsize': ps.axis_font_size,
                         'ytick.labelsize': ps.axis_font_size,
                         'axes.titlesize' : ps.title_font_size}
    plt.rcParams.update(parameters)
      
    fig_15            = plt.figure(save_filename_15)
    fig_16            = plt.figure(save_filename_16)
    fig_17            = plt.figure(save_filename_17)
    fig_18            = plt.figure(save_filename_18)
    #fig_19            = plt.figure(save_filename_19)
    
    fig_15.set_size_inches(width,height)
    fig_16.set_size_inches(width,height)
    fig_17.set_size_inches(width,height)
    fig_18.set_size_inches(width,height)
    #fig_19.set_size_inches(width,height)

    axis_15 = fig_15.add_subplot(1,1,1) 
    axis_16 = fig_16.add_subplot(1,1,1)    
    axis_17 = fig_17.add_subplot(1,1,1)    
    axis_18 = fig_18.add_subplot(1,1,1)   
    #axis_19 = fig_19.add_subplot(1,1,1)  
    #axis_20 = fig_19.add_subplot(1,1,1)      

    axis_15.plot(saf_plot_1_percent_adoption_list, saf_plot_1_land_area[0, 1, :]/1000000, color = 'g', marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"HEFA - upper 10 Airports - Soybean") 
    axis_15.plot(saf_plot_1_percent_adoption_list, saf_plot_1_land_area[1, 1, :]/1000000, color = 'b', marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"HEFA - upper 10 Airports - Canola") 
    axis_15.plot(saf_plot_1_percent_adoption_list, saf_plot_1_land_area[0, 3, :]/1000000, color = 'g', marker = ps.markers[1], linewidth = ps.line_width,markersize = ps.marker_size, label = f"HEFA - upper 50 Airports - Soybean") 
    axis_15.plot(saf_plot_1_percent_adoption_list, saf_plot_1_land_area[1, 3, :]/1000000, color = 'b', marker = ps.markers[1], linewidth = ps.line_width,markersize = ps.marker_size, label = f"HEFA - upper 50 Airports - Canola") 
    axis_15.set_ylabel(r'Land area in millions [M Acres]')
    axis_15.set_xlabel(r'Percent adoption [%]')  
    set_axes(axis_15)    
    
    axis_16.plot(saf_plot_1_percent_adoption_list, saf_plot_1_land_area[4, 1, :]/1000000, color = 'g', marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"ATJ - upper 10 Airports - Corn") 
    axis_16.plot(saf_plot_1_percent_adoption_list, saf_plot_1_land_area[5, 1, :]/1000000, color = 'b', marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"ATJ - upper 10 Airports - Wheat") 
    axis_16.plot(saf_plot_1_percent_adoption_list, saf_plot_1_land_area[4, 3, :]/1000000, color = 'g', marker = ps.markers[1], linewidth = ps.line_width,markersize = ps.marker_size, label = f"ATJ - upper 50 Airports - Corn") 
    axis_16.plot(saf_plot_1_percent_adoption_list, saf_plot_1_land_area[5, 3, :]/1000000, color = 'b', marker = ps.markers[1], linewidth = ps.line_width,markersize = ps.marker_size, label = f"ATJ - upper 50 Airports - Wheat") 
    axis_16.set_ylabel(r'Land area in millions [M Acres]')
    axis_16.set_xlabel(r'Percent adoption [%]')  
    set_axes(axis_16)
    
    index_month = 4
    
    non_zero_indices_1 = saf_plot_1_CASM_w_SAF_Aircraft[0, 1, :, index_month] != 0
    non_zero_indices_2 = saf_plot_1_CASM_w_SAF_Aircraft[1, 1, :, index_month] != 0
    non_zero_indices_3 = saf_plot_1_CASM_w_SAF_Aircraft[0, 3, :, index_month] != 0
    non_zero_indices_4 = saf_plot_1_CASM_w_SAF_Aircraft[1, 3, :, index_month] != 0 
    
    xnew_1 = np.linspace(saf_plot_1_percent_adoption_list[non_zero_indices_1][0], saf_plot_1_percent_adoption_list[non_zero_indices_1][-1], 100)
    cs_1   = CubicSpline(saf_plot_1_percent_adoption_list[non_zero_indices_1], saf_plot_1_CASM_w_SAF_Aircraft[0, 1, :, index_month][non_zero_indices_1]) 
    xnew_2 = np.linspace(saf_plot_1_percent_adoption_list[non_zero_indices_2][0], saf_plot_1_percent_adoption_list[non_zero_indices_2][-1], 100)
    cs_2   = CubicSpline(saf_plot_1_percent_adoption_list[non_zero_indices_2], saf_plot_1_CASM_w_SAF_Aircraft[1, 1, :, index_month][non_zero_indices_2])
    xnew_3 = np.linspace(saf_plot_1_percent_adoption_list[non_zero_indices_3][0], saf_plot_1_percent_adoption_list[non_zero_indices_3][-1], 100)
    cs_3   = CubicSpline(saf_plot_1_percent_adoption_list[non_zero_indices_3], saf_plot_1_CASM_w_SAF_Aircraft[0, 3, :, index_month][non_zero_indices_3])
    xnew_4 = np.linspace(saf_plot_1_percent_adoption_list[non_zero_indices_4][0], saf_plot_1_percent_adoption_list[non_zero_indices_4][-1], 100)
    cs_4   = CubicSpline(saf_plot_1_percent_adoption_list[non_zero_indices_4], saf_plot_1_CASM_w_SAF_Aircraft[1, 3, :, index_month][non_zero_indices_4])    
    
    axis_17.plot(xnew_1, cs_1(xnew_1), color = 'g', marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"HEFA - upper 10 Airports - Soybean") 
    axis_17.plot(xnew_2, cs_2(xnew_2), color = 'b', marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"HEFA - upper 10 Airports - Canola")  
    axis_17.plot(xnew_3, cs_3(xnew_3), color = 'g', marker = ps.markers[1], linewidth = ps.line_width,markersize = ps.marker_size, label = f"HEFA - upper 50 Airports - Soybean") 
    axis_17.plot(xnew_4, cs_4(xnew_4), color = 'b', marker = ps.markers[1], linewidth = ps.line_width,markersize = ps.marker_size, label = f"HEFA - upper 50 Airports - Canola")  
    axis_17.set_ylabel(r'CASM [-]')
    axis_17.set_xlabel(r'Percent adoption [%]')  
    set_axes(axis_17)  
    
    non_zero_indices_5 = saf_plot_1_CASM_w_SAF_Aircraft[4, 1, :, index_month] != 0
    non_zero_indices_6 = saf_plot_1_CASM_w_SAF_Aircraft[5, 1, :, index_month] != 0
    non_zero_indices_7 = saf_plot_1_CASM_w_SAF_Aircraft[4, 3, :, index_month] != 0
    non_zero_indices_8 = saf_plot_1_CASM_w_SAF_Aircraft[5, 3, :, index_month] != 0   
    
    xnew_5 = np.linspace(saf_plot_1_percent_adoption_list[non_zero_indices_5][0], saf_plot_1_percent_adoption_list[non_zero_indices_5][-1], 100)
    cs_5   = CubicSpline(saf_plot_1_percent_adoption_list[non_zero_indices_5], saf_plot_1_CASM_w_SAF_Aircraft[4, 1, :, index_month][non_zero_indices_5]) 
    xnew_6 = np.linspace(saf_plot_1_percent_adoption_list[non_zero_indices_6][0], saf_plot_1_percent_adoption_list[non_zero_indices_6][-1], 100)
    cs_6   = CubicSpline(saf_plot_1_percent_adoption_list[non_zero_indices_6], saf_plot_1_CASM_w_SAF_Aircraft[5, 1, :, index_month][non_zero_indices_6])
    xnew_7 = np.linspace(saf_plot_1_percent_adoption_list[non_zero_indices_7][0], saf_plot_1_percent_adoption_list[non_zero_indices_7][-1], 100)
    cs_7   = CubicSpline(saf_plot_1_percent_adoption_list[non_zero_indices_7], saf_plot_1_CASM_w_SAF_Aircraft[4, 3, :, index_month][non_zero_indices_7])
    xnew_8 = np.linspace(saf_plot_1_percent_adoption_list[non_zero_indices_8][0], saf_plot_1_percent_adoption_list[non_zero_indices_8][-1], 100)
    cs_8   = CubicSpline(saf_plot_1_percent_adoption_list[non_zero_indices_8], saf_plot_1_CASM_w_SAF_Aircraft[5, 3, :, index_month][non_zero_indices_8])    
    
    axis_18.plot(xnew_5, cs_5(xnew_5), color = 'g', marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"ATJ - upper 10 Airports - Corn") 
    axis_18.plot(xnew_6, cs_6(xnew_6), color = 'b', marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"ATJ - upper 10 Airports - Wheat") 
    axis_18.plot(xnew_7, cs_7(xnew_7), color = 'g', marker = ps.markers[1], linewidth = ps.line_width,markersize = ps.marker_size, label = f"ATJ - upper 50 Airports - Corn") 
    axis_18.plot(xnew_8, cs_8(xnew_8), color = 'b', marker = ps.markers[1], linewidth = ps.line_width,markersize = ps.marker_size, label = f"ATJ - upper 50 Airports - Wheat") 
    axis_18.set_ylabel(r'CASM [-]')
    axis_18.set_xlabel(r'Percent adoption [%]')  
    set_axes(axis_18)      
    
    # Emissions

    fig_19 = plt.figure(save_filename_19)
    fig_19.set_size_inches(width, height)

    # Define the grid for two vertically stacked plots with custom height ratios
    gs = GridSpec(2, 1, height_ratios=[1, 1.5], hspace=0)
    axis_19 = fig_19.add_subplot(gs[0, 0])  # Top plot
    axis_20 = fig_19.add_subplot(gs[1, 0])  # Bottom plot

    
    fuel_color = 'grey'
    bar_width  = 5
    saf_colors = ['b', 'g', 'c', 'r']
    
    for i, adoption in enumerate(saf_plot_2_percent_adoption_list):
        if adoption == 0:
            axis_19.bar(adoption, 100, color=fuel_color, edgecolor='black', width=bar_width)
        else:
            saf_portion = adoption / 100  
            fuel_portion = 1 - saf_portion  
            
            for j in range(4):  
                offset = (j - 1.5) * bar_width
                label = 'Fuel' if i == 1 and j == 0 else None
                axis_19.bar(adoption + offset, fuel_portion * 100, color=fuel_color, edgecolor='black', width=bar_width, label=label)

                saf_cumulative = fuel_portion
                blend_fractions = saf_plot_2_percent_fuel_use_list[j] / 100  # Normalized blend fractions
                blend_fractions_old = 0
                saf_adjusted_ratios = np.zeros(4)
   
                for w in range(len(blend_fractions) + 1):
                    if w == 0:
                        saf_adjusted_ratios[w] = (blend_fractions[w] - blend_fractions_old)* saf_portion
                        blend_fractions_old  = blend_fractions[w]
                    elif w == 1 or w == 2:
                        saf_adjusted_ratios[w] = (blend_fractions[w] - blend_fractions_old) * saf_portion
                        blend_fractions_old = blend_fractions[w]
                    else:
                        saf_adjusted_ratios[-1] = (1 - blend_fractions_old) * saf_portion
                
                for k, blend_ratio in enumerate(saf_adjusted_ratios):
                    label = f'SAF {k+1}' if i == 1 and j == 0 else None  
                    axis_19.bar(adoption + offset, blend_ratio * 100, bottom=saf_cumulative * 100,
                                color=saf_colors[k], edgecolor='black', width=bar_width, label=label)
                    saf_cumulative += blend_ratio
            
    axis_19.set_xticks(saf_plot_2_percent_adoption_list)
    axis_19.set_xticklabels(saf_plot_2_percent_adoption_list)
    axis_19.set_xlabel('SAF adoption [%]')    
    axis_19.set_ylabel('Blend [%]')
    axis_19.set_yticklabels([0, 25, 50, 75, 100])
    axis_19.set_yticks([0, 25, 50, 75, 100])
    
    emissions_sum  = np.sum(saf_plot_2_Emissions_w_SAF_Aircraft, axis=3)

    for i, adoption in enumerate(saf_plot_2_percent_adoption_list):
        if adoption == 0:
            axis_20.bar(adoption, emissions_sum[:,i, 0], width=bar_width, color='yellow', edgecolor='black')
        else:
            for j in range(4):  
                offset = (j - 1.5) * bar_width
                emission_value = emissions_sum[:,i, j]
                axis_20.bar(adoption + offset, emission_value, color='yellow', edgecolor='black', width=bar_width)
    
    # Flip the bottom plot upside down
    axis_20.invert_yaxis()
    
    # Formatting bottom plot
    axis_20.set_xticks(saf_plot_2_percent_adoption_list)
    axis_20.set_xticklabels(saf_plot_2_percent_adoption_list)
    axis_20.set_ylabel('Emissions [Mt CO2]')
    axis_20.set_xlabel('Percent adoption [%]')
    
    plt.grid(axis='y')
    
    if show_legend_S == 'Yes':    
        leg15 = axis_15.legend(loc='upper left')
        leg16 = axis_16.legend(loc='upper left')
        leg17 = axis_17.legend(loc='upper left') 
        leg18 = axis_18.legend(loc='upper left') 
        leg19 = axis_19.legend(loc='upper right', bbox_to_anchor=(1.35, 1), borderaxespad=0.)
        
    # Adjusting the sub-plots for legend 
    fig_15.tight_layout()    
    fig_16.tight_layout()    
    fig_17.tight_layout()    
    fig_18.tight_layout()    
    fig_19.tight_layout()
    
    if save_figure:
        fig_15.savefig(save_filename_15)
        fig_16.savefig(save_filename_16)
        fig_17.savefig(save_filename_17)
        fig_18.savefig(save_filename_18)
        fig_19.savefig(save_filename_19)   
        
    return     

def generate_electrification_plots(file_type,save_figure,width, height): 
    
    electric_data            = load_results('electric_data.res')
                             
    show_legend_E            = 'Yes'
    
    aircraft_capacity_E_list = ['Twin Otter', 'ATR 72', 'Embraer 190', 'Boeing 737', 'Boeing 777']
    month_list               = ['January', 'February', 'March', 'April', 'May', 'June', 'July','August', 'September', 'October', 'November', 'December']  
    
    aircraft_capacity_E      = electric_data.aircraft_capacity
    weight_fraction_E        = electric_data.weight_fraction  
    cell_e0_E                = electric_data.cell_e0          
    aircraft_range_E         = electric_data.aircraft_range   
    passenger_volume_E       = electric_data.passenger_volume 
    CASM_electric_E          = electric_data.CASM_electric    
    CASM_Jet_A_E             = electric_data.CASM_Jet_A       
                             
    save_filename_1          = "1_electrification_plot_Range_Twin_Otter"
    save_filename_2          = "2_electrification_plot_Range_ATR_72"
    save_filename_3          = "3_electrification_plot_Range_Embraer_190"
    save_filename_4          = "4_electrification_plot_Range_Boeing_737"
    save_filename_5          = "5_electrification_plot_Range_Boeing_777"
    save_filename_6          = "6_electrification_plot_Pax_volume_weight_fr" 
    save_filename_7          = "7_electrification_plot_Pax_volume_spec_en"   
    save_filename_8          = "8_electrification_plot_CASM_Twin_Otter" 
    save_filename_9          = "9_electrification_plot_CASM_ATR_72" 
    save_filename_10         = "10_electrification_plot_CASM_Embraer_190" 
    save_filename_11         = "11_electrification_plot_CASM_Boeing_737" 
    save_filename_12         = "12_electrification_plot_CASM_Boeing_777"  
    
    # get plotting style 
    ps                       = plot_style()  
                             
    parameters               = {'axes.labelsize' : ps.axis_font_size,
                                'xtick.labelsize': ps.axis_font_size,
                                'ytick.labelsize': ps.axis_font_size,
                                'axes.titlesize' : ps.title_font_size}
    plt.rcParams.update(parameters)
      
    fig_1                    = plt.figure(save_filename_1)
    fig_2                    = plt.figure(save_filename_2)
    fig_3                    = plt.figure(save_filename_3)
    fig_4                    = plt.figure(save_filename_4)
    fig_5                    = plt.figure(save_filename_5)
    fig_6                    = plt.figure(save_filename_6)
    fig_7                    = plt.figure(save_filename_7)
    fig_8                    = plt.figure(save_filename_8)
    fig_9                    = plt.figure(save_filename_9)
    fig_10                   = plt.figure(save_filename_10)
    fig_11                   = plt.figure(save_filename_11)
    fig_12                   = plt.figure(save_filename_12)
    
    fig_1.set_size_inches(width,height)
    fig_2.set_size_inches(width,height)
    fig_3.set_size_inches(width,height)
    fig_4.set_size_inches(width,height)
    fig_5.set_size_inches(width,height)
    fig_6.set_size_inches(width,height)
    fig_7.set_size_inches(width,height)
    fig_8.set_size_inches(width,height)
    fig_9.set_size_inches(width,height)
    fig_10.set_size_inches(width,height)
    fig_11.set_size_inches(width,height)
    fig_12.set_size_inches(width,height)

    axis_1                   = fig_1.add_subplot(1,1,1) 
    axis_2                   = fig_2.add_subplot(1,1,1)    
    axis_3                   = fig_3.add_subplot(1,1,1)    
    axis_4                   = fig_4.add_subplot(1,1,1)    
    axis_5                   = fig_5.add_subplot(1,1,1)    
    axis_6                   = fig_6.add_subplot(1,1,1)    
    axis_7                   = fig_7.add_subplot(1,1,1)    
    axis_8                   = fig_8.add_subplot(1,1,1)    
    axis_9                   = fig_9.add_subplot(1,1,1)    
    axis_10                  = fig_10.add_subplot(1,1,1)    
    axis_11                  = fig_11.add_subplot(1,1,1)    
    axis_12                  = fig_12.add_subplot(1,1,1)
    
    # ------------------------------------------------
         
    cell_e0_3d, weight_fraction_3d = np.meshgrid(cell_e0_E, weight_fraction_E)
    
    sum_over_month_E         = np.sum(aircraft_range_E,   axis=1)
    passenger_volume_sum     = np.sum(passenger_volume_E, axis=1)
    
    vmin_1 = np.min(sum_over_month_E[0,:,:]/1000)
    vmax_1 = np.max(sum_over_month_E[0,:,:]/1000) 
    vmin_2 = np.min(sum_over_month_E[1,:,:]/1000)
    vmax_2 = np.max(sum_over_month_E[1,:,:]/1000) 
    vmin_3 = np.min(sum_over_month_E[2,:,:]/1000)
    vmax_3 = np.max(sum_over_month_E[2,:,:]/1000) 
    vmin_4 = np.min(sum_over_month_E[3,:,:]/1000)
    vmax_4 = np.max(sum_over_month_E[3,:,:]/1000) 
    vmin_5 = np.min(sum_over_month_E[4,:,:]/1000)
    vmax_5 = np.max(sum_over_month_E[4,:,:]/1000)     
    
    # Twin Otter 
    contour = axis_1.contourf(cell_e0_3d, weight_fraction_3d, sum_over_month_E[0,:,:]/1000, cmap='viridis', vmin=vmin_1, vmax=vmax_1)
    axis_1.set_xlabel(r'Specific energy density [Wh/kg]')
    axis_1.set_ylabel(r'Battery weight fraction [%]')
    set_axes(axis_1)
    cbar = fig_1.colorbar(contour, ax=axis_1, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label(r'Range in thousand miles [K mi]')  

    specific_energy_density = []
    weight_fraction = []
    for j in range(passenger_volume_sum.shape[1]):  
        for k in range(passenger_volume_sum.shape[2]): 
            if passenger_volume_sum[0, j, k] < passenger_volume_sum[0, j - 1, k]:  
                specific_energy_density.append(cell_e0_E[k])  
                weight_fraction.append(weight_fraction_E[j])
                break
    
    if len(specific_energy_density) == len(weight_fraction):
        axis_1.plot(specific_energy_density, weight_fraction, color='white', linestyle = '-', label="Passenger Reduction")
        axis_1.legend()    
        
    # ATR 72
    contour = axis_2.contourf(cell_e0_3d, weight_fraction_3d, sum_over_month_E[1,:,:]/1000, cmap='viridis', vmin=vmin_2, vmax=vmax_2)
    axis_2.set_xlabel(r'Specific energy density [Wh/kg]')
    axis_2.set_ylabel(r'Battery weight fraction [%]')
    set_axes(axis_2)  
    cbar = fig_2.colorbar(contour, ax=axis_2, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label(r'Range in thousand miles [K mi]')  
    
    specific_energy_density = []
    weight_fraction = []
    for j in range(passenger_volume_sum.shape[1]):  
        for k in range(passenger_volume_sum.shape[2]): 
            if passenger_volume_sum[1, j, k] < passenger_volume_sum[1, j - 1, k]:  
                specific_energy_density.append(cell_e0_E[k])  
                weight_fraction.append(weight_fraction_E[j])
                break
    
    if len(specific_energy_density) == len(weight_fraction):
        axis_2.plot(specific_energy_density, weight_fraction, color='white', linestyle = '-', label="Passenger Reduction")
        axis_2.legend()        
    
    # Embraer 190
    contour = axis_3.contourf(cell_e0_3d, weight_fraction_3d, sum_over_month_E[2,:,:]/1000, cmap='viridis', vmin=vmin_3, vmax=vmax_3)
    axis_3.set_xlabel(r'Specific energy density [Wh/kg]')
    axis_3.set_ylabel(r'Battery weight fraction [%]')
    set_axes(axis_3) 
    cbar = fig_3.colorbar(contour, ax=axis_3, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label(r'Range in thousand miles [K mi]')  
    
    specific_energy_density = []
    weight_fraction = []
    for j in range(passenger_volume_sum.shape[1]):  
        for k in range(passenger_volume_sum.shape[2]): 
            if passenger_volume_sum[2, j, k] < passenger_volume_sum[2, j - 1, k]:  
                specific_energy_density.append(cell_e0_E[k])  
                weight_fraction.append(weight_fraction_E[j])
                break
    
    if len(specific_energy_density) == len(weight_fraction):
        axis_3.plot(specific_energy_density, weight_fraction, color='white', linestyle = '-', label="Passenger Reduction")
        axis_3.legend()    
    
    # Boeing 737
    contour = axis_4.contourf(cell_e0_3d, weight_fraction_3d, sum_over_month_E[3,:,:]/1000, cmap='viridis', vmin=vmin_4, vmax=vmax_4)
    axis_4.set_xlabel(r'Specific energy density [Wh/kg]')
    axis_4.set_ylabel(r'Battery weight fraction [%]')
    set_axes(axis_4) 
    cbar = fig_4.colorbar(contour, ax=axis_4, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label(r'Range in thousand miles [K mi]')   
    
    specific_energy_density = []
    weight_fraction = []
    for j in range(passenger_volume_sum.shape[1]):  
        for k in range(passenger_volume_sum.shape[2]): 
            if passenger_volume_sum[3, j, k] < passenger_volume_sum[3, j - 1, k]:  
                specific_energy_density.append(cell_e0_E[k])  
                weight_fraction.append(weight_fraction_E[j])
                break
    
    if len(specific_energy_density) == len(weight_fraction):
        axis_4.plot(specific_energy_density, weight_fraction, color='white', linestyle = '-', label="Passenger Reduction")
        axis_4.legend()     
    
    # Boeing 777
    contour = axis_5.contourf(cell_e0_3d, weight_fraction_3d, sum_over_month_E[4,:,:]/1000, cmap='viridis', vmin=vmin_5, vmax=vmax_5)
    axis_5.set_xlabel(r'Specific energy density [Wh/kg]')
    axis_5.set_ylabel(r'Battery weight fraction [%]')
    set_axes(axis_5)
    cbar = fig_5.colorbar(contour, ax=axis_5, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label(r'Range in thousand miles [K mi]') 
    
    specific_energy_density = []
    weight_fraction = []
    for j in range(passenger_volume_sum.shape[1]):  
        for k in range(passenger_volume_sum.shape[2]): 
            if passenger_volume_sum[4, j, k] < passenger_volume_sum[4, j - 1, k]:  
                specific_energy_density.append(cell_e0_E[k])  
                weight_fraction.append(weight_fraction_E[j])
                break
    
    if len(specific_energy_density) == len(weight_fraction):
        axis_5.plot(specific_energy_density, weight_fraction, color='white', linestyle = '-', label="Passenger Reduction")
        axis_5.legend()     
    
    # ------------------------------------------------
    
    sum_over_month_E         = np.sum(passenger_volume_E, axis=1)   
    line_colors_capacity_E   = cm.inferno(np.linspace(0,0.9,len(aircraft_capacity_E)))     
    index                    = 3
    xnew                     = np.linspace(weight_fraction_E[0], weight_fraction_E[-1], 100)
        
    for i in range(len(aircraft_capacity_E)):
        cs                            = CubicSpline(weight_fraction_E, sum_over_month_E [i,:,index])
        ynew                          = cs(xnew)        
        axis_6.plot(xnew, ynew/1000000, color = line_colors_capacity_E[i], marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"{aircraft_capacity_E_list[i]}") 
    axis_6.set_ylabel(r'Passenger volume in millions [M]')
    axis_6.set_xlabel(r'Battery weight fraction [%]')  
    axis_6.set_title(f"Specific energy density: {cell_e0_E[index]} Wh/kg")    
    set_axes(axis_6) 
    
    # ------------------------------------------------
    
    index                    = 3
    xnew                     = np.linspace(cell_e0_E[0], cell_e0_E[-1], 100)
        
    for i in range(len(aircraft_capacity_E)):
        cs                            = CubicSpline(cell_e0_E, sum_over_month_E[i,index,:])
        ynew                          = cs(xnew)                
        axis_7.plot(xnew, ynew/1000000, color = line_colors_capacity_E[i], marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"{aircraft_capacity_E_list[i]}") 
    axis_7.set_xlabel(r'Specific energy density [Wh/kg]')
    axis_7.set_ylabel(r'Passenger volume in millions [M]')
    axis_7.set_title(f"Weight fraction: {weight_fraction_E[index]} %")     
    set_axes(axis_7) 
    
    # ------------------------------------------------
    
    index_month             = 4
    line_colors_capacity1   = cm.inferno(np.linspace(0,0.9,len(weight_fraction_E)))     
    sm                      = mpl.cm.ScalarMappable(cmap=cm.inferno, norm=plt.Normalize(vmin=min(weight_fraction_E), vmax=max(weight_fraction_E)))
    sm.set_array([])        
    
    # Twin Otter 
    axis = axis_8 
    for i in range(len(weight_fraction_E)):
        non_zero_indices = CASM_electric_E[0, index_month, i, :] != 0
        if np.sum(non_zero_indices) >= 2:
            xnew = np.linspace(cell_e0_E[non_zero_indices][0], cell_e0_E[non_zero_indices][-1], 100)
            cs   = CubicSpline(cell_e0_E[non_zero_indices], CASM_electric_E[0, index_month, i, :][non_zero_indices])
            axis.plot(xnew, cs(xnew), 
                        color=line_colors_capacity1[i], 
                        marker=ps.markers[0], 
                        linewidth=ps.line_width, 
                        markersize=ps.marker_size, 
                        label=f"{weight_fraction_E[i]}%") 
    axis.set_xlabel(r'Specific energy density [Wh/kg]')
    axis.set_ylabel(r'CASM [-]')  
    axis.set_title(f"Month: {month_list[index_month]}")
    set_axes(axis) 
    cbar = fig_8.colorbar(sm, ax=axis_8, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label('Weight Fraction [%]')      
    
    # ATR 72
    axis = axis_9 
    for i in range(len(weight_fraction_E)):
        non_zero_indices = CASM_electric_E[1, index_month, i, :] != 0
        if np.sum(non_zero_indices) >= 2:
            xnew = np.linspace(cell_e0_E[non_zero_indices][0], cell_e0_E[non_zero_indices][-1], 100)
            cs   = CubicSpline(cell_e0_E[non_zero_indices], CASM_electric_E[1, index_month, i, :][non_zero_indices])
            axis.plot(xnew, cs(xnew), 
                        color=line_colors_capacity1[i], 
                        marker=ps.markers[0], 
                        linewidth=ps.line_width, 
                        markersize=ps.marker_size, 
                        label=f"{weight_fraction_E[i]}%") 
    axis.set_xlabel(r'Specific energy density [Wh/kg]')
    axis.set_ylabel(r'CASM [-]')  
    axis.set_title(f"Month: {month_list[index_month]}")
    set_axes(axis) 
    cbar = fig_9.colorbar(sm, ax=axis_9, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label('Weight Fraction [%]') 
    
    # Embraer 190
    axis = axis_10 
    for i in range(len(weight_fraction_E)):
        non_zero_indices = CASM_electric_E[2, index_month, i, :] != 0
        if np.sum(non_zero_indices) >= 2:
            xnew = np.linspace(cell_e0_E[non_zero_indices][0], cell_e0_E[non_zero_indices][-1], 100)
            cs   = CubicSpline(cell_e0_E[non_zero_indices], CASM_electric_E[2, index_month, i, :][non_zero_indices])
            axis.plot(xnew, cs(xnew), 
                        color=line_colors_capacity1[i], 
                        marker=ps.markers[0], 
                        linewidth=ps.line_width, 
                        markersize=ps.marker_size, 
                        label=f"{weight_fraction_E[i]}%") 
    axis.set_xlabel(r'Specific energy density [Wh/kg]')
    axis.set_ylabel(r'CASM [-]')  
    axis.set_title(f"Month: {month_list[index_month]}")
    set_axes(axis) 
    cbar = fig_10.colorbar(sm, ax=axis_10, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label('Weight Fraction [%]') 
    
    # Boeing 737
    axis = axis_11 
    for i in range(len(weight_fraction_E)):
        non_zero_indices = CASM_electric_E[3, index_month, i, :] != 0
        if np.sum(non_zero_indices) >= 2:
            xnew = np.linspace(cell_e0_E[non_zero_indices][0], cell_e0_E[non_zero_indices][-1], 100)
            cs   = CubicSpline(cell_e0_E[non_zero_indices], CASM_electric_E[3, index_month, i, :][non_zero_indices])
            axis.plot(xnew, cs(xnew), 
                        color=line_colors_capacity1[i], 
                        marker=ps.markers[0], 
                        linewidth=ps.line_width, 
                        markersize=ps.marker_size, 
                        label=f"{weight_fraction_E[i]}%") 
    axis.set_xlabel(r'Specific energy density [Wh/kg]')
    axis.set_ylabel(r'CASM [-]')  
    axis.set_title(f"Month: {month_list[index_month]}")
    set_axes(axis) 
    cbar = fig_11.colorbar(sm, ax=axis_11, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label('Weight Fraction [%]')    
    
    # Boeing 777
    axis = axis_12 
    for i in range(len(weight_fraction_E)):
        non_zero_indices = CASM_electric_E[4, index_month, i, :] != 0
        if np.sum(non_zero_indices) >= 2:
            xnew = np.linspace(cell_e0_E[non_zero_indices][0], cell_e0_E[non_zero_indices][-1], 100)
            cs   = CubicSpline(cell_e0_E[non_zero_indices], CASM_electric_E[4, index_month, i, :][non_zero_indices])
            axis.plot(xnew, cs(xnew), 
                        color=line_colors_capacity1[i], 
                        marker=ps.markers[0], 
                        linewidth=ps.line_width, 
                        markersize=ps.marker_size, 
                        label=f"{weight_fraction_E[i]}%") 
    axis.set_xlabel(r'Specific energy density [Wh/kg]')
    axis.set_ylabel(r'CASM [-]')  
    axis.set_title(f"Month: {month_list[index_month]}")
    set_axes(axis) 
    cbar = fig_12.colorbar(sm, ax=axis_12, orientation='vertical', fraction=0.02, pad=0.15)
    cbar.set_label('Weight Fraction [%]')         
    
    # ------------------------------------------------
         
    if show_legend_E == 'Yes':    

        leg6 =  axis_6.legend(loc='upper right')
        leg7 =  axis_7.legend(loc='upper left') 
        
    # Adjusting the sub-plots for legend 
    fig_1.tight_layout()    
    fig_2.tight_layout()    
    fig_3.tight_layout()    
    fig_4.tight_layout()     
    fig_5.tight_layout()     
    fig_6.tight_layout()     
    fig_7.tight_layout()      
    fig_8.tight_layout()      
    fig_9.tight_layout()      
    fig_10.tight_layout()      
    fig_11.tight_layout()      
    fig_12.tight_layout()  
    
    if save_figure:
        fig_1.savefig(save_filename_1)
        fig_2.savefig(save_filename_2)
        fig_3.savefig(save_filename_3)
        fig_4.savefig(save_filename_4)
        fig_5.savefig(save_filename_5)
        fig_6.savefig(save_filename_6)
        fig_7.savefig(save_filename_7) 
        fig_8.savefig(save_filename_8) 
        fig_9.savefig(save_filename_9) 
        fig_10.savefig(save_filename_10) 
        fig_11.savefig(save_filename_11) 
        fig_12.savefig(save_filename_12)        
    return    


def generate_hydrogen_plots(file_type,save_figure,width, height): 
    hydrogen_data         =  load_results('hydrogen_data.res')   
    
    aircraft_capacity_H_list = ['Twin Otter', 'ATR 72', 'Embraer 190', 'Boeing 737', 'Boeing 777']    
                          
    show_legend_H         = 'Yes'
                          
    aircraft_capacity_H   = hydrogen_data.aircraft_capacity
    volume_fraction_H     = hydrogen_data.volume_fraction  
    CO2e_w_H2_Aircraft_H  = hydrogen_data.CO2e_w_H2_Aircraft
    CO2e_wo_H2_Aircraft_H = hydrogen_data.CO2e_wo_H2_Aircraft
    aircraft_range_H      = hydrogen_data.aircraft_range   
    passenger_volume_H    = hydrogen_data.passenger_volume 
    CASM_w_H2_Aircraft_H  = hydrogen_data.CASM_w_H2_Aircraft    
    CASM_wo_H2_Aircraft_H = hydrogen_data.CASM_wo_H2_Aircraft       
        
    save_filename_13   = "13_hydrogen_plot_Pax_VolumeFr"
    save_filename_14   = "14_hydrogen_plot_CASM_VolumeFr"
    
    # get plotting style 
    ps                = plot_style()  
                      
    parameters        = {'axes.labelsize' : ps.axis_font_size,
                         'xtick.labelsize': ps.axis_font_size,
                         'ytick.labelsize': ps.axis_font_size,
                         'axes.titlesize' : ps.title_font_size}
    plt.rcParams.update(parameters)
      
    fig_13             = plt.figure(save_filename_13)
    fig_14             = plt.figure(save_filename_14)
    
    fig_13.set_size_inches(width,height)
    fig_14.set_size_inches(width,height)

    axis_13 = fig_13.add_subplot(1,1,1) 
    axis_14 = fig_14.add_subplot(1,1,1)  
    
    # ------------------------------------------------
    
    line_colors_capacity   = cm.inferno(np.linspace(0,0.9,len(aircraft_capacity_H)))
    sum_over_month_H       = np.sum(passenger_volume_H, axis=1)
    
    xnew                   = np.linspace(volume_fraction_H[0], volume_fraction_H[-1], 100)

    for i in range(len(aircraft_capacity_H)):
        cs                            = CubicSpline(volume_fraction_H, sum_over_month_H[i,:])
        ynew                          = cs(xnew)        
        axis_13.plot(xnew, ynew/1000000, color = line_colors_capacity[i], marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"{aircraft_capacity_H_list[i]}") 
    axis_13.set_xlabel(r'Volume fraction [%]')
    axis_13.set_ylabel(r'Passenger volume in millions [M]')
    set_axes(axis_13)    
    
    # ------------------------------------------------
    
    index_month            = 4
    
    line_colors_capacity1  = cm.inferno(np.linspace(0,0.9,len(aircraft_capacity_H)))     
        
    for i in range(len(aircraft_capacity_H)):
        cs                            = CubicSpline(volume_fraction_H, CASM_w_H2_Aircraft_H[i,index_month,:])
        ynew                          = cs(xnew)         
        axis_14.plot(xnew, ynew, color = line_colors_capacity1[i], marker = ps.markers[0], linewidth = ps.line_width,markersize = ps.marker_size, label = f"{aircraft_capacity_H_list[i]}") 
    axis_14.set_xlabel(r'Volume fraction [%]')
    axis_14.set_ylabel(r'CASM [-]')
    set_axes(axis_14)         
    
    if show_legend_H == 'Yes':    
        leg13 =  axis_13.legend(loc='upper left')
        leg14 =  axis_14.legend(loc='upper left')
        
    # Adjusting the sub-plots for legend 
    fig_13.tight_layout()    
    fig_14.tight_layout()    
    
    if save_figure:
        fig_13.savefig(save_filename_13)
        fig_14.savefig(save_filename_14)  
        
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
    parameters = {'axes.labelsize' : 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize' : 18,
                  'figure.dpi'     : 200
                  }
    #parameters = {'axes.labelsize' : 6,
                  #'xtick.labelsize': 4,
                  #'ytick.labelsize': 4,
                  #'axes.titlesize' : 6,
                  #'figure.dpi'     : 200
                  #}    


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
    #plot_parameters.marker_size            = 1
    #plot_parameters.legend_fontsize        = '1'
    #plot_parameters.legend_title_font_size = 4
    #plot_parameters.axis_font_size         = 6
    #plot_parameters.title_font_size        = 6     

    return plot_parameters


if __name__ == '__main__':
    main()
    plt.show()
