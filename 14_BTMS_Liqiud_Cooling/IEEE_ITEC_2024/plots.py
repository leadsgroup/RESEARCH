from  RCAIDE.Core import Data, Units  

# Python Imports  
import numpy as np   
import pickle 
import matplotlib.pyplot as plt  
import matplotlib.cm as cm  

import sys 
 
# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():  

    save_figure   = True
    file_type     = '.png'
    fig_name      = 'Figure_Name'
    
    # mass flow rates 
    m_dot_range   = np.array([1,2,3,4])
    
    # define ploting parameters 
    PP           = define_plotting_parameters(data_range_length=len(m_dot_range)) 
      
    # initialize plots 
    fig = plt.figure(fig_name)    
    fig.set_size_inches(PP.figure_width,PP.figure_height)  
    axis = fig.add_subplot(1,1,1)    
    
    # plot 
    for i in range(m_dot_range):
        axis.plot( x , y,markersize = PP.marker_size,  marker = PP.markers[i], linewidth = PP.line_width , color = PP.colors[i], label = 'mdot = ' + str(m_dot_range[i]))   
   
    # legend 
    axis.legend(loc='upper right', ncol = 2, prop={'size': PP.legend_font_size})  
    
    # tight layout 
    plt.tight_layout()
    
    # save results 
    if save_figure: 
        plt.savefig( fig_name+ file_type)   
         
    
    return  

    
def define_plotting_parameters(data_range_length):
    
    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 28,
                  'xtick.labelsize': 24,
                  'ytick.labelsize': 24,
                  'axes.titlesize': 28}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 2 
    plot_parameters.line_style       = '-' 
    plot_parameters.figure_width     = 8
    plot_parameters.figure_height    = 6
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20 
    plot_parameters.plot_grid        = True   
    plot_parameters.markers          = ['o','v','s','^','p','^','D','X','*']
    plot_parameters.colors           = cm.inferno(np.linspace(0,1,data_range_length))                
    plot_parameters.legend_font      = 20                             # legend_font_size          
    
    return 
 
 
if __name__ == '__main__': 
    main()    
    plt.show()   
 
