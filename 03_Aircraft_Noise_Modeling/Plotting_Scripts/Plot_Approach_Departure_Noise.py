# RCAIDE Imports 
from  RCAIDE.Framework.Core import Data, Units
from RCAIDE.Library.Plots  import *
import plotly.graph_objects as go 
from RCAIDE import  load 
from RCAIDE import  save

# Python Imports  
import numpy as np   
import pickle 
import matplotlib.pyplot as plt  
import matplotlib.cm as cm  
import  os
import  sys

local_path_1 =  os.path.split(sys.path[0])[0]

sys.path.append( os.path.join(local_path_1, 'Approach_Departure_Noise'))
# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():  

    save_figure = True
    show_figure = True
    filetype    = '.png'
    
    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 24,
                  'xtick.labelsize': 22,
                  'ytick.labelsize': 22,
                  'axes.titlesize': 24}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 3 
    plot_parameters.line_style       = '-' 
    plot_parameters.figure_width     = 7 
    plot_parameters.figure_height    = 5 
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20 
    plot_parameters.plot_grid        = True   
    plot_parameters.markers          = ['o','v','s','^','p','^','D','X','*']
    plot_parameters.colors           = cm.inferno(np.linspace(0,1,5))     
    plot_parameters.lw               = 2                              # line_width                
    plot_parameters.legend_font      = 20                             # legend_font_size          
    
  
    aircraft_tags  = ['HC','TR','SR']       

    filename_list = ['HC_Processed_Noise', 'TR_Processed_Noise', 'TSR_Processed_Noise'] 
    for i, postprocessed_filename in enumerate(filename_list):
        
        # --------------------------------------------------------------------------------------    
        # CONTOUR PLOT
        # -------------------------------------------------------------------------------------- 
        separator       = os.path.sep
        relative_path   = sys.path[-1] + separator       
        
        results = load(relative_path + postprocessed_filename + '.res')
        
        fig_name = aircraft_tags[i] + '_Approach_Contour'
        fig      = plt.figure(fig_name)    
        fig.set_size_inches(plot_parameters.figure_width,plot_parameters.figure_height)  
        axis     = fig.add_subplot(1,1,1)       
            
        SPL = results.L_max           
        Y   = results.microphone_locations[:, 1].reshape(300,300)
        X   = results.microphone_locations[:, 0].reshape(300,300)
          
        levs                = np.linspace(30,100,15)   
        CS                  = axis.contourf(X/Units.feet , Y/Units.feet, SPL, levels  = levs, cmap=plt.cm.jet)
        
        levs2 = np.array([45,65])
        CS2 = axis.contour(X/Units.feet , Y/Units.feet, SPL, levels=levs2, linewidths=1, colors='k') 
            
            
        axis.clabel(CS2, fontsize=20)
        cbar                = fig.colorbar(CS)
        cbar.ax.set_ylabel('L$_{Amax}$ [dBA]', rotation =  90)      
        axis.set_xlabel('x [ft]')
        axis.set_ylabel('y [ft]',labelpad = 15) 
        
        axis.grid(False)  
        axis.minorticks_on()     
        plt.tight_layout() 

         
        # calculate position and direction vectors:
        x_pos = results.trajectory[:,0] / Units.feet
        y_pos = results.trajectory[:,1] / Units.feet
        
        x0_pos = x_pos[:-1] 
        x1_pos = x_pos[1:] 
        y0_pos = y_pos[:-1] 
        y1_pos = y_pos[1:] 
        xpos = (x0_pos+x1_pos)/2
        ypos = (y0_pos+y1_pos)/2
        xdir = x1_pos-x0_pos
        ydir = y1_pos-y0_pos  
        axis.plot(x_pos,y_pos ,  color = 'w', linewidth = 1)
        axis.set_xlim([-7500, 7500])
        axis.set_ylim([-7500, 7500])
        
        # plot arrow on each line:
        arrow_i =  1
        for X_pos,Y_pos,dX_pos,dY_pos in zip(xpos, ypos, xdir, ydir):
            if (arrow_i % 50) ==  0:
                axis.arrow(X_pos,Y_pos, dX_pos,dY_pos, head_width=500, head_length=500, fc='w', ec='w')
            arrow_i += 1
 
        if save_figure:
            safe_filename = fig_name
            plt.savefig(safe_filename + filetype, dpi=600)  

def colorax(vmin, vmax):
    return dict(cmin=vmin,  cmax=vmax)       
 
 

if __name__ == '__main__': 
    main()    
    plt.show()   
 
