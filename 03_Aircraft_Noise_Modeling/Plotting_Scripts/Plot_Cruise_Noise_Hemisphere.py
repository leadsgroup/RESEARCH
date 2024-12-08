# RCAIDE Imports 
from  RCAIDE.Framework.Core import Data, Units
from RCAIDE.Library.Plots  import *
import plotly.graph_objects as go 
from RCAIDE import  load 
from RCAIDE import  save

# Python Imports  
import numpy as np
import  pandas as pd
import pickle 
import matplotlib.pyplot as plt  
import matplotlib.cm as cm  
import  os
import  sys

local_path_1 =  os.path.split(sys.path[0])[0]

sys.path.append( os.path.join(local_path_1, 'Hemisphere_Noise'))



local_path_1 =  os.path.split(sys.path[0])[0]

local_path_2 =  os.path.split(os.path.split(sys.path[0])[0])[0]   
sys.path.append(os.path.join(local_path_2, 'Aircraft' + os.path.sep + 'Hexacopter'))
sys.path.append(os.path.join(local_path_2, 'Aircraft' + os.path.sep + 'Tiltrotor'))
sys.path.append(os.path.join(local_path_2, 'Aircraft' + os.path.sep + 'Tilt_Stopped_Rotor'))
from Hexacopter                       import vehicle_setup as  HC_vehicle_setup
from Hexacopter                       import configs_setup as  HC_configs_setup 
from Tiltrotor                        import vehicle_setup as  TR_vehicle_setup
from Tiltrotor                        import configs_setup as  TR_configs_setup 
from Tilt_Stopped_Rotor_V_Tail        import vehicle_setup as  TSR_vehicle_setup
from Tilt_Stopped_Rotor_V_Tail        import configs_setup as  TSR_configs_setup

# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():
        

    save_figure = True
    show_figure = True
    plot_aircraft = False 
    filetype    = '.png'
    
    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 28,
                  'xtick.labelsize': 24,
                  'ytick.labelsize': 24,
                  'axes.titlesize': 28}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 3 
    plot_parameters.line_style       = '-' 
    plot_parameters.figure_width     = 10 
    plot_parameters.figure_height    = 7 
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20 
    plot_parameters.plot_grid        = True   
    plot_parameters.markers          = ['o','v','s','^','p','^','D','X','*']
    plot_parameters.colors           = cm.inferno(np.linspace(0,1,5))     
    plot_parameters.lw               = 2                              # line_width                
    plot_parameters.legend_font      = 20                             # legend_font_size          
    
  
    aircraft_tags  = ['HC' ,'TR','TSR']       

    filename_list = ['HC_Hemisphere_Processed_Noise', 'TR_Hemisphere_Processed_Noise', 'TSR_Hemisphere_Processed_Noise'] 
    for i in range(len(aircraft_tags)):
        
        # --------------------------------------------------------------------------------------    
        # CONTOUR PLOT
        # -------------------------------------------------------------------------------------- 
        separator       = os.path.sep
        relative_path   = sys.path[-4] + separator       
        
        results = load(relative_path + filename_list[i] + '.res')  
        

        r     = 20     
        X     = r * np.outer(np.cos(results.theta),np.sin(results.phi))
        Y     = r * np.outer(np.sin(results.theta),np.sin(results.phi))
        Z     = r * np.outer(np.cos(results.phi), np.ones(np.size(results.theta))).T
         
        hemisphere_noise         = results.hemisphere_noise.reshape( (len(results.theta),len(results.phi))) 
        colormap                 = 'jet'
        file_type                = ".png"
        background_color         = 'white'
        grid_color               = 'white'
        noise_scale_label        = 'dBA' 
        width                    = 1400
        height                   = 800
        hemisphere_filename      = aircraft_tags[i] + '_Hemisphere_Contour'

        # ---------------------------------------------------------------------------
        # TRHEE DIMENSIONAL NOISE CONTOUR
        # ---------------------------------------------------------------------------
    
        plot_data   = []
        
        if plot_aircraft: 
            if aircraft_tags[i] == 'HC':
                vec     = HC_vehicle_setup(redesign_rotors = False)     
                configs  = HC_configs_setup(vec)
                vehicle  = configs.forward_flight
            elif aircraft_tags[i] == 'TSR':    
                vec  = TSR_vehicle_setup(redesign_rotors = False)     
                configs  = TSR_configs_setup(vec)
                vehicle  = configs.forward_flight
            elif aircraft_tags[i] == 'TR':
                vec      = TR_vehicle_setup(redesign_rotors = False)     
                configs  = TR_configs_setup(vec)
                vehicle  = configs.forward_flight
     
            plot_data,_,_,_,_,_,_, = generate_3d_vehicle_geometry_data(plot_data,vehicle)        
            
        # TERRAIN CONTOUR     
        hemisphere_contour   = go.Surface(x=X,y=-Y,z=Z, surfacecolor= hemisphere_noise, colorscale  = colormap, cmin=45, cmax=100, ) #colorbar   = dict(title = noise_scale_label, titleside = "right", orientation = "v"))
        plot_data.append(hemisphere_contour) 
        hemisphere_contour   = go.Surface(x=X,y=Y,z=Z, surfacecolor= hemisphere_noise, colorscale  = colormap, showscale   = True, cmin=45, cmax=100,  colorbar   = dict(title = noise_scale_label, titleside = "right", orientation = "v"))        
       
        plot_data.append(hemisphere_contour) 
    
        # Define Colorbar Bounds 
        fig_3d = go.Figure(data=plot_data) 
                             
        fig_3d.update_layout( 
                 title_x                                = 0.5,
                 width                                  = width,
                 height                                 = height,  
                 font=dict(
                     family="Times New Roman",
                     size=28,
                     color="black"
                 ),           
                 scene_aspectmode                       = 'auto', 
                 scene                                  = dict(xaxis = dict(visible=False),
                                                            yaxis = dict(visible=False),
                                                            zaxis =dict(visible=False)), 
                 scene_camera=dict(up    = dict(x=0, y=0, z=1),
                                   center= dict(x=0.4, y=0, z=0.2),
                                   eye   = dict(x=1.0, y=1.0, z=1.))   
        )   
        if save_figure:
            fig_3d.write_image(hemisphere_filename + ".png")
            
        
        if show_figure:
            fig_3d.write_html( hemisphere_filename + '.html', auto_open=True)              
         

if __name__ == '__main__': 
    main()    
    plt.show()   
 
