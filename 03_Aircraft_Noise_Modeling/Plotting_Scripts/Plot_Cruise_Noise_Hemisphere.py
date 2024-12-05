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

sys.path.append( os.path.join(local_path_1, 'Approach_Departure_Noise'))
# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():
    
    df = pd.DataFrame.from_dict({'x' : [0,3,8,7,5,3,2,1],
                                 'y' : [0,1,3,5,9,8,7,5]})
    x = df['x']
    y = df['y']
    # calculate position and direction vectors:
    x0 = x.iloc[range(len(x)-1)].values
    x1 = x.iloc[range(1,len(x))].values
    y0 = y.iloc[range(len(y)-1)].values
    y1 = y.iloc[range(1,len(y))].values
    xpos = (x0+x1)/2
    ypos = (y0+y1)/2
    xdir = x1-x0
    ydir = y1-y0
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.plot(x,y)
    # plot arrow on each line:
    for X,Y,dX,dY in zip(xpos, ypos, xdir, ydir):
        ax.annotate("", xytext=(X,Y),xy=(X+0.001*dX,Y+0.001*dY), 
        arrowprops=dict(arrowstyle="->", color='k'), size = 20)
        

    #save_figure = True
    #show_figure = True
    #filetype    = '.png'
    
    ## Universal Plot Settings 
    #plt.rcParams['axes.linewidth'] = 1.
    #plt.rcParams["font.family"] = "Times New Roman"
    #parameters = {'axes.labelsize': 28,
                  #'xtick.labelsize': 24,
                  #'ytick.labelsize': 24,
                  #'axes.titlesize': 28}
    #plt.rcParams.update(parameters)
    #plot_parameters                  = Data()
    #plot_parameters.line_width       = 3 
    #plot_parameters.line_style       = '-' 
    #plot_parameters.figure_width     = 10 
    #plot_parameters.figure_height    = 7 
    #plot_parameters.marker_size      = 10 
    #plot_parameters.legend_font_size = 20 
    #plot_parameters.plot_grid        = True   
    #plot_parameters.markers          = ['o','v','s','^','p','^','D','X','*']
    #plot_parameters.colors           = cm.inferno(np.linspace(0,1,5))     
    #plot_parameters.lw               = 2                              # line_width                
    #plot_parameters.legend_font      = 20                             # legend_font_size          
    
  
    #aircraft_tags  = ['HC','TR','SR']       

    #filename_list = ['HC_Processed_Noise', 'TR_Processed_Noise', 'TSR_Processed_Noise'] 
    #for i, postprocessed_filename in enumerate(filename_list):
        
        ## --------------------------------------------------------------------------------------    
        ## CONTOUR PLOT
        ## -------------------------------------------------------------------------------------- 
        #separator       = os.path.sep
        #relative_path   = sys.path[-1] + separator       
        
        #results = load(relative_path + postprocessed_filename + '.res')
        
        ### @ingroup Library-Plots-Noise 
        #min_noise_level          = 35  
        #max_noise_level          = 90
        #PHI ,  THETA             = np.meshgrid(results.phi, results.theta)
        
    
        #r     = 20             
        #phi   = results.phi
        #theta = results.theta
     
        #X                        = r * np.outer(np.sin(phi), np.cos(theta))
        #Y                        = r * np.outer(np.sin(phi), np.sin(theta))
        #Z                        = r * np.outer(np.cos(phi), np.ones(np.size(theta))) 
        #hemisphere_noise         = results.hemisphere_noise.reshape(len(results.phi),len(results.theta)) 
        #colormap                 = 'jet'
        #file_type                = ".png"
        #background_color         = 'white'
        #grid_color               = 'white'
        #noise_scale_label        = r'L$_{A}$ [dBA]' 
        #width                    = 1400
        #height                   = 800
        #hemisphere_filename      = aircraft_tags[i] + '_Hemisphere_Contour'

        ## ---------------------------------------------------------------------------
        ## TRHEE DIMENSIONAL NOISE CONTOUR
        ## ---------------------------------------------------------------------------
        
        ##fig = plt.figure(figsize = (8, 6))
        ##ax = fig.add_subplot(projection='3d')
        ###N = 100
        ###r = np.linspace(0, 1, N)
        ###z = np.sqrt(1 - r**2)
        ###intensity = np.linspace(hemisphere_noise, hemisphere_noise, N).reshape(1, -1)
        
        ##surf = ax.plot_surface(X, Y, -Z, facecolors=cm.jet(hemisphere_noise))
        ###ax.axes.set_zlim3d(-1, 1)
        ##plt.show()
        
        
    
        #plot_data   = []      
            
        ## TERRAIN CONTOUR 
        #hemisphere_contour   = contour_surface_slice(X,Y,-Z,hemisphere_noise,color_scale=colormap, showscale= True, colorbar_title = noise_scale_label, colorbar_location = 'right', colorbar_orientation = 'v' )
        #plot_data.append(hemisphere_contour) 
    
        ## Define Colorbar Bounds 
        #fig_3d = go.Figure(data=plot_data)  
            
                             
        #fig_3d.update_layout(
                 #title_text                             = hemisphere_filename, 
                 #title_x                                = 0.5,
                 #width                                  = width,
                 #height                                 = height, 
                 #font_size                              = 12,
                 #scene_aspectmode                       = 'auto', 
                 #scene                                  = dict(xaxis = dict(visible=False),
                                                            #yaxis = dict(visible=False),
                                                            #zaxis =dict(visible=False)), 
                 #scene_camera=dict(up    = dict(x=0, y=0, z=1),
                                   #center= dict(x=-0.05, y=0, z=-0.0),
                                   #eye   = dict(x=-1.0, y=-1.0, z=.4))   
        #)  
        
        #if save_figure:
            #fig_3d.write_image(hemisphere_filename + ".png")
            
        
        #if show_figure:
            #fig_3d.write_html( hemisphere_filename + '.html', auto_open=True)              
         



def colorax(vmin, vmax):
    return dict(cmin=vmin,  cmax=vmax)       
 
 

if __name__ == '__main__': 
    main()    
    plt.show()   
 
