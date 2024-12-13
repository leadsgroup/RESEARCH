import numpy as np
from RCAIDE.Framework.Core import Units, Data 
import matplotlib.pyplot as plt 
from scipy.interpolate import griddata  
import matplotlib.colors
import matplotlib.colors as colors  
from geopy.distance import geodesic as GD
import matplotlib.ticker as ticker
import pickle
from RCAIDE import  load
import  os
import  sys
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

local_path_1 =  os.path.split(sys.path[0])[0]

sys.path.append( os.path.join(local_path_1, 'Post_Processing_Functions'))
from Aircraft_Noise_Emissions   import generate_terrain_microphone_locations


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():  
    
    
    
    
    
    
    
    
    
    
    
    # Plot noise on top of topography
    
    # Universal Plot Settings
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"]    = "Times New Roman"
    parameters = {'axes.labelsize': 26,
                  'xtick.labelsize': 24,
                  'ytick.labelsize': 24,
                  'axes.titlesize': 26}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 3 
    plot_parameters.line_style       = '-' 
    plot_parameters.figure_width     = 11 
    plot_parameters.figure_height    = 6
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20 
    
    # topography data  
    use_lat_long                            = True
    save_figure                             = True
    file_type                               = '.png'
                 
    # tag for city              
    city                                    = 'LA' 
    aircraft_models                         = ['TRS']# ,'SR','TR']
    altitudes                               = ['1000']
    flight_intervals                        = ['60','30','10']
    microphone_x_resolution                 = 1200 
    microphone_y_resolution                 = 2700
    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator 
    topography_file = relative_path +  '..' + separator +  'City_Simulations' + separator + 'Los_Angeles' + separator + 'Topography' + separator + 'LA_Metropolitan_Area.txt'    
    microphone_locations , microphone_coordinates= generate_terrain_microphone_locations(topography_file, microphone_x_resolution, microphone_y_resolution)
    mic_loc =    np.reshape(microphone_locations, (microphone_x_resolution,microphone_y_resolution,3))
    mic_coord =     np.reshape(microphone_coordinates,   (microphone_x_resolution,microphone_y_resolution,3))
    for i in  range(len(aircraft_models)):
        # Read in results
        filename = relative_path +  '..' + separator +  'City_Simulations' + separator + 'Los_Angeles' + separator +'Tilt_Stopped_Rotor' + separator + 'Cumulative_' + aircraft_models[i] + '_'+ city +'_' +altitudes[0] +'ft.res'
        noise_data = load(filename)
        plot_2D_noise_contour(mic_coord, mic_loc, topography_file, 
                          noise_level              = noise_data.Total_L_dn, # d
                          min_noise_level          = 30,  
                          max_noise_level          = 90, 
                          noise_scale_label        = r'$L_{dn} [dBA]$',
                          save_figure              = False,
                          show_figure              = True,
                          save_filename            = "2D_Noise_Contour",
                          show_elevation           = False,
                          use_lat_long_coordinates = True,  
                          colormap                 = 'jet',
                          file_type                = ".png",
                          width = 10)  
    
    return 

def plot_2D_noise_contour(microphone_coordinates,
                          microphone_locations, 
                          topography_file, 
                          noise_level              = None ,
                          min_noise_level          = 35,  
                          max_noise_level          = 90, 
                          noise_scale_label        = None,
                          save_figure              = True,
                          show_figure              = True,
                          save_filename            = "2D_Noise_Contour",
                          show_elevation           = False,
                          use_lat_long_coordinates = True,  
                          colormap                 = 'jet',
                          file_type                = ".png",
                          width                    = 10, 
                          height                   = 7,
                          *args, **kwargs): 
    """This plots a 2D noise contour of a noise level 

    Assumptions:
    None

    Source:
    None

    Inputs: 
       noise_data        - noise data structure 
       noise_level       - noise level (dBA, DNL, SENEL etc)
       min_noise_level   - minimal noise level 
       max_noise_level   - maximum noise level 
       noise_scale_label - noise level label 
       save_figure       - save figure 
       show_figure       - show figure 
       save_filename     - save file flag
       show_trajectory   - plot aircraft trajectory flag
       show_microphones  - show microhpone flag 

    Outputs:
       Plots

    Properties Used:
    N/A
    """
        # settings

    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator  
    topography_file = relative_path +  '..' + separator +  'City_Simulations' + separator +  'Los_Angeles'    + separator +  'Topography' + separator + 'LA_Metropolitan_Area.txt'
    
    
    Textsize               = 20
    mic_x_res              = 1200
    mic_y_res              = 2700
    opacity                = 0.5
    save_figure            = True
    use_lat_long           = True    
    file_type              = '.png'   
    city                   = 'Los_Angeles'
    

    #microphone_locations ,microphone_coordinates  = generate_terrain_microphone_locations(topography_file, mic_x_res, mic_y_res)  
    
     
    lat_dist   = microphone_locations[:, :,0] #.reshape((mic_x_res, mic_y_res))
    long_dist  = microphone_locations[:,:, 1] #.reshape((mic_x_res, mic_y_res))
    lat_deg    = microphone_coordinates[:,:, 0] #.reshape((mic_x_res, mic_y_res))
    long_deg   = microphone_coordinates[:,:, 1] #.reshape((mic_x_res, mic_y_res))
    elevation  = microphone_coordinates[:,:, 2] #.reshape((mic_x_res, mic_y_res))
      
    elevation = elevation/Units.ft 
    colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56)) 
    colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200))  
    
    # combine them and build a new colormap
    colors          = np.vstack((colors_undersea, colors_land))
    cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors) 
    norm            = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation))   
    
    fig = plt.figure()
    fig.set_size_inches(9,5)
    axis = fig.add_subplot(1,1,1) 
    
    if use_lat_long: 
        CS  = axis.contourf(long_deg,lat_deg,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=opacity)  
        cbar = fig.colorbar(CS, ax=axis, orientation='horizontal',  shrink=0.6, pad=0.01, anchor=(0.76, 0))     
        cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation = 0,  fontsize=Textsize,  loc='center')
        cbar.ax.yaxis.set_label_coords(-0.32, 0)
        cbar.ax.tick_params(colors='black', labelsize=16)  
        axis.set_xlabel('Longitude [°]', fontsize=Textsize)
        axis.set_ylabel('Latitude [°]', fontsize=Textsize)   
    else:      
     
        CS   = axis.contourf(long_dist/Units.nmi,lat_dist/Units.nmi,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=opacity)     
        cbar = fig.colorbar(CS, ax=axis)        
        cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90) 
        axis.set_xlabel('Longitudinal Distance [nmi]')
        axis.set_ylabel('Latitudinal Distance [nmi]') 
    #axis.axis('equal')  


    
    elevation       = microphone_locations[:,2]/Units.ft      
    colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
    colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200))  
    colors          = np.vstack((colors_undersea, colors_land))
    cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors) 
    norm = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation))  
    parameters = {
                  'figure.dpi': 500
                  }
    plt.rcParams.update(parameters)
    #fig = plt.figure(save_filename)
    fig.set_size_inches(width,height)
    
    #axis = fig.add_subplot(1,1,1) 
    
    noise_levels   = np.linspace(35,100,256) #np.linspace(min_noise_level,max_noise_level,10)  
    noise_cmap     = plt.get_cmap('turbo')
    noise_new_cmap = truncate_colormap(noise_cmap,0.0, 1.0) 
     
    if use_lat_long_coordinates and (topography_file != None ):
        LAT  = microphone_coordinates[:,:,0]
        LONG = microphone_coordinates[:,:,1]
    else:
        LAT  = microphone_locations[:,:,0]
        LONG = microphone_locations[:,:,1]
        axis.set_xlabel('x [m]')
        axis.set_ylabel('y [m]')  
    
    if show_elevation:
        CS_1  = axis.contourf(LONG,LAT,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=0.5)  
        cbar = fig.colorbar(CS_1, ax=axis)     
        cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)  

    # plot aircraft noise levels   
    CS_2    = axis.contourf(LONG,LAT,noise_level ,noise_levels,cmap = noise_new_cmap)
    axis.set_ylim(LAT[0, 0], LAT[-1, 0])
    axis.set_xlim(LONG[0, 0], LONG[-1, -1])
    cbar    = fig.colorbar(CS_2, ax=axis,  orientation='vertical')        
    cbar.ax.set_ylabel(noise_scale_label, rotation =  90,  fontsize=Textsize)
    cbar.ax.tick_params(colors='black', labelsize=16)
    #plt.xlabel("Longitude", fontsize=16)
   # plt.ylabel("Latitude", fontsize=16)
    axis.tick_params(colors='black', labelsize=Textsize)
    # -------------------------
    # Create vertiport inset
    # -------------------------
    vertiport_coord = np.array([242.13719, 33.67619])
    
    index1 = np.argmin(abs(LAT[:, 0]-vertiport_coord[1]))
    index2 = np.argmin(abs(LONG[0, :]-vertiport_coord[0]))
    diff1 = 24
    diff2 = 54
    # Filter the x and y values for the zoomed range
    zoomed_lat = LAT[(index1-diff1):(index1+diff1), (index2-diff2):(index2+diff2)]#LONG[(LONG >= 242.35) & (LONG <= 242.45)]
    zoomed_long = LONG[(index1-diff1):(index1+diff1), (index2-diff2):(index2+diff2)] #LAT[np.where(LONG==zoomed_long[0]) : np.where(LONG==zoomed_long[-1])+1]
    zoomed_noise_level = noise_level[(index1-diff1):(index1+diff1), (index2-diff2):(index2+diff2)]
    
    y_size = abs(zoomed_lat[0, 0] - zoomed_lat[-1, -1])
    x_size = abs(zoomed_long[0, 0] - zoomed_long[-1, -1])
    rect = patches.Rectangle((zoomed_long[0, 0] , zoomed_lat[0, 0] ), x_size, y_size, linewidth=2, edgecolor='black', facecolor='red', alpha=0.5)
    axis.add_patch(rect)

    # Add inset
    inset_ax = fig.add_axes([0.14, 0.25, 0.2, 0.2])  # [x, y, width, height]
    
    # Access x and y data for the current case
    
    # Plot the zoomed data
    inset_ax.contourf(zoomed_long,  zoomed_lat,  zoomed_noise_level ,noise_levels, cmap = noise_new_cmap)    
    inset_ax.tick_params(colors='black', labelsize=7)
    
    #  Optionally, format tick labels as needed
    inset_ax.ticklabel_format(useOffset=False) 
    fig.tight_layout()

    
    #if save_figure: 
    figure_title  = save_filename
    plt.savefig(figure_title + file_type )
        
    
         
    return fig       

# ------------------------------------------------------------------ 
# Truncate colormaps
# ------------------------------------------------------------------  
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


class FixPointNormalize(matplotlib.colors.Normalize):
    """ 
    Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    Subclassing Normalize to obtain a colormap with a fixpoint 
    somewhere in the middle of the colormap.
    This may be useful for a `terrain` map, to set the "sea level" 
    to a color in the blue/turquise range. 
    """
    def __init__(self, vmin=None, vmax=None, sealevel=0, col_val = 0.21875, clip=False):
        # sealevel is the fix point of the colormap (in data units)
        self.sealevel = sealevel
        # col_val is the color value in the range [0,1] that should represent the sealevel.
        self.col_val = col_val
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
        return np.ma.masked_array(np.interp(value, x, y))   

def colorax(vmin, vmax):
    return dict(cmin=vmin,cmax=vmax)

def plot_topography(): 
    # settings

    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator  
    topography_file = relative_path +  '..' + separator +  'City_Simulations' + separator +  'Los_Angeles'    + separator +  'Topography' + separator + 'LA_Metropolitan_Area.txt'
    
    

    mic_x_res              = 1200
    mic_y_res              = 2700
    opacity                = 0.5
    save_figure            = True
    use_lat_long           = True    
    file_type              = '.png'   
    city                   = 'Los_Angeles'
    

    microphone_locations ,microphone_coordinates  = generate_terrain_microphone_locations(topography_file, mic_x_res, mic_y_res)  
    
     
    lat_dist   = microphone_locations[:,0].reshape((mic_x_res, mic_y_res))
    long_dist  = microphone_locations[:,1].reshape((mic_x_res, mic_y_res))
    lat_deg    = microphone_coordinates[:,0].reshape((mic_x_res, mic_y_res))
    long_deg   = microphone_coordinates[:,1].reshape((mic_x_res, mic_y_res))
    elevation  = microphone_coordinates[:,2].reshape((mic_x_res, mic_y_res))
      
    elevation = elevation/Units.ft 
    colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56)) 
    colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200))  
    
    # combine them and build a new colormap
    colors          = np.vstack((colors_undersea, colors_land))
    cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors) 
    norm            = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation))   
    
    fig = plt.figure()
    fig.set_size_inches(9,5)
    axis = fig.add_subplot(1,1,1) 
    
    if use_lat_long: 
        CS  = axis.contourf(long_deg,lat_deg,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=opacity)  
        cbar = fig.colorbar(CS, ax=axis)     
        cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)  
        axis.set_xlabel('Longitude [°]')
        axis.set_ylabel('Latitude [°]')   
    else:      
     
        CS   = axis.contourf(long_dist/Units.nmi,lat_dist/Units.nmi,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=opacity)     
        cbar = fig.colorbar(CS, ax=axis)        
        cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90) 
        axis.set_xlabel('Longitudinal Distance [nmi]')
        axis.set_ylabel('Latitudinal Distance [nmi]') 
    #axis.axis('equal') 
    fig.tight_layout()    

    return fig

 
if __name__ == '__main__': 
    main()   
    plt.show()
