# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
import MARC
from MARC.Core import Units
import pickle 
import time 
import numpy as np 
from MARC.Core import Data  
import matplotlib.cm as cm
import matplotlib.pyplot as plt  
import numpy as np 
from matplotlib.animation import FuncAnimation 
import matplotlib.animation as animation 
import matplotlib.colors
import matplotlib.colors as colors  

import matplotlib
matplotlib.use("Agg") 
from matplotlib.animation import FFMpegWriter

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(): 
    ti                         = time.time()

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
      
    # tag for city   
    city               = 'LA' 
    aircraft_models    = ['HC'] #['HC','SR','TR'] 
    altitudes          = ['1000'] #['1000','2000'] 
    flight_intervals   = ['60'] #['10','30','60']
    use_lat_long       = False 
    frames_per_sec     = 100 
    
    generate_video(frames_per_sec,city,aircraft_models,altitudes,use_lat_long,flight_intervals,plot_parameters)
    

    tf = time.time()
    print ('Time taken: ' + str(round(((tf-ti)/60),4)) + ' mins')      
    
    return 

def generate_video(frames_per_sec,city,aircraft_models,altitudes,use_lat_long,flight_intervals,PP):
     
   
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):
            for f_i in range(len(flight_intervals)):
                
                aircraft  = aircraft_models[m_i] 
                altitude  = altitudes[alt_i]
                frequency = flight_intervals[f_i]
                
                # load data  
                processed_results_filename = 'Raw_Data_' + aircraft + '/'+ aircraft + '_' + altitude + 'ft_' + city + '_' +  frequency +  'min_All'   
                Results         = load_results(processed_results_filename)  
                SPL_dBA         = Results.cumulative_noise_exposure
                timesteps       = len(SPL_dBA)     
                
                # ===========================================================

                elevation       = Results.elevation/Units.ft      
                colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
                colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200))  
                colors          = np.vstack((colors_undersea, colors_land))
                cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors) 
                norm = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation)) 
                
                altitude      = altitudes[alt_i]
                aircraft      = aircraft_models[m_i] 
                
                figure_title  = aircraft + '_'+ altitude + '_ft_mission_' + city 
    
              
                fig = plt.figure(figure_title)
                fig.set_size_inches(PP.figure_width,PP.figure_height)
                
                #axis           = fig.add_subplot(1,1,1)  
                noise_levels   = np.linspace(45,90,10) # np.array([45,50,55,60,70,80,90])   
                noise_cmap     = plt.get_cmap('turbo')
                noise_new_cmap = truncate_colormap(noise_cmap,0.0, 1.0) 
                 
                if use_lat_long: 
                    LAT  = Results.lat_deg
                    LONG = Results.long_deg
                    plt.xlabel('Longitude [°]')
                    plt.ylabel('Latitude [°]') 
                else:
                    LAT  = Results.lat_dist/Units.nmi
                    LONG = Results.long_dist/Units.nmi 
                    plt.xlabel('x [nmi]')
                    plt.ylabel('y [nmi]') 
                
                CS_1    = plt.contourf(LONG,LAT,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=0.5)   
                cbar_1  = plt.colorbar() 
                cbar_1.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)   
                
                # plot aircraft noise levels   
                CS_2     = plt.contourf(LONG,LAT,SPL_dBA[0],noise_levels,cmap = noise_new_cmap)     
                cbar_2     = plt.colorbar()      
                cbar_2.ax.set_ylabel('L$_{Amax}$ [dBA]', rotation =  90) 
                
                # ===========================================================
                
                
                writer = FFMpegWriter(fps=frames_per_sec) 
                with writer.saving(fig,  figure_title + '.mp4', 200):
                    for ts_idx in range(timesteps):  
                           
                        ## plot aircraft noise levels   
        
                        if use_lat_long: 
                            LAT  = Results.lat_deg
                            LONG = Results.long_deg
                            plt.xlabel('Longitude [°]')
                            plt.ylabel('Latitude [°]') 
                        else:
                            LAT  = Results.lat_dist/Units.nmi
                            LONG = Results.long_dist/Units.nmi 
                            plt.xlabel('x [nmi]')
                            plt.ylabel('y [nmi]') 
                            
                        for c_1 in CS_1.collections:
                            c_1.remove()  # removes only the contours, leaves the rest intact 
                        CS_1    = plt.contourf(LONG,LAT,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=0.5)  
                        #cbar_1  = plt.colorbar() 
                        #cbar_1.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)   
                        #plt.tight_layout()
                        

                        for c_2 in CS_2.collections:
                            c_2.remove()  # removes only the contours, leaves the rest intact                         
                        CS_2    = plt.contourf(LONG,LAT,SPL_dBA[ts_idx],noise_levels,cmap = noise_new_cmap)      
                        #cbar_2  = plt.colorbar()      
                        #cbar_2.ax.set_ylabel('L$_{Amax}$ [dBA]', rotation =  90) 
                        
                        plt.tight_layout()
                    
                        # ===========================================================
                         
                        writer.grab_frame() 
            
    return      


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