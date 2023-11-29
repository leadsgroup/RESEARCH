import numpy as np
from MARC.Core import Units, Data 
import matplotlib.pyplot as plt 
from scipy.interpolate import griddata  
import matplotlib.colors
import matplotlib.colors as colors  
from geopy.distance import geodesic as GD
import pickle


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():
    
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
    use_lat_long                = True
    save_figure                 = False
    file_type                   = '.pdf'
     
    # tag for city  
    city                        = 'LA' 
    aircraft_models             = ['HC','SR','TR'] 
    altitudes                   = ['1000','2000']
    flight_frequency            = '15' 
    
    plot_direct_flight_noise_LAmax(city,aircraft_models,altitudes,
                              flight_frequency, use_lat_long,save_figure,file_type,plot_parameters)  
    
    plot_direct_flight_noise_LAeqT(city,aircraft_models,altitudes,
                                  flight_frequency, use_lat_long,save_figure,file_type,plot_parameters)    
                          
    return 


 
# ---------------------------------------------
# PLOT DIRECT FLIGHT AIRCRAFT NOISE
# ---------------------------------------------    
def plot_direct_flight_noise_LAmax(city,aircraft_models,altitudes,
                              flight_frequency, use_lat_long,save_figure,file_type,PP): 
    
    
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):
            
            aircraft = aircraft_models[m_i] 
            altitude =  altitudes[alt_i]
            
            # load data  
            processed_results_filename = 'Raw_Data_' + aircraft + '/'+ aircraft + '_' + altitude + 'ft_' + city + '_' + flight_frequency +  'min_All'   
            Results = load_results(processed_results_filename) 

            elevation       = Results.elevation/Units.ft      
            colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
            colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200))  
            colors          = np.vstack((colors_undersea, colors_land))
            cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors) 
            norm = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation)) 
            
            altitude      = altitudes[alt_i]
            aircraft      = aircraft_models[m_i] 
            
            figure_title  = aircraft + '_'+ altitude + ' ft_mission_' + city + '_Lmax' 

            fig = plt.figure(figure_title)
            fig.set_size_inches(PP.figure_width,PP.figure_height)
            
            axis = fig.add_subplot(1,1,1) 
            
            noise_levels   = np.linspace(45,90,10) # np.array([45,50,55,60,70,80,90])   
            noise_cmap     = plt.get_cmap('turbo')
            noise_new_cmap = truncate_colormap(noise_cmap,0.0, 1.0) 
             
            if use_lat_long: 
                LAT  = Results.lat_deg
                LONG = Results.long_deg
                axis.set_xlabel('Longitude [°]')
                axis.set_ylabel('Latitude [°]') 
            else:
                LAT  = Results.lat_dist/Units.nmi
                LONG = Results.long_dist/Units.nmi 
                axis.set_xlabel('x [nmi]')
                axis.set_ylabel('y [nmi]') 
            
            CS_1  = axis.contourf(LONG,LAT,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=0.5)  
            cbar = fig.colorbar(CS_1, ax=axis)     
            cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)   
            
            # plot aircraft noise levels   
            CS_2     = axis.contourf(LONG,LAT,Results.L_Amax ,noise_levels,cmap = noise_new_cmap)     
            cbar     = fig.colorbar(CS_2, ax=axis)        
            cbar.ax.set_ylabel('L$_{Amax}$ [dBA]', rotation =  90) 
             
            fig.tight_layout()  
            if save_figure:
                plt.savefig(figure_title + file_type)    
 
    return  
 

def plot_direct_flight_noise_LAeqT(city,aircraft_models,altitudes,
                              flight_frequency, use_lat_long,save_figure,file_type,PP):  

    
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):
            aircraft = aircraft_models[m_i] 
            altitude =  altitudes[alt_i]
            
            # load data  
            processed_results_filename = 'Raw_Data_' + aircraft + '/'+ aircraft + '_' + altitude + 'ft_' + city + '_' + flight_frequency +  'min_All'   
            Results = load_results(processed_results_filename)

            elevation       = Results.elevation/Units.ft      
            colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
            colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200))  
            colors          = np.vstack((colors_undersea, colors_land))
            cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors) 
            norm = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation)) 
            
            altitude      = altitudes[alt_i]
            aircraft      = aircraft_models[m_i] 
            
            figure_title  = aircraft + '_'+ altitude + ' ft_mission_' + city + '_L_AeqT' 

            fig = plt.figure(figure_title)
            fig.set_size_inches(PP.figure_width,PP.figure_height)
            
            axis = fig.add_subplot(1,1,1) 
            
            noise_levels   = np.linspace(45,90,10) # np.array([45,50,55,60,70,80,90])   
            noise_cmap     = plt.get_cmap('turbo')
            noise_new_cmap = truncate_colormap(noise_cmap,0.0, 1.0) 
             
            if use_lat_long: 
                LAT  = Results.lat_deg
                LONG = Results.long_deg
                axis.set_xlabel('Longitude [°]')
                axis.set_ylabel('Latitude [°]') 
            else:
                LAT  = Results.lat_dist/Units.nmi
                LONG = Results.long_dist/Units.nmi 
                axis.set_xlabel('x [nmi]')
                axis.set_ylabel('y [nmi]') 
            
            CS_1  = axis.contourf(LONG,LAT,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=0.5)  
            cbar = fig.colorbar(CS_1, ax=axis)     
            cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)   
            

            # plot aircraft noise levels   
            CS_2  = axis.contourf(LONG,LAT,Results.L_Aeq_total ,noise_levels,cmap = noise_new_cmap)     
            cbar    = fig.colorbar(CS_2, ax=axis)        
            cbar.ax.set_ylabel('24-hr L$_{AeqT}$ [dBA]', rotation =  90)   
             
            fig.tight_layout()  
            if save_figure:
                plt.savefig(figure_title + file_type)    
 
  
    return  

# ------------------------------------------------------------------
# SUPPORTING FUNCTIONS 
# ------------------------------------------------------------------ 

def get_topography_data(topography_file,N_lat,N_long ): 
 
    data = np.loadtxt(topography_file)
    Long = data[:,0]
    Lat  = data[:,1]
    Elev = data[:,2]  
    
    x_min_coord = np.min(Lat)
    x_max_coord = np.max(Lat)
    y_min_coord = np.min(Long)
    y_max_coord = np.max(Long)
    if np.min(Long)>180: 
        y_min_coord = np.min(Long)-360
    if np.max(Long)>180:
        y_max_coord = np.max(Long)-360  
    
    top_left_map_coords      = np.array([x_max_coord,y_min_coord])
    bottom_left_map_coords   = np.array([x_min_coord,y_min_coord])  
    bottom_right_map_coords  = np.array([x_min_coord,y_max_coord]) 
    
    x_dist_max = GD(top_left_map_coords,bottom_left_map_coords).m 
    y_dist_max = GD(bottom_right_map_coords,bottom_left_map_coords).m  
     
    [long_dist,lat_dist]  = np.meshgrid(np.linspace(0,y_dist_max,N_long),np.linspace(0,x_dist_max,N_lat))
    [long_deg,lat_deg]    = np.meshgrid(np.linspace(np.min(Long),np.max(Long),N_long),np.linspace(np.min(Lat),np.max(Lat),N_lat)) 
    elevation             = griddata((Lat,Long), Elev, (lat_deg, long_deg), method='linear')     
     
    return long_dist,lat_dist,long_deg,lat_deg,elevation

 
# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def save_results(results,filename): 
    pickle_file  =  filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return   


# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results  


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
    
if __name__ == '__main__': 
    main()   
    plt.show()
