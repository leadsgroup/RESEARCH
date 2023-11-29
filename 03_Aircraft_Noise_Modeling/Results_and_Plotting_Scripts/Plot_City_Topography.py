import numpy as np
from MARC.Core import Units
import matplotlib.pyplot as plt
from matplotlib.image import imread
from scipy.interpolate import griddata  
from MARC.Methods.Missions.compute_point_to_point_geospacial_data           import compute_point_to_point_geospacial_data
import matplotlib.colors
import matplotlib.colors as colors  
from geopy.distance import geodesic as GD
from MARC.Visualization.Topography                                          import * 
import scipy.ndimage as ndimage
import pickle


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main(): 
    # settings 
    topography_file        = '../Maps_and_Scales/Los_Angeles/LA_Metropolitan_Area_5.txt'   
    N_lat                  = 225  
    N_long                 = 390    
    save_figure            = True
    use_lat_long           = True    
    file_type              = '.pdf'   
    city                   = 'Los_Angeles'  
          
    long_dist,lat_dist,long_deg,lat_deg,elevation  = get_topography_data(topography_file,N_lat,N_long)  
      
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
        CS  = axis.contourf(long_deg,lat_deg,elevation,cmap =cut_terrain_map,norm=norm,levels = 20)  
        cbar = fig.colorbar(CS, ax=axis)     
        cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)  
        axis.set_xlabel('Longitude [°]')
        axis.set_ylabel('Latitude [°]')   
    else:      
     
        CS   = axis.contourf(long_dist/Units.nmi,lat_dist/Units.nmi,elevation,cmap =cut_terrain_map,norm=norm,levels = 20)     
        cbar = fig.colorbar(CS, ax=axis)        
        cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90) 
        axis.set_xlabel('Longitudinal Distance [nmi]')
        axis.set_ylabel('Latitudinal Distance [nmi]') 
    #axis.axis('equal') 
    fig.tight_layout()    

    if save_figure:
    
        save_filename = '../Papers_and_Presentation/Images/'  + city + '_Topography'        
        plt.savefig(save_filename + file_type )
    return  

    
def get_topography_data(topography_file,N_lat,N_long ): 

    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams.update({'axes.labelsize': 18,
                  'xtick.labelsize': 18,
                  'ytick.labelsize': 18,
                  'axes.titlesize': 32,
                  'font.size': 10}) 
    
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
    elevation                 = griddata((Lat,Long), Elev, (lat_deg, long_deg), method='linear')     
     
    return long_dist,lat_dist,long_deg,lat_deg,elevation 
 
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
