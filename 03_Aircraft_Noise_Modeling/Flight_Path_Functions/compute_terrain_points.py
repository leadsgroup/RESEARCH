# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------   
from RCAIDE.Framework.Core    import Units,  Data  

from RCAIDE.Framework.Analyses.Geodesics.Geodesics import Calculate_Distance
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors 
import numpy as np 

# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Library-Plots-Topograpgy 
def compute_terrain_points(topography_file,
                            number_of_latitudinal_points  = 100,
                            number_of_longitudinal_points = 100,   ): 
     
  
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
    top_right_map_coords     = np.array([x_max_coord,y_max_coord])
    bottom_right_map_coords  = np.array([x_min_coord,y_max_coord]) 
    
    x_dist_max = Calculate_Distance(top_left_map_coords,bottom_left_map_coords) * Units.kilometers
    y_dist_max = Calculate_Distance(bottom_right_map_coords,bottom_left_map_coords) * Units.kilometers
    
    [long_dist,lat_dist]  = np.meshgrid(np.linspace(0,y_dist_max,number_of_longitudinal_points),np.linspace(0,x_dist_max,number_of_latitudinal_points))
    [long_deg,lat_deg]    = np.meshgrid(np.linspace(np.min(Long),np.max(Long),number_of_longitudinal_points),np.linspace(np.min(Lat),np.max(Lat),number_of_latitudinal_points)) 
    elevation             = griddata((Lat,Long), Elev, (lat_deg, long_deg), method='linear')     
    elevation             = elevation/Units.feet
     
    terrain_data = Data( 
        top_left_map_coordinates      = top_left_map_coords,     
        bottom_left_map_coordinates   = bottom_left_map_coords,  
        top_right_map_coordinates     = top_right_map_coords,    
        bottom_right_map_coordinates  = bottom_right_map_coords,   
        longitudinal_distance         = long_dist,
        latitudinal_distance          = lat_dist,
        longitudinal_coordinates      = long_deg,
        latitudinal_coordinates       = lat_deg,
        elevation                     = elevation)
    
    return terrain_data

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
    
  