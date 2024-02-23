import numpy as np
from RCAIDE.Core import Units , Data 
import matplotlib.pyplot as plt 
import matplotlib.colors
import matplotlib.colors as colors   
from RCAIDE.Visualization.Topography                                          import *  
import pickle


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():
     
    # topography data   
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
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
    
    save_figure     = True  
    file_type       = '.pdf'  
    use_lat_long    = True   
    city            = 'LA'  
    airport_tags    = ['LAX','LGB','BUR','DIS','SNA','ONT','SBD']
     
    # coordinates for departure airport 
    airport_coor    = np.array([[33.9348,-118.39678], [33.82347,-118.15542], [34.1843,-118.36587], 
                                           [33.81151,-117.9179886], [33.6919,-117.886], [34.04920,-117.59908], [34.0758,-117.2512]])
      

    # ---------------------------------------------
    # PLOT AIRPORTS WITH NOISE 
    # ---------------------------------------------
    results_filename = 'Raw_Data_HC/HC_1000ft_LA_10min_All' 
    Results          = load_results(results_filename)   
 
    elevation       = Results.elevation/Units.ft      
    colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
    colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200))  
    colors          = np.vstack((colors_undersea, colors_land))
    cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors) 
    norm = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation)) 
     
    
    figure_title  = city + '_Airport_Noise' 
 
    fig = plt.figure(figure_title)
    fig.set_size_inches(plot_parameters.figure_width,plot_parameters.figure_height)
    
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
    CS_2  = axis.contourf(LONG,LAT,Results.L_Aeq_jetliner ,noise_levels,cmap = noise_new_cmap)     
    cbar    = fig.colorbar(CS_2, ax=axis)        
    cbar.ax.set_ylabel('24-hr L$_{AeqT}$ [dBA]', rotation =  90)   
    
    # get aircraft noise and annotate  
    for airport in range(len(airport_coor)):   
        
        a_long    = airport_coor[airport,1]
        a_lat     = airport_coor[airport,0] 
        lat_deg   = Results.lat_deg   
        long_deg  = Results.long_deg  
        lat_dist  = Results.lat_dist  
        long_dist = Results.long_dist 
        
    
        if a_long>180: 
            a_long = a_long-360  
        elif a_long<0:
            a_long = a_long + 360
            
        a_dis_x   = ((a_lat - lat_deg[0,0])/(lat_deg[-1,0] - lat_deg[0,0]))*lat_dist[-1,0]/Units.nmi 
        a_dis_y   = ((a_long -long_deg[0,0]) /(long_deg[0,-1] - long_deg[0,0]))*long_dist[0,-1]/Units.nmi  
         
        if use_lat_long: 
            axis.annotate(airport_tags[airport], (a_long*0.9996,a_lat*0.9992), size=14)  
        else: 
            axis.annotate(airport_tags[airport], (a_dis_y*0.9996,a_dis_x*0.9992), size=14 )  
                                   
    # plot airports 
    airports_xy, airports_latlong = get_airport_runway_coordinates(long_dist,lat_dist,long_deg,lat_deg)
    if use_lat_long:
        for idx in range(len(airports_xy)): 
            axis.plot([airports_latlong[idx,0],airports_latlong[idx ,2]],[airports_latlong[idx ,1],airports_latlong[idx ,3]], 'k-', linewidth = 2.5)
    else:
        for idx in range(len(airports_xy)):
            axis.plot([airports_xy[idx,0],airports_xy[idx ,2]],[airports_xy[idx ,1],airports_xy[idx ,3]], 'k-', linewidth = 2.5)
         
        
    fig.tight_layout()  
    if save_figure:
        plt.savefig(figure_title + file_type)  
        
    return  

def get_airport_runway_coordinates(long_dist,lat_dist,long_deg,lat_deg):
    #x1,y1,x2,y2
    
    coords_dist = np.array([[8058 ,39600 , 12330,40216  ],
                            [9553 ,38339 , 14658,39057 ],
                            [36101,25148 , 34336,26853  ],
                            [34476,25511 , 36163,25504  ],
                            [30485,26430 , 35988,26450  ],
                            [60653,9194  , 61949,11473  ],
                            [16371,66265 , 15976,68042  ],
                            [16527,66678 , 15137,66681  ],
                            [3915.6 ,66700 , 3837.2 , 68685 ],
                            [119048,54907, 120743, 55323 ],
                            [ 88048 ,51480 ,84990 ,51460 ],
                            [82027 , 42771,84049 , 42772 ],
                            [82835,42524 , 83592, 43252 ]]) 
    
    # to do in the future (convert to coods_latlong)
    coords_latlong = np.zeros_like(coords_dist)
    for i in range(len(coords_dist)):
        x1  = coords_dist[i,0]/long_dist[0,-1]
        y1  = coords_dist[i,1]/lat_dist[-1,0]
        x2  = coords_dist[i,2]/long_dist[0,-1]
        y2  = coords_dist[i,3]/lat_dist[-1,0]

        x1_deg = x1*(long_deg[0,-1] - long_deg[0,0]) + long_deg[0,0]
        y1_deg = y1*(lat_deg[-1,0]  - lat_deg[0,0])  + lat_deg[0,0] 
        x2_deg = x2*(long_deg[0,-1] - long_deg[0,0]) + long_deg[0,0]
        y2_deg = y2*(lat_deg[-1,0]  - lat_deg[0,0])  + lat_deg[0,0]      
        
        coords_latlong[i,0] = x1_deg
        coords_latlong[i,1] = y1_deg
        coords_latlong[i,2] = x2_deg
        coords_latlong[i,3] = y2_deg
        
    return coords_dist/Units.nmi , coords_latlong



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
