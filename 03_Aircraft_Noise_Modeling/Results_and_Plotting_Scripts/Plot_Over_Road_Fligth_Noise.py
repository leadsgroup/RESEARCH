import numpy as np
from RCAIDE.Core import Units
import matplotlib.pyplot as plt
from matplotlib.image import imread
from scipy.interpolate import griddata 
import matplotlib.colors
import matplotlib.colors as colors  
from geopy.distance import geodesic as GD
import pickle


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():
     
    # topography data 
    topography_file        = '../Maps_and_Scales/Los_Angeles/LA_Metropolitan_Area_5.txt'   
    N_lat                  = 225  
    N_long                 = 390    
    use_lat_long           = False
    save_figure            = True 
    
    include_aircraft_noise = True 
    include_car_noise      = True 
 
    
    long_dist,lat_dist,long_deg,lat_deg,elevation  = get_topography_data(topography_file,N_lat,N_long)
    
    # community noise data 
    total_noise_file    = '../Maps_and_Scales/Los_Angeles/LA_Aircraft_Noise.png'
    car_noise_file      = '../Maps_and_Scales/Los_Angeles/LA_Aircraft_Noise.png'
    train_noise_file    = '../Maps_and_Scales/Los_Angeles/LA_Aircraft_Noise.png'
    aircraft_noise_file = '../Maps_and_Scales/Los_Angeles/LA_Aircraft_Noise_Modified.png'
    aircraft_noise,car_noise, train_noise = vehicle_noise_map(total_noise_file,car_noise_file,train_noise_file,
                                                              aircraft_noise_file,long_dist,lat_dist,long_deg,lat_deg,N_lat,N_long) 
     
    # tag for city 
    city                        = 'Los_Angeles' 
    city_acronym                = 'LA' 
    aircraft_models             = ['HC']
    
    # altitude 
    altitudes                   = ['1000']
    
    # departure airport location  
    departure_location         = ['LAX','LGB','BUR','LAX','BUR','LAX','DIS','SNA','SNA']
    
    # coordinates for departure airport 
    departure_coord            = np.array([[33.9348,-118.39678],
                                           [33.81347,-118.15542],
                                           [34.1843,-118.36587],
                                           [33.9348,-118.39678],
                                           [34.1843,-118.36587],
                                           [33.9348,-118.39678],
                                           [33.81151,-117.9179886],
                                           [33.6719,-117.886],
                                           [33.6719,-117.886]])
    
    # destination airport location
    destination_location       = ['SBD','ONT','ONT','BUR','DIS','DIS','SBD','BUR','SBD']
    
    # coordinates for destination airport 
    destination_coord          = np.array([[34.0758,-117.2512],  
                                           [34.04920,-117.59908],
                                           [34.04920,-117.59908],
                                           [34.1843,-118.36587],
                                           [33.81151,-117.9179886],
                                           [33.81151,-117.9179886],
                                           [34.0758,-117.2512],
                                           [34.1843,-118.36587],
                                           [34.0758,-117.2512]])
               
    # ---------------------------------------------      
    # PLOT AIRCRAFT NOISE OVER ROADWAY
    # ---------------------------------------------  
    
    LAX_to_DIS_aircraft_trajectory()
    
    # departure airport location  
    departure_location         = ['LAX']
    
    # coordinates for departure airport 
    departure_coord            = np.array([[33.9348,-118.39678]])
    
    # destination airport location
    destination_location       = ['DIS']
    
    # coordinates for destination airport 
    destination_coord          = np.array([[33.81151,-117.9179886]])  
    
    plot_roadway_noise(topography_file,city,city_acronym,aircraft_models,altitudes,departure_location,departure_coord ,
                                 destination_location,destination_coord,long_dist,lat_dist,long_deg,lat_deg,elevation,
                                 include_aircraft_noise,include_car_noise,aircraft_noise,car_noise,train_noise,use_lat_long,save_figure)
    
    
    plot_max_aircraft_noise_over_roadway(topography_file,city,city_acronym,aircraft_models,altitudes,departure_location,departure_coord ,
                              destination_location,destination_coord,long_dist,lat_dist,long_deg,lat_deg,elevation,
                              include_aircraft_noise,include_car_noise,aircraft_noise,car_noise,train_noise,use_lat_long,save_figure)     
    return 



 
 
  

# ---------------------------------------------      
# PLOT AIRCRAFT NOISE OVER ROADWAY
# ---------------------------------------------  
def plot_roadway_noise(topography_file,city,city_acronym,aircraft_models,altitudes,departure_location,departure_coord ,
                              destination_location,destination_coord,long_dist,lat_dist,long_deg,lat_deg,elevation,
                              include_aircraft_noise,include_car_noise,aircraft_noise,car_noise,train_noise,use_lat_long,save_figure):  
     

    elevation = elevation/Units.ft     
    colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
    colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200)) 
    
    # combine them and build a new colormap
    colors          = np.vstack((colors_undersea, colors_land))
    cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors)
    
    norm = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation)) 
    file_type = '.pdf'
    
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):
            altitude      = altitudes[alt_i]
            aircraft      = aircraft_models[m_i] 
            figure_title  = aircraft + '_'+ altitude + ' ft_mission_' + city
            save_filename = aircraft + '_'+ altitude + ' ft_mission_' + city 

            fig = plt.figure(figure_title)
            fig.set_size_inches(11,5)
            axis = fig.add_subplot(1,1,1) 
            
            noise_levels   = np.linspace(45,90,10) # np.array([45,50,55,60,70,80,90])   
            noise_cmap     = plt.get_cmap('turbo')
            noise_new_cmap = truncate_colormap(noise_cmap,0.0, 1.0) 
             
            if use_lat_long: 
                LAT  = lat_deg
                LONG = long_deg
                axis.set_xlabel('Longitude [°]')
                axis.set_ylabel('Latitude [°]') 
            else:
                LAT  = lat_dist/Units.nmi
                LONG = long_dist/Units.nmi 
                axis.set_xlabel('Longitude Distance [nmi]')
                axis.set_ylabel('Latitude Distance [nmi]') 
            
            CS_1  = axis.contourf(LONG,LAT,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=0.5)  
            cbar = fig.colorbar(CS_1, ax=axis)     
            cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)   
                 
            # plot car noise level  
            long_dist_high_res,lat_dist_high_res,long_de_high_resg,lat_deg_high_res,road_noise = LAX_to_DIS_car_noise()
            if use_lat_long:
                CS_2     = axis.contourf(long_de_high_resg,lat_deg_high_res,road_noise ,noise_levels,cmap = noise_new_cmap)  
            else:
                CS_2     = axis.contourf(long_dist_high_res/Units.nmi,lat_dist_high_res/Units.nmi,road_noise,noise_levels,cmap = noise_new_cmap) 
                     
            cbar     = fig.colorbar(CS_2, ax=axis)        
            cbar.ax.set_ylabel('LA$_{max}$ [dBA]', rotation =  90)  
             
            fig.tight_layout()  
            if save_figure:
                plt.savefig(save_filename + file_type)    
 
    return  


def count_between(X, lower_bound, upper_bound):
    count = 0
    x = X.flatten()
    for i in range(len(x)):
        if x[i] >= lower_bound and x[i] < upper_bound: 
            count += 1  
    return count  


# ---------------------------------------------      
# PLOT AIRCRAFT NOISE OVER ROADWAY
# ---------------------------------------------  
def plot_max_aircraft_noise_over_roadway(topography_file,city,city_acronym,aircraft_models,altitudes,departure_location,departure_coord ,
                              destination_location,destination_coord,long_dist,lat_dist,long_deg,lat_deg,elevation,
                              include_aircraft_noise,include_car_noise,aircraft_noise,car_noise,train_noise,use_lat_long,save_figure):  
     

    elevation = elevation/Units.ft     
    colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
    colors_land     = plt.cm.terrain(np.linspace(0.25, 1, 200)) 
    
    # combine them and build a new colormap
    colors          = np.vstack((colors_undersea, colors_land))
    cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors)
    
    norm = FixPointNormalize(sealevel=0,vmax=np.max(elevation),vmin=np.min(elevation)) 
    file_type = '.pdf'
    
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):
            altitude      = altitudes[alt_i]
            aircraft      = aircraft_models[m_i] 
            figure_title  = aircraft + '_'+ altitude + ' ft_mission_' + city + 'and_Car'
            save_filename = aircraft + '_'+ altitude + ' ft_mission_' + city + 'and_Car'

            fig = plt.figure(figure_title)
            fig.set_size_inches(11,5)
            axis = fig.add_subplot(1,1,1) 
            
            noise_levels   = np.linspace(45,90,10) # np.array([45,50,55,60,70,80,90])   
            noise_cmap     = plt.get_cmap('turbo')
            noise_new_cmap = truncate_colormap(noise_cmap,0.0, 1.0) 
             
            if use_lat_long: 
                LAT  = lat_deg
                LONG = long_deg
                axis.set_xlabel('Longitude [°]')
                axis.set_ylabel('Latitude [°]') 
            else:
                LAT  = lat_dist/Units.nmi
                LONG = long_dist/Units.nmi  
                axis.set_xlabel('Longitude Distance [nmi]')
                axis.set_ylabel('Latitude Distance [nmi]') 
                
            
            CS_1  = axis.contourf(LONG,LAT,elevation,cmap =cut_terrain_map,norm=norm,levels = 20, alpha=0.5)  
            cbar = fig.colorbar(CS_1, ax=axis)     
            cbar.ax.set_ylabel('Elevation above sea level [ft]', rotation =  90)  
             
            
            # get aircraft noise and annotate  
            flight_no = 0
            filename      =  aircraft + '_1000ft_mission_LA_LAX_to_DIS_Over_Road'
            results = load_results(filename)         
        
            noise_data      = post_process_noise_data(results)  
            SPL_contour_gm  = noise_data.SPL_dBA_ground_mic
            
            # plot car noise level 
            if include_car_noise:
                long_dist_high_res,lat_dist_high_res,long_de_high_resg,lat_deg_high_res,road_noise = LAX_to_DIS_car_noise()
                if use_lat_long:
                    CS_2     = axis.contourf(long_de_high_resg,lat_deg_high_res,road_noise ,noise_levels,cmap = noise_new_cmap)  
                else:
                    CS_2     = axis.contourf(long_dist_high_res/Units.nmi,lat_dist_high_res/Units.nmi,road_noise,noise_levels,cmap = noise_new_cmap) 
             
            CS_3     = axis.contourf(LONG,LAT,np.max(SPL_contour_gm,axis = 0) ,noise_levels,cmap = noise_new_cmap)     
            cbar     = fig.colorbar(CS_3, ax=axis)        
            cbar.ax.set_ylabel('LA$_{max}$ [dBA]', rotation =  90) 
             
            fig.tight_layout()  
            if save_figure:
                plt.savefig(save_filename + file_type)    
 
    return   
 


# ------------------------------------------------------------------
# SUPPORTING FUNCTIONS 
# ------------------------------------------------------------------
def vehicle_noise_map(total_noise_file,car_noise_file,train_noise_file, aircraft_noise_file,long_dist,lat_dist,long_deg,lat_deg,N_lat,N_long):
    
    # import image 
    #total_noise_raw_data          = imread(total_noise_file ) 
    car_noise_raw_data            = imread(car_noise_file )
    train_noise_raw_data          = imread(train_noise_file )
    aircraft_noise_raw_data       = imread(aircraft_noise_file )
    
    # import test samples of colors 
    blue            = imread('../Maps_and_Scales/Color_Scales/Color_Blue.png' )   # > 90
    magenta         = imread('../Maps_and_Scales/Color_Scales/Color_Magenta.png' ) # 80 - 89.0
    purple          = imread('../Maps_and_Scales/Color_Scales/Color_Purple.png' ) #  70 - 79.9
    pink            = imread('../Maps_and_Scales/Color_Scales/Color_Pink.png' )   #  60 - 64.9 
    red             = imread('../Maps_and_Scales/Color_Scales/Color_Red.png' )    #  55 - 59.9
    orange          = imread('../Maps_and_Scales/Color_Scales/Color_Orange.png') #   50 - 54.9
    yellow_scale    = imread('../Maps_and_Scales/Color_Scales/Color_Yellow.png') #   45 - 49.9
    grey_background = imread('../Maps_and_Scales/Color_Scales/Color_Gray.png' )   # background
    black_sea       = imread('../Maps_and_Scales/Color_Scales/Color_Black.png' )  # sea  
    
    # associate color with scale 
    color_scales        = np.zeros((9,4))
    color_scales[0,0:3] = np.mean(np.mean(blue[:,:,0:3], axis = 0 ), axis = 0 )      
    color_scales[1,0:3] = np.mean(np.mean(magenta[:,:,0:3], axis = 0 ), axis = 0 )         
    color_scales[2,0:3] = np.mean(np.mean(purple[:,:,0:3], axis = 0 ), axis = 0 )          
    color_scales[3,0:3] = np.mean(np.mean(pink[:,:,0:3], axis = 0 ), axis = 0 )            
    color_scales[4,0:3] = np.mean(np.mean(red[:,:,0:3], axis = 0 ), axis = 0 )            
    color_scales[5,0:3] = np.mean(np.mean(orange[:,:,0:3], axis = 0 ), axis = 0 )         
    color_scales[6,0:3] = np.mean(np.mean(yellow_scale[:,:,0:3], axis = 0 ), axis = 0 )    
    color_scales[7,0:3] = np.mean(np.mean(grey_background[:,:,0:3], axis = 0 ), axis = 0 )  
    color_scales[8,0:3] = np.mean(np.mean(black_sea[:,:,0:3], axis = 0 ), axis = 0 )         
     
    color_scales[0,3] =  95 
    color_scales[1,3] = (80 + 89.9)/2 # medium value
    color_scales[2,3] = (70 + 79.9)/2 # medium value
    color_scales[3,3] = (60 + 69.9)/2 # medium value
    color_scales[4,3] = (55 + 59.9)/2 # medium value
    color_scales[5,3] = (50 + 54.9)/2 # medium value
    color_scales[6,3] = (45 + 49.9)/2 # medium value
    color_scales[7,3] = 32 # background noise 
    color_scales[8,3] = 0  # water   
    
     
    aircraft_noise = vectorize_image(aircraft_noise_raw_data ,color_scales,N_lat,N_long) 
    car_noise      = vectorize_image(car_noise_raw_data ,color_scales,N_lat,N_long) 
    train_noise    = vectorize_image(train_noise_raw_data ,color_scales,N_lat,N_long)  
    
    return aircraft_noise, car_noise, train_noise

def vectorize_image(raw_data,color_scales,N_lat,N_long): 
    
    num_pixels_in_col = raw_data.shape[0]
    num_pixels_in_row = raw_data.shape[1]

    tol_1 = 0.1   # 0.095
    tol_2 = 0.05  # 0.05
    tol_3 = 0.25  # 0.075
    vectorized_data = np.zeros((num_pixels_in_col,num_pixels_in_row)) 

    data = []
    x_data    = []
    y_data    = []
    for s in range(len(color_scales)):  
        for col_i in range(num_pixels_in_col):
            for row_i in range(num_pixels_in_row):
                pixel_data = raw_data[col_i,row_i]
                if  abs(color_scales[s,0] - pixel_data[0]) < tol_1:
                    if abs(color_scales[s,1] -  pixel_data[1]) < tol_2:
                        if abs(color_scales[s,2] - pixel_data[2]) < tol_3:
                            vectorized_data[col_i,row_i] = color_scales[s,3]
                            data.append(color_scales[s,3])
                            y_data.append(col_i)
                            x_data.append(row_i)

    x     = np.array(x_data)
    y     = np.array(y_data) 
    z     = np.array(data)    

    xi    = np.linspace(0,num_pixels_in_row,N_long) # longitude direction 
    yi    = np.linspace(0,num_pixels_in_col,N_lat)  # latitude direction  
    
    # grid the data.
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic') 
    
    return np.flip(zi,axis = 0)
    

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


def LAX_to_DIS_car_noise():

    # topography data 
    topography_file        = '../Maps_and_Scales/Los_Angeles/LA_Metropolitan_Area_5.txt'   
    N_lat                  = 225*6   
    N_long                 = 390*6    
    
    long_dist,lat_dist,long_deg,lat_deg,elevation  = get_topography_data(topography_file,N_lat,N_long)
    
    road_coordinates = np.array([[12020.26460318356  , 38706.4610664855],
                                 [12191.027431763929 , 37013.4541263636],
                                 [14844.019703490529 , 37016.0139323786],
                                 [15981.638721251782 , 37379.8720730845],
                                 [17622.3077240734   , 36414.09383227629],
                                 [21034.33500418593  , 37021.98681308036],
                                 [27098.109896268204 , 36906.9174379291],
                                 [31263.98610244645  , 35097.13458531274],
                                 [39980.753266255335 , 34984.6250161765],
                                 [54113.03503212024  , 25082.80776799602]])
    
    # Refine Points
    interval = 10
    X = np.array([road_coordinates[0][0]]) 
    Y = np.array([road_coordinates[0][1]])
    for i in range(1, len(road_coordinates)):
        # find distance
        distance = ((road_coordinates[i,0] - road_coordinates[i-1,0])**2 + (road_coordinates[i,1] - road_coordinates[i-1,1])**2)**0.5
        num_pts  = int(np.floor(distance/interval))
        
        # interpolate
        x0 = road_coordinates[i-1:i+1,0]
        y0 = road_coordinates[i-1:i+1,1] 
        xinterp  = np.linspace(x0[0],x0[1],num_pts)
        yinterp = np.interp(xinterp, x0 , y0)
        
        # append 
        X = np.append(X,xinterp[1:])  
        Y = np.append(Y,yinterp[1:])  
    
    # Loop through gird and store noise 
    stencil_size = 5
    road_noise   = np.ones_like(lat_dist)*32
    for j in range(len(X)):
        # find index of nearest value in x  
        x_center_idx  = abs(long_dist[0,:] - X[j]).argmin() 
        y_center_idx  = abs(lat_dist[:,0] - Y[j]).argmin()     

        # find index of nearest value in y  
        for k in range((x_center_idx-stencil_size),(x_center_idx+stencil_size)):
            for l in range((y_center_idx-stencil_size),(y_center_idx+stencil_size)): 
                del_x           = abs(long_dist[l,k] - X[j]) 
                del_y           = abs(lat_dist[l,k] - Y[j])
                dist            = np.sqrt(del_x**2 + del_y**2)  
                road_noise[l,k] = np.maximum(road_noise[l,k],np.maximum(1E-04*dist**2 - 0.15*dist+ 67.656,32) ) 
    
    return long_dist,lat_dist,long_deg,lat_deg,road_noise
 

def LAX_to_DIS_aircraft_trajectory():
 
    
    road_coordinates = np.array([[12020.26460318356, 38706.4610664855],
                                 [12191.027431763929, 37013.4541263636], 
                                 [27098.109896268204, 36906.9174379291],
                                 [31263.98610244645, 35097.13458531274],
                                 [39980.753266255335, 34984.6250161765],
                                 [54113.03503212024, 25082.80776799602]]) 
    
    # Refine Points 
    distances = [0]
    angles    = [0]
    for i in range(1, len(road_coordinates)):
        # find distance
        delta_x = (road_coordinates[i,0] - road_coordinates[i-1,0])/Units.nmi
        delta_y = (road_coordinates[i,1] - road_coordinates[i-1,1])/Units.nmi
        distances.append( distances[i-1] + (delta_x**2 + delta_y**2)**0.5)
        
        if delta_x>0 and delta_y>0: 
            angle = np.arctan(abs(delta_x/delta_y))
        elif delta_x>0 and delta_y<0: 
            angle = np.pi - np.arctan(abs(delta_x/delta_y)) 
        elif delta_x<0 and delta_y<0: 
            angle = np.pi + np.arctan(abs(delta_x/delta_y))  
        elif delta_x<0 and delta_y>0: 
            angle = 2*np.pi - np.arctan(abs(delta_x/delta_y))         
        angles.append(angle/Units.degrees)
    
    print('Distances: ', distances)
    print('True Course Angles: ', angles)
    
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
    coords_latlong = None 
    return coords_dist/Units.nmi , coords_latlong


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
