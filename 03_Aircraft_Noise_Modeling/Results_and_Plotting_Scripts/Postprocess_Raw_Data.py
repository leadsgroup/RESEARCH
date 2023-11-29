import numpy as np
from MARC.Core import Units , Data 
import matplotlib.pyplot as plt
from matplotlib.image import imread
from scipy.interpolate import griddata 
from MARC.Visualization.Performance.Common.post_process_noise_data import post_process_noise_data 
from MARC.Methods.Noise.Fidelity_Zero.Noise_Tools.decibel_arithmetic         import SPL_arithmetic 
import matplotlib.colors
import matplotlib.colors as colors  
from geopy.distance import geodesic as GD
import time  
import pickle


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():
    #test              = np.array([[10,10]])
    #res                = SPL_arithmetic(test, sum_axis=1)
    
    #A_SPL             = np.ones((2,4))*10
    #B_SPL             = np.ones((2,4))*10
    #B_SPL[0,1]      = np.nan
    #B_SPL[1,1]      = np.nan
    #B_SPL[0,3]      = np.nan
    
    
    #p_prefs_A         =  10**(A_SPL[:,:,None]/10)
    #p_prefs_B         =  10**(B_SPL[:,:,None]/10)
    #C                 =  np.concatenate((p_prefs_A,p_prefs_B),axis = 2)
    #D                 =  10*np.log10(np.nansum(C,axis=2))  
        

    ti                         = time.time()
    
    # topography data 
    topography_file        = '../Maps_and_Scales/Los_Angeles/LA_Metropolitan_Area_5.txt'   
    N_lat                  = 225  
    N_long                 = 390     
    
    long_dist,lat_dist,long_deg,lat_deg,elevation  = get_topography_data(topography_file,N_lat,N_long)
    
    # community noise data 
    total_noise_file    = '../Maps_and_Scales/Los_Angeles/LA_Aircraft_Noise.png'
    car_noise_file      = '../Maps_and_Scales/Los_Angeles/LA_Aircraft_Noise.png'
    train_noise_file    = '../Maps_and_Scales/Los_Angeles/LA_Aircraft_Noise.png'
    aircraft_noise_file = '../Maps_and_Scales/Los_Angeles/LA_Aircraft_Noise_Modified_2.png'
    jetliner_L_AeqT_24hr,_,_ = community_noise_map(total_noise_file,car_noise_file,train_noise_file,aircraft_noise_file,long_dist,lat_dist,long_deg,lat_deg,N_lat,N_long)  
    
    # tag for city 
    city                        = 'Los_Angeles' 
    city_acronym                = 'LA' 
    aircraft_models             =  ['SR']
    
    # altitude 
    altitudes                   = ['1000']
    
    # departure airport location  
    dep_loc         = ['LAX','LGB','BUR','LAX','BUR','LAX','DIS','SNA','SNA','SBD','ONT','ONT','BUR','DIS','DIS','SBD','BUR','SBD'] 
    
    # destination airport location
    des_loc       =   ['SBD','ONT','ONT','BUR','DIS','DIS','SBD','BUR','SBD','LAX','LGB','BUR','LAX','BUR','LAX','DIS','SNA','SNA'] 

    # 1hr min interval 
    #flight_times = np.array(['06:00:00',
                             #'07:00:00',
                             #'08:00:00',
                             #'09:00:00',
                             #'10:00:00',
                             #'11:00:00',
                             #'12:00:00',
                             #'13:00:00',
                             #'14:00:00',
                             #'15:00:00',
                             #'16:00:00',
                             #'17:00:00',
                             #'18:00:00',
                             #'19:00:00',]) 
                             
    #30 min interval 
    #flight_times = np.array(['06:00:00','06:30:00',
                             #'07:00:00','07:30:00',
                             #'08:00:00','08:30:00',
                             #'09:00:00','09:30:00',
                             #'10:00:00','10:30:00',
                             #'11:00:00','11:30:00',
                             #'12:00:00','12:30:00',
                             #'13:00:00','13:30:00',
                             #'14:00:00','14:30:00',
                             #'15:00:00','15:30:00',
                             #'16:00:00','16:30:00',
                             #'17:00:00','17:30:00',
                             #'18:00:00','18:30:00',
                             #'19:00:00','19:30:00',]) 
    
    
    # 10 min minute interval
    #flight_times = np.array(['06:00:00','06:10:00','06:20:00','06:30:00','06:40:00','06:50:00',
                                 #'07:00:00','07:10:00','07:20:00','07:30:00','07:40:00','07:50:00',
                                 #'08:00:00','08:10:00','08:20:00','08:30:00','08:40:00','08:50:00',
                                 #'09:00:00','09:10:00','09:20:00','09:30:00','09:40:00','09:50:00',
                                 #'10:00:00','10:10:00','10:20:00','10:30:00','10:40:00','10:50:00',
                                 #'11:00:00','11:10:00','11:20:00','11:30:00','11:40:00','11:50:00',
                                 #'12:00:00','12:10:00','12:20:00','12:30:00','12:40:00','12:50:00',
                                 #'13:00:00','13:10:00','13:20:00','13:30:00','13:40:00','13:50:00',
                                 #'14:00:00','14:10:00','14:20:00','14:30:00','14:40:00','14:50:00',
                                 #'15:00:00','15:10:00','15:20:00','15:30:00','15:40:00','15:50:00',
                                 #'16:00:00','16:10:00','16:20:00','16:30:00','16:40:00','16:50:00',
                                 #'17:00:00','17:10:00','17:20:00','17:30:00','17:40:00','17:50:00',
                                 #'18:00:00','18:10:00','18:20:00','18:30:00','18:40:00','18:50:00',
                                 #'19:00:00','19:10:00','19:20:00','19:30:00','19:40:00','19:50:00',])     
    
        
    postprocess_direct_flight_noise(city,city_acronym,aircraft_models,altitudes,dep_loc,des_loc,flight_times,jetliner_L_AeqT_24hr,
                                    N_lat,N_long,long_dist,lat_dist,long_deg,lat_deg,elevation) 
    
    #flight_times = np.array(['00:00:00','00:10:00','00:20:00','00:30:00','00:40:00','00:50:00'])     
    
    #postprocess_direct_flight_noise_2hr(city,city_acronym,aircraft_models,altitudes,dep_loc,des_loc,flight_times,jetliner_L_AeqT_24hr,
    #                                N_lat,N_long,long_dist,lat_dist,long_deg,lat_deg,elevation)     
    

    tf = time.time() 
    print ('time taken: '+ str(round(((tf-ti)/60),3)) + ' mins')     
    return 

def postprocess_direct_flight_noise_2hr(city,city_acronym,aircraft_models,altitudes,dep_loc,des_loc,flight_times,jetliner_L_AeqT_24hr,
                                    N_lat,N_long,long_dist,lat_dist,long_deg,lat_deg,elevation):   
    
    Results           = Data() 
    Results.lat_deg   = lat_deg
    Results.long_deg  = long_deg 
    Results.lat_dist  = lat_dist 
    Results.long_dist = long_dist 
    Results.elevation = elevation 
    
    time_step         = 2 # seconds  
    number_of_flights = len(flight_times) 
    T                 = 2*Units.hours
    num_time_steps    = int(T/time_step) 
    
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):
            altitude      = altitudes[alt_i]
            aircraft      = aircraft_models[m_i]  
            
            # get aircraft noise and annotate 
            L_Amax            = np.zeros((len(dep_loc) ,N_lat,N_long))     # maximum A-Weighted sound pressure level 
            CME               = np.ones((num_time_steps,N_lat,N_long))*np.nan   # cumulative noise exposure
            routes            = []
            mission_times     = np.zeros(len(dep_loc))

            for flight_no in range(len(dep_loc)):  
                
                # create file name path
                filename          = 'Raw_Data_' + aircraft + '/' + aircraft + '_' + altitude + 'ft_mission_' + city_acronym + '_' +  dep_loc[flight_no] + '_to_' + des_loc[flight_no] 
                results           = load_results(filename)         
                 
                # post process noise data  
                noise_data  = post_process_noise_data(results,time_step)   
                t           = noise_data.time                     
                time_step   = t[1]-t[0]  
                for i in range(number_of_flights): 
                    # get start time of flight
                    t0  = int((np.float(flight_times[i].split(':')[0])*60*60 + \
                              np.float(flight_times[i].split(':')[1])*60 + \
                              np.float(flight_times[i].split(':')[2]))/time_step)  
                    p_prefs_A                   =  10**(CME[t0:t0+len(t)][:,:,:,None]/10)
                    p_prefs_B                   =  10**(noise_data.SPL_dBA[:,:,:,None]/10)
                    C                           =  np.concatenate((p_prefs_A,p_prefs_B),axis = 3)
                    CME[t0:t0+len(t)]           =  10*np.log10(np.nansum(C,axis=3))  
                
    
                # 1. L_Amax
                L_Amax[flight_no] = np.max(noise_data.SPL_dBA,axis=0)  

                routes.append(dep_loc[flight_no] + '_to_' + des_loc[flight_no])  
                mission_times[flight_no]      = results.segments[-1].conditions.frames.inertial.time[-1,0]
            
            Results.cumulative_noise_exposure = CME
            frequency                         = str(int(np.float(flight_times[1].split(':')[1]) - np.float(flight_times[0].split(':')[1]))) 
            if frequency == '0':
                frequency = '60' 
            
            processed_results_filename        = 'Raw_Data_' + aircraft + '/'+ aircraft + '_' + altitude + 'ft_' + city_acronym + '_' + frequency +  '2_hrs' 
            save_results(Results,processed_results_filename) 
    
    return  



 
def postprocess_direct_flight_noise(city,city_acronym,aircraft_models,altitudes,dep_loc,des_loc,flight_times,jetliner_L_AeqT_24hr,
                                    N_lat,N_long,long_dist,lat_dist,long_deg,lat_deg,elevation):   
    
    Results           = Data() 
    Results.lat_deg   = lat_deg
    Results.long_deg  = long_deg 
    Results.lat_dist  = lat_dist 
    Results.long_dist = long_dist 
    Results.elevation = elevation 
    
    time_step         = 20 # seconds  
    number_of_flights = len(flight_times) 
    T                 = 15*Units.hours
    num_time_steps    = int(T/time_step) 
    
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):
            altitude      = altitudes[alt_i]
            aircraft      = aircraft_models[m_i]  
            
            # get aircraft noise and annotate 
            L_Amax            = np.zeros((len(dep_loc) ,N_lat,N_long))     # maximum A-Weighted sound pressure level 
            CME               = np.ones((num_time_steps,N_lat,N_long))*np.nan   # cumulative noise exposure
            routes            = []
            mission_times     = np.zeros(len(dep_loc))
            energy_per_flight = np.zeros(len(dep_loc))
            energy_per_route  = np.zeros(len(dep_loc))

            for flight_no in range(len(dep_loc)):  
                
                # create file name path
                filename          = 'Raw_Data_' + aircraft + '/' + aircraft + '_' + altitude + 'ft_mission_' + city_acronym + '_' +  dep_loc[flight_no] + '_to_' + des_loc[flight_no] 
                results           = load_results(filename)         
                # post process noise data  
                noise_data  = post_process_noise_data(results,time_step)   
                t           = noise_data.time                     
                time_step   = t[1]-t[0]  
                for i in range(number_of_flights): 
                    # get start time of flight
                    t0  = int((np.float(flight_times[i].split(':')[0])*60*60 + \
                              np.float(flight_times[i].split(':')[1])*60 + \
                              np.float(flight_times[i].split(':')[2]) - 6*Units.hours)/time_step)  
                    p_prefs_A                   =  10**(CME[t0:t0+len(t)][:,:,:,None]/10)
                    noise_data.SPL_dBA[noise_data.SPL_dBA == 28.] =  np.nan
                    p_prefs_B                   =  10**(noise_data.SPL_dBA[:,:,:,None]/10)
                    C                           =  np.concatenate((p_prefs_A,p_prefs_B),axis = 3)
                    CME[t0:t0+len(t)]           =  10*np.log10(np.nansum(C,axis=3))  
                
    
                # 1. L_Amax
                L_Amax[flight_no] = np.max(noise_data.SPL_dBA,axis=0)  

                routes.append(dep_loc[flight_no] + '_to_' + des_loc[flight_no])  
                mission_times[flight_no]      = results.segments[-1].conditions.frames.inertial.time[-1,0]
                energy_per_flight[flight_no]  = results.segments[0].conditions.propulsion.battery.pack.energy[0,0] - results.segments[-1].conditions.propulsion.battery.pack.energy[-1,0]
                energy_per_route[flight_no]   = energy_per_flight[flight_no]*number_of_flights
            
            Results.UAM_routes                  = routes 
            Results.UAM_route_times             = mission_times
            Results.energy_per_flight_per_route = energy_per_flight
            Results.energy_per_route            = energy_per_route
            Results.total_energy_consumed       = np.sum(energy_per_route)
            
            Results.cumulative_noise_exposure = CME
                            
            # 2. L_AeqT_24hr  
            delta_t                = time_step*np.ones((num_time_steps,N_lat,N_long))   
            L_AeqT_24hr            = np.zeros_like(CME)
            L_AeqT_24hr[:,:,:]     = CME   
            p_i                    = 10**(L_AeqT_24hr/10)   
            Results.L_AeqT_24hr    = 10*np.log10((1/(24*Units.hours))*np.nansum(p_i*delta_t, axis = 0)) 
            
            # 2. L_AeqT
            delta_t                = time_step*np.ones((num_time_steps,N_lat,N_long))   
            L_AeqT                 = np.zeros_like(CME)
            L_AeqT[:,:,:]          = CME   
            p_i                    = 10**(L_AeqT/10)   
            Results.L_AeqT         = 10*np.log10((1/T)*np.nansum(p_i*delta_t, axis = 0))             
            
            # 3. L_AeqT  
            Results.SEL            =  SPL_arithmetic(CME,sum_axis=0)
                 
            # 4. L_AeqT 
            idx_7am                = int(1*Units.hours/time_step)  
            DNL                    = np.zeros_like(CME)
            DNL[:,:,:]             = CME
            DNL[0:idx_7am]         = DNL[0:idx_7am] + 10   
            p_dn_i                 = 10**(DNL/10)    
            Results.L_dn           = 10*np.log10((1/(24*Units.hours))*np.nansum(p_dn_i*delta_t, axis = 0))   
                            
            Results.L_Aeq_jetliner    = jetliner_L_AeqT_24hr  
            Results.L_AeqT_24hr_total = 10*np.log10(np.nansum(np.concatenate((10**(jetliner_L_AeqT_24hr[None,:,:]/10),10**(Results.L_AeqT_24hr[None,:,:]/10))),axis=0))    
            Results.L_Amax            = np.max(L_Amax,axis=0)  
            
            jetliner_band_45_50    = count_between(jetliner_L_AeqT_24hr, 45,50)  
            jetliner_band_50_55    = count_between(jetliner_L_AeqT_24hr, 50,55)  
            jetliner_band_55_60    = count_between(jetliner_L_AeqT_24hr, 55,60)  
            jetliner_band_60_70    = count_between(jetliner_L_AeqT_24hr, 60,70) 
            jetliner_band_70_80    = count_between(jetliner_L_AeqT_24hr, 70,80) 
            jetliner_band_80_90    = count_between(jetliner_L_AeqT_24hr, 80,90)   
            total_band_45_50       = count_between(Results.L_AeqT_24hr_total, 45,50)   
            total_band_50_55       = count_between(Results.L_AeqT_24hr_total, 50,55)  
            total_band_55_60       = count_between(Results.L_AeqT_24hr_total, 55,60)  
            total_band_60_70       = count_between(Results.L_AeqT_24hr_total, 60,70) 
            total_band_70_80       = count_between(Results.L_AeqT_24hr_total, 70,80)  
            total_band_80_90       = count_between(Results.L_AeqT_24hr_total, 80,90)     
            
            # percent increases 
            Results.increase_band_45_50  = 100*(total_band_45_50 - jetliner_band_45_50    )/jetliner_band_45_50   
            Results.increase_band_50_55  = 100*(total_band_50_55 - jetliner_band_50_55    )/jetliner_band_50_55   
            Results.increase_band_55_60  = 100*(total_band_55_60 - jetliner_band_55_60    )/jetliner_band_55_60   
            Results.increase_band_60_70  = 100*(total_band_60_70 - jetliner_band_60_70    )/jetliner_band_60_70  
            Results.increase_band_70_80  = 100*(total_band_70_80 - jetliner_band_70_80    )/jetliner_band_70_80   
            Results.increase_band_80_90  = 100*(total_band_80_90 - jetliner_band_80_90    )/jetliner_band_80_90     
             
            frequency                         = str(int(np.float(flight_times[1].split(':')[1]) - np.float(flight_times[0].split(':')[1]))) 
            if frequency == '0':
                frequency = '60' 
            
            processed_results_filename        = 'Raw_Data_' + aircraft + '/'+ aircraft + '_' + altitude + 'ft_' + city_acronym + '_' + frequency +  'min_All' 
            save_results(Results,processed_results_filename) 
    
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

def community_noise_map(total_noise_file,car_noise_file,train_noise_file, aircraft_noise_file,long_dist,lat_dist,long_deg,lat_deg,N_lat,N_long):
    
    # import image 
    #total_noise_raw_data          = imread(total_noise_file ) 
    car_noise_raw_data            = imread(car_noise_file )
    train_noise_raw_data          = imread(train_noise_file )
    aircraft_noise_raw_data       = imread(aircraft_noise_file )
    
    # import test samples of colors 
    blue            = imread('../Maps_and_Scales/Color_Scales/Color_Blue.png' )    #  > 90
    magenta         = imread('../Maps_and_Scales/Color_Scales/Color_Magenta.png' ) #  80 - 89.0
    purple          = imread('../Maps_and_Scales/Color_Scales/Color_Purple.png' )  #  70 - 79.9
    pink            = imread('../Maps_and_Scales/Color_Scales/Color_Pink.png' )    #  60 - 69.9 
    red             = imread('../Maps_and_Scales/Color_Scales/Color_Red.png' )     #  55 - 59.9
    orange          = imread('../Maps_and_Scales/Color_Scales/Color_Orange.png')   #  50 - 54.9
    yellow_scale    = imread('../Maps_and_Scales/Color_Scales/Color_Yellow.png')   #  45 - 49.9
    grey_background = imread('../Maps_and_Scales/Color_Scales/Color_Gray.png' )    # background
    black_sea       = imread('../Maps_and_Scales/Color_Scales/Color_Black.png' )   # sea  
    
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
    
     
    airport_noise  = vectorize_image(aircraft_noise_raw_data ,color_scales,N_lat,N_long) 
    car_noise      = vectorize_image(car_noise_raw_data ,color_scales,N_lat,N_long) 
    train_noise    = vectorize_image(train_noise_raw_data ,color_scales,N_lat,N_long)  
    
    return airport_noise, car_noise, train_noise

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


def count_between(X, lower_bound, upper_bound):
    count = 0
    x = X.flatten()
    for i in range(len(x)):
        if x[i] >= lower_bound and x[i] < upper_bound: 
            count += 1  
    return count   

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


