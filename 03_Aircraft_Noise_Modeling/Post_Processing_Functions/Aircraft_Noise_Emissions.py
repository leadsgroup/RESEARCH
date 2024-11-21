# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  
# noise imports
import RCAIDE
from RCAIDE.Framework.Core import Data , Units 
from RCAIDE.Library.Methods.Noise.Common.background_noise     import background_noise
from RCAIDE.Library.Methods.Noise.Metrics                     import *    
from RCAIDE.Library.Methods.Noise.Common.background_noise     import background_noise 
from RCAIDE.Framework.Analyses.Geodesics.Geodesics            import Calculate_Distance
from RCAIDE.Library.Plots import * 

# package imports 
from scipy.interpolate import griddata
from scipy.interpolate  import RegularGridInterpolator     
import numpy as np  
 


def read_flight_simulation_results(results, baseline_results, aircraft_origin_coordinates, aircraft_destination_coordinates):
    # Unpack settings      
    settings   = baseline_results.segments[0].analyses.noise.settings
 
    # Create data strutures for storing noise data 
    noise_results          =  Data()
    noise_results.settings =  Data() 
    noise_results.settings.aircraft_origin_coordinates      = aircraft_origin_coordinates      
    noise_results.settings.aircraft_destination_coordinates = aircraft_destination_coordinates  
    noise_results.settings.number_of_microphone_in_stencil  = settings.number_of_microphone_in_stencil 
    noise_results.settings.noise_hemisphere_phi_angles      = settings.noise_hemisphere_phi_angles
    noise_results.settings.noise_hemisphere_theta_angles    = settings.noise_hemisphere_theta_angles
    noise_results.settings.mean_sea_level_altitude          = settings.mean_sea_level_altitude
    noise_results.settings.noise_hemisphere_radius          = settings.noise_hemisphere_radius  
                           
    
    # Step 3: Create empty arrays to store noise data 
    N_segs = len(results.segments) 
    N_cpts = results.segments[0].state.numerics.number_of_control_points

    noise_results.segment_name            = []
    noise_results.time                    = np.zeros((N_segs,N_cpts,1 )) 
    noise_results.position_vector         = np.zeros((N_segs,N_cpts,3 )) 
    noise_results.hemisphere_SPL_dBA      = np.zeros((N_segs,N_cpts, 72)) 
    
    # Step 5: loop through segments and store noise 
    for seg in range(N_segs):  
        segment                                     = results.segments[seg] 
        noise_results.segment_name.append(segment.tag)
        noise_results.time[seg]                     = segment.state.conditions.frames.inertial.time 
        noise_results.position_vector[seg]          = segment.state.conditions.frames.inertial.position_vector
        noise_results.hemisphere_SPL_dBA[seg]       = segment.state.conditions.noise.hemisphere_SPL_dBA 
         
    Total_Energy = 0
    # get total energy comsumed for each route and append it onto noise 
    for network in results.segments[0].analyses.energy.vehicle.networks: 
        busses  = network.busses 
        for bus in busses:
            for battery_module in bus.battery_modules: 
                E_start_module  = results.segments[0].conditions.energy[bus.tag].battery_modules[battery_module.tag].energy 
                E_end_module    = results.segments[-1].conditions.energy[bus.tag].battery_modules[battery_module.tag].energy
                Total_Energy    += (E_start_module - E_end_module)  
    noise_results.energy_consumed = Total_Energy
        
    return noise_results

# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def post_process_noise_data(noise_results, topography_file ,  flight_times, time_period , noise_times_steps,  N_gm_x, N_gm_y): 
    """This translates all noise data into metadata for plotting 
    
    Assumptions:
    None
    
    Source: 
 
    Inputs: results 
         
    Outputs: noise_data
    
    Properties Used:
    N/A
    """
  
    # Step 1: Unpack settings       
    noise_data = Data()   
     
    # Step 2: Determing microhpone points where noise is to be computed 
    noise_results.settings.topography_file  = topography_file 
    aircraft_origin_coordinates             = noise_results.settings.aircraft_origin_coordinates      
    aircraft_destination_coordinates        = noise_results.settings.aircraft_destination_coordinates  
    n                                       = noise_results.settings.number_of_microphone_in_stencil  

    aircraft_origin_location ,  aircraft_destination_location =  compute_point_to_point_geospacial_data(topography_file,aircraft_origin_coordinates,aircraft_destination_coordinates)

    microphone_locations ,microphone_coordinates  = generate_terrain_microphone_locations(topography_file, N_gm_x, N_gm_y) 
    noise_data.topography_file                    = noise_results.settings.topography_file  
    noise_data.microphone_coordinates             = microphone_coordinates.reshape(N_gm_x,N_gm_y,3)     
    noise_data.aircraft_origin_coordinates        = noise_results.settings.aircraft_origin_coordinates          
    noise_data.aircraft_destination_coordinates   = noise_results.settings.aircraft_destination_coordinates     
    noise_data.microphone_y_resolution            = N_gm_y
    noise_data.microphone_x_resolution            = N_gm_x              
    noise_data.microphone_locations               = microphone_locations.reshape(N_gm_x,N_gm_y,3)
    noise_data.energy_consumed                    = noise_results.energy_consumed
    
    # Step 3: Create empty arrays to store noise data 
    N_segs = len(noise_results.time[:, 0, 0])
    num_gm_mic      = len(microphone_locations) 
    
    # Step 4: Initalize Arrays 
    N_ctrl_pts            = ( N_segs-1) * (noise_times_steps -1) + noise_times_steps # ensures that noise is computed continuously across segments 
    SPL_dBA               = np.ones((N_ctrl_pts,N_gm_x,N_gm_y))*background_noise()  
    Aircraft_pos          = np.empty((0,3))
    Time                  = np.empty((0))
    mic_locs              = np.zeros((N_ctrl_pts,n))   
 
    idx =  0
    
    # Step 5: loop through segments and store noise 
    for seg in range(N_segs):   
        phi                      = noise_results.settings.noise_hemisphere_phi_angles
        theta                    = noise_results.settings.noise_hemisphere_theta_angles
        mean_sea_level_altitude  = noise_results.settings.mean_sea_level_altitude
        time                     = noise_results.time[seg][:,0]
        position_vector          = noise_results.position_vector[seg]
        hemisphere_SPL_dBA       = noise_results.hemisphere_SPL_dBA[seg]  
        
        # Step 5.1 : Compute relative microhpone locations 
        noise_time,noise_pos,RML,PHI,THETA,num_gm_mic  = compute_relative_noise_evaluation_locations(mean_sea_level_altitude,noise_times_steps,aircraft_origin_location,microphone_locations,position_vector,time) 
         
        # Step 5.2: Compute aircraft position and npose at interpolated hemisphere locations
        cpt   = 0 
        if seg == (N_segs - 1):
            noise_time_ = noise_time 
        else:
            noise_time_ = noise_time[:-1]
             
        Aircraft_pos = np.vstack((Aircraft_pos,noise_pos))
        Time         = np.hstack((Time,noise_time_))
        
        for i in range(len(noise_time_)):
            # Step 5.2.1 :Noise interpolation 
            delta_t         = (noise_time[i] -time[cpt]) / (time[cpt+1] - time[cpt])
            SPL_lower       = hemisphere_SPL_dBA[cpt].reshape(len(phi),len(theta))
            SPL_uppper      = hemisphere_SPL_dBA[cpt+1].reshape(len(phi),len(theta))
            SPL_gradient    = SPL_uppper -  SPL_lower
            SPL_interp      = SPL_lower + SPL_gradient *delta_t
     
            #  Step 5.2.2 Create surrogate   
            SPL_dBA_surrogate = RegularGridInterpolator((phi, theta),SPL_interp  ,method = 'linear',   bounds_error=False, fill_value=None)       
            
            #  Step 5.2.3 Query surrogate
            R                 = np.linalg.norm(RML[i], axis=1) 
            locs              = np.argsort(R)[:n]
            pts               = (PHI[i][locs],THETA[i][locs]) 
            SPL_dBA_unscaled  = SPL_dBA_surrogate(pts) 
            
            #  Step 5.2.4 Scale data using radius  
            R_ref                = noise_results.settings.noise_hemisphere_radius  
            SPL_dBA_scaled       = SPL_dBA_unscaled - 20*np.log10(R[locs]/R_ref)
            
            SPL_dBA_temp         = SPL_dBA[idx].flatten()
            SPL_dBA_temp[locs]   = SPL_dBA_scaled
            SPL_dBA[idx]         = SPL_dBA_temp.reshape(N_gm_x,N_gm_y) 
            mic_locs[idx]        = locs 
            idx += 1
            
            if noise_time[i] >= time[cpt+1]:
                cpt += 1             
                
    # Step 6: Make any readings less that background noise equal to background noise
    SPL_dBA                             = np.nan_to_num(SPL_dBA) 
    SPL_dBA[SPL_dBA<background_noise()] = background_noise()  
     
    # Step 7: Store data 
    noise_data.SPL_dBA               = SPL_dBA
    noise_data.time                  = Time 
    noise_data.aircraft_position     = Aircraft_pos
    noise_data.microhpone_locations  = mic_locs
    
    L_max = np.max(SPL_dBA, axis= 0) 
    
    # Step 8: Perform noise metric calculations 
    L_eq,L_eq_24hr, L_dn =  compute_noise_metrics(noise_data, flight_times, time_period)
    
    processed_noise_data = Data( 
         energy_consumed = noise_data.energy_consumed,      
         L_max      = L_max,
         L_eq       = L_eq,
         L_eq_24hr  = L_eq_24hr,
         L_dn       = L_dn)
    return  processed_noise_data


# ----------------------------------------------------------------------
#  Compute Point to Point Geospacial Data
# ---------------------------------------------------------------------
## @ingroup Library-Missions
def compute_point_to_point_geospacial_data(topography_file,aircraft_origin_coordinates,aircraft_destination_coordinates):
    """This computes the absolute microphone/observer locations on a defined topography
            
    Assumptions: 
        topography_file is a text file obtained from https://topex.ucsd.edu/cgi-bin/get_data.cgi
    
    Source:
        N/A  

    Inputs:   
        topography_file                        - file of lattide, longitude and elevation points     
        origin_coordinates                     - coordinates of origin location                                              [degrees]
        destination_coordinates                - coordinates of destimation location                                            [degrees]  
        
    Outputs:                                   
        latitude_longitude_micrphone_locations - latitude-longitude and elevation coordinates of all microphones in domain      [deg,deg,m]  
        flight_range                           - gound distance between origin and destination location                      [meters]              
        true_course                            - true course angle measured clockwise from true north                     [radians]                      
        origin_location                        - cartesial coordinates of origin location relative to computational domain   [meters]                   
        destination_xyz_location               - cartesial coordinates of destination location relative to computational domain [meters]    
    
    Properties Used:
        N/A       
    """     
    # convert cooordinates to array 
    origin_coordinates      = np.asarray(aircraft_origin_coordinates)
    destination_coordinates = np.asarray(aircraft_destination_coordinates)
    
    # extract data from file 
    data  = np.loadtxt(topography_file)
    Long  = data[:,0]
    Lat   = data[:,1]
    Elev  = data[:,2] 

    x_min_coord = np.min(Lat)
    y_min_coord = np.min(Long)
    dep_lat     = origin_coordinates[0]
    dep_long    = origin_coordinates[1]
    des_lat     = destination_coordinates[0]
    des_long    = destination_coordinates[1]
    if dep_long < 0: 
        dep_long = 360 + dep_long
    if des_long< 0:
        des_long =360 +  des_long 
    
    bottom_left_map_coords   = np.array([x_min_coord,y_min_coord])  
    x0_coord                 = np.array([dep_lat,y_min_coord])
    y0_coord                 = np.array([x_min_coord,dep_long])
    x1_coord                 = np.array([des_lat,y_min_coord])
    y1_coord                 = np.array([x_min_coord,des_long])  
    
    x0 = RCAIDE.Framework.Analyses.Geodesics.Geodesics.Calculate_Distance(x0_coord,bottom_left_map_coords) * Units.kilometers
    y0 = RCAIDE.Framework.Analyses.Geodesics.Geodesics.Calculate_Distance(y0_coord,bottom_left_map_coords) * Units.kilometers
    x1 = RCAIDE.Framework.Analyses.Geodesics.Geodesics.Calculate_Distance(x1_coord,bottom_left_map_coords) * Units.kilometers
    y1 = RCAIDE.Framework.Analyses.Geodesics.Geodesics.Calculate_Distance(y1_coord,bottom_left_map_coords) * Units.kilometers
    
    lat_flag             = np.where(origin_coordinates<0)[0]
    origin_coordinates[lat_flag]  = origin_coordinates[lat_flag] + 360 
    long_flag            = np.where(destination_coordinates<0)[0]
    destination_coordinates[long_flag] = destination_coordinates[long_flag] + 360 
    z0                   = griddata((Lat,Long), Elev, (np.array([origin_coordinates[0]]),np.array([origin_coordinates[1]])), method='nearest')[0]
    z1                   = griddata((Lat,Long), Elev, (np.array([destination_coordinates[0]]),np.array([destination_coordinates[1]])), method='nearest')[0] 
    dep_loc              = np.array([x0,y0,z0])
    des_loc              = np.array([x1,y1,z1])
    
    # pack data 
    aircraft_origin_location      = dep_loc
    aircraft_destination_location = des_loc 
        
    return aircraft_origin_location ,  aircraft_destination_location


# ---------------------------------------------------------------------------------------------------------------------- 
#  generate_terrain_microphone_locations
# ----------------------------------------------------------------------------------------------------------------------  
def generate_terrain_microphone_locations(topography_file, microphone_x_resolution, microphone_y_resolution):
    """This computes the absolute microphone/observer locations on a defined topography
            
    Assumptions: 
        topography_file is a text file obtained from https://topex.ucsd.edu/cgi-bin/get_data.cgi 
    """      
    # extract data from file 
    data  = np.loadtxt(topography_file) # settings.topography_file) CHANGED 10-15-2024
    Long  = data[:,0]
    Lat   = data[:,1]
    Elev  = data[:,2] 
    
    x_min_coord = np.min(Lat)
    x_max_coord = np.max(Lat)
    y_min_coord = np.min(Long)
    y_max_coord = np.max(Long) 

    if y_min_coord < 0: 
        y_min_coord = 360 + y_min_coord
    if y_max_coord< 0:
        y_max_coord=360 + y_max_coord 
    
    top_left_map_coords      = np.array([x_max_coord,y_min_coord])
    bottom_left_map_coords   = np.array([x_min_coord,y_min_coord])  
    bottom_right_map_coords  = np.array([x_min_coord,y_max_coord]) 
    
    x_dist_max = Calculate_Distance(top_left_map_coords,bottom_left_map_coords) * Units.kilometers
    y_dist_max = Calculate_Distance(bottom_right_map_coords,bottom_left_map_coords) * Units.kilometers
    
    [y_pts,x_pts]      = np.meshgrid(np.linspace(0,y_dist_max,microphone_y_resolution),np.linspace(0,x_dist_max,microphone_x_resolution))
    [long_deg,lat_deg] = np.meshgrid(np.linspace(np.min(Long),np.max(Long),microphone_y_resolution),np.linspace(np.min(Lat),np.max(Lat),microphone_x_resolution)) 
    z_deg              = griddata((Lat,Long), Elev, (lat_deg, long_deg), method='linear')        
    cartesian_pts      = np.dstack((np.dstack((x_pts[:,:,None],y_pts[:,:,None] )),z_deg[:,:,None])).reshape(microphone_x_resolution*microphone_y_resolution,3)
    lat_long_pts       = np.dstack((np.dstack((lat_deg[:,:,None],long_deg[:,:,None] )),z_deg[:,:,None])).reshape(microphone_x_resolution*microphone_y_resolution,3)  
    return cartesian_pts , lat_long_pts

  
# ----------------------------------------------------------------------------------------------------------------------  
#  Relative Noise Evaluatation Locations
# ----------------------------------------------------------------------------------------------------------------------       
def compute_relative_noise_evaluation_locations(mean_sea_level_altitude,noise_times_steps, aircraft_origin_location,microphone_locations,position_vector,time):
    """This computes the relative locations on the surface in the computational domain where the 
    propogated sound is computed. Vectors point from observer/microphone to aircraft/source  
            
    Assumptions: 
        Acoustic scattering is not modeled

    Source:
        N/A  

    Inputs:  
        settings.microphone_locations                - array of microphone locations on the ground  [meters] 
        segment.conditions.frames.inertial.position_vector  - position of aircraft                         [boolean]                                          

    Outputs: 
    GM_THETA   - angle measured from ground microphone in the x-z plane from microphone to aircraft 
    GM_PHI     - angle measured from ground microphone in the y-z plane from microphone to aircraft 
    RML        - relative microphone locations  
    num_gm_mic - number of ground microphones
 
    Properties Used:
        N/A       
    """       
  
    MSL_altitude      = mean_sea_level_altitude
    N                 = noise_times_steps
    
    # rediscretize time and aircraft position to get finer resolution  
    noise_time        = np.linspace(time[0], time[-1], N) 
    noise_pos         = np.zeros((N,3)) 
    noise_pos[:,0]    = np.interp(noise_time,time,position_vector[:,0])
    noise_pos[:,1]    = np.interp(noise_time,time,position_vector[:,1])
    noise_pos[:,2]    = np.interp(noise_time,time,position_vector[:,2])
    
    num_gm_mic        = len(microphone_locations)  
    RML               = np.zeros((N,num_gm_mic,3)) 
    PHI               = np.zeros((N,num_gm_mic))
    THETA             = np.zeros((N,num_gm_mic)) 
    
    for cpt in range(N):  
        relative_locations         = np.zeros((num_gm_mic,3))
        relative_locations[:,0]    = microphone_locations[:,0] - (aircraft_origin_location[0] + noise_pos[cpt,0])    
        relative_locations[:,1]    = microphone_locations[:,1] - (aircraft_origin_location[1] + noise_pos[cpt,1]) 
        if MSL_altitude:
            relative_locations[:,2]    = -(noise_pos[cpt,2])  - microphone_locations[:,2] 
        else:
            relative_locations[:,2]    = -(noise_pos[cpt,2])
            
        RML[cpt,:,:]   = relative_locations 
        PHI[cpt,:]     =  np.arctan2(np.sqrt(np.square(relative_locations[:, 0]) + np.square(relative_locations[:, 1])),  relative_locations[:, 2])  
        THETA[cpt,:]   =  np.arctan2(relative_locations[:, 1], relative_locations[:, 0]) 
    
    return noise_time,noise_pos,RML,PHI,THETA,num_gm_mic 
  
  
      
# ----------------------------------------------------------------------------------------------------------------------  
#  compute_noise_metrics
# ----------------------------------------------------------------------------------------------------------------------     
def compute_noise_metrics(noise_data, flight_times,time_period):   
    """This method calculates the Average A-weighted Sound Level (LAeqT), the Day-Night Average Sound Level and the
    Single Event Noise Exposure Level (SENEL) or Sound Exposure Level (SEL)
    
    Assumptions:  
        Flights occure between 6:00 and 9:00 pm (i.e. a 15 hour window)

    Source:
        None

    Inputs:
       noise_data  - post-processed noise data structure 

    Outputs: [dB]
       noise_data  - post-processed noise data structure 

    Properties Used:
        N/A  
    """
    
    # determine start and end of time period

    DNL_time_period           = ['07:00:00','20:00:00']
    t_7am                     = float(DNL_time_period[0].split(':')[0])*60*60 + float(DNL_time_period[0].split(':')[1])*60 +  float(DNL_time_period[0].split(':')[2])
    t_10pm                    = float(DNL_time_period[1].split(':')[0])*60*60 + float(DNL_time_period[1].split(':')[1])*60 +  float(DNL_time_period[1].split(':')[2])            
    t_start                   = float(time_period[0].split(':')[0])*60*60 + float(time_period[0].split(':')[1])*60 +  float(time_period[0].split(':')[2])
    t_end                     = float(time_period[1].split(':')[0])*60*60 + float(time_period[1].split(':')[1])*60 +  float(time_period[1].split(':')[2])     
    SPL                       = noise_data.SPL_dBA    
    N_gm_y                    = noise_data.microphone_y_resolution   
    N_gm_x                    = noise_data.microphone_x_resolution 
    flight_time               = noise_data.time    
    time_step                 = flight_time[1]-flight_time[0] 
    number_of_flights         = len(flight_times)    
    p_div_p_ref_sq_L_eq       = np.zeros((N_gm_x,N_gm_y)) 
    p_div_p_ref_sq_L_24hr     = np.zeros((N_gm_x,N_gm_y))
    p_div_p_ref_sq_L_dn       = np.zeros((N_gm_x,N_gm_y)) 
    
    ambient_noise_duration      = t_end - t_start  
    ambient_noise_duration_24hr = 24 * Units.hrs 
    for i in range(number_of_flights):   
        t_flight_during_day  = float(flight_times[i].split(':')[0])*60*60 +  float(flight_times[i].split(':')[1])*60 +  float(flight_times[i].split(':')[2]) + flight_time
        
        # compute ambient noise duration 
        ambient_noise_duration      -= flight_time[-1]
        ambient_noise_duration_24hr -= flight_time[-1]

        # create noise penalty 
        noise_penality          = np.zeros((len(flight_time),N_gm_x,N_gm_y))         
        noise_penality[t_flight_during_day<t_7am]  = 10
        noise_penality[t_flight_during_day>t_10pm] = 10         
        
        # convert SPL to pressure and multiply by duration 
        p_sq_ref_flight_sq     = np.nansum(time_step * (10**(SPL/10)), axis=0)   
        p_sq_ref_flight_sq_dn  = np.nansum(time_step * (10**( (noise_penality + SPL)/10)), axis=0)  
    
        # add to current  
        p_div_p_ref_sq_L_eq    = np.nansum(np.concatenate((p_sq_ref_flight_sq[:,:,None],p_div_p_ref_sq_L_eq[:,:,None]),axis = 2), axis =2)
        p_div_p_ref_sq_L_24hr  = np.nansum(np.concatenate((p_sq_ref_flight_sq[:,:,None],p_div_p_ref_sq_L_24hr[:,:,None]),axis = 2), axis =2)
        p_div_p_ref_sq_L_dn    = np.nansum(np.concatenate((p_sq_ref_flight_sq_dn[:,:,None],p_div_p_ref_sq_L_dn[:,:,None]),axis = 2), axis =2) 
    
    # add on background noise for remainder of time
    ambient_noise_duration       = np.maximum(ambient_noise_duration,0)
    ambient_noise_duration_24hr  = np.maximum(ambient_noise_duration_24hr,0)
    p_div_p_ref_sq_L_eq         += ambient_noise_duration *  (10**(background_noise()/10))
    p_div_p_ref_sq_L_24hr       += ambient_noise_duration_24hr *  (10**(background_noise()/10))
    p_div_p_ref_sq_L_dn         += ambient_noise_duration *  (10**(background_noise()/10))
    
    L_eq              = 10*np.log10((1/(t_end-t_start))*p_div_p_ref_sq_L_eq)
    L_eq_24hr         = 10*np.log10((1/(24*Units.hours))*p_div_p_ref_sq_L_24hr)   
    L_dn              = 10*np.log10((1/(t_end-t_start))*p_div_p_ref_sq_L_dn)
        
    return L_eq,L_eq_24hr, L_dn
 