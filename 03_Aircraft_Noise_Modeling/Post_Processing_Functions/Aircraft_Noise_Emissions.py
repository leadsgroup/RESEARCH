# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ----------------------------------------------------------------------------------------------------------------------  
# noise imports
import RCAIDE
from RCAIDE.Framework.Core import Data , Units 
from RCAIDE.Library.Methods.Noise.Common.background_noise     import background_noise
from RCAIDE.Library.Methods.Noise.Metrics import * 
from RCAIDE.Library.Methods.Noise.Common.generate_zero_elevation_microphone_locations import generate_zero_elevation_microphone_locations 
from RCAIDE.Library.Methods.Noise.Common.generate_terrain_microphone_locations        import generate_terrain_microphone_locations     
from RCAIDE.Library.Methods.Noise.Common.compute_relative_noise_evaluation_locations  import compute_relative_noise_evaluation_locations
from RCAIDE.Library.Methods.Geodesics.compute_point_to_point_geospacial_data          import compute_point_to_point_geospacial_data   
from RCAIDE.Library.Methods.Noise.Common.background_noise     import background_noise  

# Python package imports   
import numpy as np  
from RCAIDE.Library.Plots import * 

# package imports
import numpy as np
from scipy.interpolate                                           import RegularGridInterpolator  
import os 
import pickle
import sys 
import pandas as pd 
from RCAIDE import  load 
from RCAIDE import  save  
# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def post_process_noise_data(results, topography_file ,  flight_times, time_period , noise_times_steps,  N_gm_x, N_gm_y): 
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
    settings   = results.segments[0].analyses.noise.settings
    n          = settings.number_of_microphone_in_stencil 
    noise_data = Data()   
     
    # Step 2: Determing microhpone points where noise is to be computed
    microphone_coordinates   = None
    settings.topography_file = topography_file
    if settings.topography_file !=  None:
        compute_point_to_point_geospacial_data(settings)
        microphone_locations ,microphone_coordinates = generate_terrain_microphone_locations(settings) 
        noise_data.topography_file                    = settings.topography_file  
        noise_data.microphone_coordinates             = microphone_coordinates.reshape(N_gm_x,N_gm_y,3)     
        noise_data.aircraft_origin_coordinates        = settings.aircraft_origin_coordinates          
        noise_data.aircraft_destination_coordinates   = settings.aircraft_destination_coordinates         
    else:    
        microphone_locations =  generate_zero_elevation_microphone_locations(settings)   
    noise_data.microphone_y_resolution       = N_gm_y
    noise_data.microphone_x_resolution       = N_gm_x              
    noise_data.microphone_locations          = microphone_locations.reshape(N_gm_x,N_gm_y,3)         
    
    # Step 3: Create empty arrays to store noise data 
    num_fligth_segs = len(results.segments)
    num_gm_mic      = len(microphone_locations) 
    
    # Step 4: Initalize Arrays 
    N_ctrl_pts            = ( num_fligth_segs-1) * (noise_times_steps -1) + noise_times_steps # ensures that noise is computed continuously across segments 
    SPL_dBA               = np.ones((N_ctrl_pts,N_gm_x,N_gm_y))*background_noise()  
    Aircraft_pos          = np.empty((0,3))
    Time                  = np.empty((0))
    mic_locs              = np.zeros((N_ctrl_pts,n))   
 
    idx =  0
    
    # Step 5: loop through segments and store noise 
    for seg in range(num_fligth_segs):  
        segment    = results.segments[seg]
        settings   = segment.analyses.noise.settings  
        phi        = settings.noise_hemisphere_phi_angles
        theta      = settings.noise_hemisphere_theta_angles
        conditions = segment.state.conditions  
        time       = conditions.frames.inertial.time[:,0]
        
        # Step 5.1 : Compute relative microhpone locations 
        noise_time,noise_pos,RML,PHI,THETA,num_gm_mic  = compute_relative_noise_evaluation_locations(settings, microphone_locations,segment) 
         
        # Step 5.2: Compute aircraft position and npose at interpolated hemisphere locations
        cpt   = 0 
        if seg == (num_fligth_segs - 1):
            noise_time_ = noise_time 
        else:
            noise_time_ = noise_time[:-1]
             
        Aircraft_pos = np.vstack((Aircraft_pos,noise_pos))
        Time         = np.hstack((Time,noise_time_))
        
        for i in range(len(noise_time_)):
            # Step 5.2.1 :Noise interpolation 
            delta_t         = (noise_time[i] -time[cpt]) / (time[cpt+1] - time[cpt])
            SPL_lower       = conditions.noise.hemisphere_SPL_dBA[cpt].reshape(len(phi),len(theta))
            SPL_uppper      = conditions.noise.hemisphere_SPL_dBA[cpt+1].reshape(len(phi),len(theta))
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
            R_ref                = settings.noise_hemisphere_radius  
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
    L_eq,L_eq_24hr, L_dn =  compute_noise_metrics(noise_data, flight_times)
    
    processed_noise_data = Data(
         L_max      = L_max,
         L_eq       = L_eq,
         L_eq_24hr  = L_eq_24hr,
         L_dn       = L_dn)
    return  processed_noise_data

    
# ----------------------------------------------------------------------------------------------------------------------  
#  compute_noise_metrics
# ----------------------------------------------------------------------------------------------------------------------     
## @ingroup Methods-Noise-Metrics  
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
 