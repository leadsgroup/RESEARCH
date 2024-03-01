# RCAIDE Imports 
from  RCAIDE.Core import Data, Units  

# Python Imports  
import numpy as np   
import pickle 
import matplotlib.pyplot as plt  
import matplotlib.cm as cm 
from RCAIDE.Visualization.Noise.post_process_noise_data import post_process_noise_data

import sys 
 
# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():  

    save_figure = True
    file_type   = '.png'
    
    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 28,
                  'xtick.labelsize': 24,
                  'ytick.labelsize': 24,
                  'axes.titlesize': 28}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 2 
    plot_parameters.line_style       = '-' 
    plot_parameters.figure_width     = 8
    plot_parameters.figure_height    = 6
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20 
    plot_parameters.plot_grid        = True   
    plot_parameters.markers          = ['o','v','s','^','p','^','D','X','*']
    plot_parameters.colors           = cm.inferno(np.linspace(0,1,5))                
    plot_parameters.legend_font      = 20                             # legend_font_size          
    
    #flyover_noise_LAmax(plot_parameters,save_figure,file_type)
    flyover_noise_LAeqT(plot_parameters,save_figure,file_type)
    return 


def flyover_noise_LAmax(PP,save_figure,file_type): 
    idx            = 5
    aircraft_tags  = ['SR','TR','TW','HC'] 
    altitudes      = np.array([200,500,1000,1500,2000])
    PP.colors      = cm.inferno(np.linspace(0.2,0.8,5))[::-1]     
    num_aircrafts  = len(aircraft_tags)
     
      
    car_distances       = np.array([0, 50, 71, 111,130,180,230,500,750,1000])  
    
    for a_i in range(num_aircrafts):
        
        # initialize plots
        fig_name  = aircraft_tags[a_i]  + '_Flyover_Noise' 
        fig = plt.figure(fig_name)    
        fig.set_size_inches(PP.figure_width,PP.figure_height)  
        axes = fig.add_subplot(1,1,1)       
        
        for alt_i in range(len(altitudes)): 
        
            # load results 
            filename     = 'Raw_Data_' + aircraft_tags[a_i]  + '/'+   aircraft_tags[a_i]  + '_Flyover_' + str(altitudes[alt_i]) + 'ft' 
            results      = load_results(filename) 
            
             
            # format results 
            noise_data   = post_process_noise_data(results) 
            N_gm_y       = noise_data.N_gm_y
            SPL          = noise_data.SPL_dBA      
            gm           = noise_data.receptor_locations
            gm_y         = gm[:,:,1]  
            max_SPL      = np.max(SPL,axis=0)          
            semi_pts    = int(np.floor(len(gm_y[idx,:])/2))
            
            # plot results  
            axes.plot(gm_y[idx,semi_pts:]/Units.feet, max_SPL[idx,semi_pts:],markersize = PP.marker_size,  marker = PP.markers[alt_i], linewidth = PP.line_width , color = PP.colors[alt_i], label= str(altitudes[alt_i]) + ' ft' ) 
            #axes.fill_between(gm_y[idx,semi_pts:]/Units.feet, max_SPL[idx,semi_pts:], color = PP.colors[alt_i] , alpha=1.0, label= str(altitudes[alt_i]) + ' ft' )
            
            axes.set_ylabel(r'L$_{Amax}$ [dBA]')
            axes.set_xlabel('Distance from highway [ft]')    
            axes.set_xlim([0,500/Units.feet])  
            axes.set_ylim([32,85])  
    
        car_noise = 1E-04*car_distances**2 - 0.15*car_distances + 67.656 
        #axes.plot(car_distances/Units.feet, np.maximum(car_noise,35), marker = PP.markers[alt_i],markersize = PP.marker_size  , linewidth = PP.line_width , color = 'black', label= 'highway traffic' ) 
        axes.plot(car_distances/Units.feet, np.maximum(car_noise,35) , linewidth = PP.line_width , color = 'black', label= 'highway' )  
        axes.legend(loc='upper right', ncol = 2, prop={'size': PP.legend_font_size})  
        plt.tight_layout()
        # save results 
        if save_figure:
            save_filename  =  '../Papers_and_Presentation/Images/'  + aircraft_tags[a_i]  + '_Flyover_' + str(altitudes[alt_i]) + 'ft' 
            plt.savefig( save_filename+ file_type)   
         
    
    return  



def flyover_noise_LAeqT(PP,save_figure,file_type): 
    

    idx            = 0
    aircraft_tags  = ['SR','TR','TW','HC'] 
    altitudes      = np.array([200,500,1000,1500,2000])
    PP.colors      = cm.inferno(np.linspace(0.2,0.8,5))[::-1]     
    num_aircrafts  = len(aircraft_tags)    
    
    flight_times = np.array(['00:00:00','00:05:00',
                             '00:10:00','00:15:00',
                             '00:20:00','00:25:00',
                             '00:30:00','00:35:00',
                             '00:40:00','00:45:00',
                             '00:50:00','00:55:00', ])   
    time_step         = 20   # seconds  
    number_of_flights = len(flight_times) 
    T                 = 2*Units.hours
    num_time_steps    = int(T/time_step) 
    
    for a_i in range(num_aircrafts):

        # initialize plots
        fig_name  = aircraft_tags[a_i]  + '_Flyover_Noise' 
        fig = plt.figure(fig_name)    
        fig.set_size_inches(PP.figure_width,PP.figure_height)  
        axes = fig.add_subplot(1,1,1)       
        
        for alt_i in range(len(altitudes)): 
            
            # get aircraft noise and annotate  
 
            # create file name path

            # load results 
            filename     = 'Raw_Data_' + aircraft_tags[a_i]  + '/'+   aircraft_tags[a_i]  + '_Flyover_' + str(altitudes[alt_i]) + 'ft' 
            results      = load_results(filename) 
    
            # post process noise data  
            noise_data  = post_process_noise_data(results,time_step)   
            t           = noise_data.time                     
            
            CME         = np.ones((num_time_steps,len(noise_data.receptor_locations[:,0]),len(noise_data.receptor_locations[0,:])))*np.nan   # cumulative noise exposure 
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
            
            # 2. L_AeqT
            delta_t                = time_step*np.ones((num_time_steps,len(noise_data.receptor_locations[:,0]),len(noise_data.receptor_locations[0,:])))   
            L_AeqT                 = np.zeros_like(CME)
            L_AeqT[:,:,:]          = CME   
            p_i                    = 10**(L_AeqT/10)   
            L_AeqT                 = 10*np.log10((1/(1*Units.hours))*np.nansum(p_i[0:3600]*delta_t[0:3600], axis = 0))             
        
            # format results 
            noise_data   = post_process_noise_data(results)   
            gm           = noise_data.receptor_locations
            gm_y         = gm[:,:,1]   
            semi_pts    = int(np.floor(len(gm_y[idx,:])/2))
            
            # plot results  
            axes.plot(gm_y[idx,semi_pts:]/Units.feet, L_AeqT[idx,semi_pts:],markersize = PP.marker_size,  marker = PP.markers[alt_i], linewidth = PP.line_width , color = PP.colors[alt_i], label= str(altitudes[alt_i]) + ' ft' ) 
            #axes.fill_between(gm_y[idx,semi_pts:]/Units.feet, max_SPL[idx,semi_pts:], color = PP.colors[alt_i] , alpha=1.0, label= str(altitudes[alt_i]) + ' ft' )
            
            axes.set_ylabel(r'1hr - L$_{AeqT}$ [dBA]')
            axes.set_xlabel('Distance from highway [ft]')    
            axes.set_xlim([0,500/Units.feet])  
            axes.set_ylim([20,85])  
        
        highway_noise, distances = get_highway_noise()
        axes.plot(distances /Units.feet, highway_noise , linewidth = PP.line_width , color = 'black', label= 'highway' )  
        axes.legend(loc='upper right', ncol = 2, prop={'size': PP.legend_font_size})  
        plt.tight_layout()
        # save results 
        if save_figure:
            save_filename  =  '../Papers_and_Presentation/Images/'  + aircraft_tags[a_i]  + '_Flyover'
            plt.savefig( save_filename+ file_type)   
         
    
    return  


def get_highway_noise():
    

    D  = np.linspace(15,800,100)    
  
    # A = automobile, MT = medium truck , HT = Heavy truck
    
    # number of vehicles in hour 
    N_A  = 1000
    N_MT = 25
    N_HT = 5
    
    # average speed in km/hr
    S_A  = 20  
    S_MT = 20  
    S_HT = 20  
    
    # reference energy 
    Lmean_E_A  = 38.1*np.log10(S_A) - 2.4
    Lmean_E_MT = 33.9*np.log10(S_MT) + 16.4
    Lmean_E_HT = 24.6*np.log10(S_HT) + 38.5
    
    # site parater 
    alpha = 0  
    
    # sheilding adjustment 
    delta_s = 0 
    
    # refenence distance 
    D_0 = 15
     
    # noise level greater or equal to 15 meters 
    Leq_A_ge_15  = Lmean_E_A  + 10*np.log10((N_A*np.pi*D_0)/(S_A)) + 10*np.log10(D_0/D)**(1 + alpha)  + delta_s - 25 
    Leq_MT_ge_15 = Lmean_E_MT + 10*np.log10((N_MT*np.pi*D_0)/(S_MT)) + 10*np.log10(D_0/D)**(1 + alpha) + delta_s - 25 
    Leq_HT_ge_15 = Lmean_E_HT + 10*np.log10((N_HT*np.pi*D_0)/(S_HT)) + 10*np.log10(D_0/D)**(1 + alpha)   + delta_s -25 
     
    Leq_ge_15   = 10*np.log10(  10**(Leq_A_ge_15/10) + 10**(Leq_MT_ge_15/10)  + 10**(Leq_HT_ge_15/10) ) 
    
    distances = np.hstack((np.array([0,5]),D))
    Leq       = np.hstack((np.array([80,66.18]),Leq_ge_15))
    return Leq ,distances

 
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
 
