# MARC Imports 
from  MARC.Core import Data, Units  

# Python Imports  
import numpy as np   
import pickle 
import matplotlib.pyplot as plt  
import matplotlib.cm as cm 
from MARC.Visualization.Performance.Common.post_process_noise_data import post_process_noise_data

# import propeller/rotors geometries  
import sys
sys.path.append('../Results')
 
# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():  

    save_figure = False 
    
    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 2.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 32,
                  'xtick.labelsize': 28,
                  'ytick.labelsize': 28,
                  'axes.titlesize': 32}
    plt.rcParams.update(parameters)
    plot_parameters                  = Data()
    plot_parameters.line_width       = 3 
    plot_parameters.line_style       = '-' 
    plot_parameters.figure_width     = 10 
    plot_parameters.figure_height    = 7 
    plot_parameters.marker_size      = 10 
    plot_parameters.legend_font_size = 20 
    plot_parameters.plot_grid        = True   
    plot_parameters.markers          = ['o','v','s','^','p','^','D','X','*']
    plot_parameters.colors           = cm.inferno(np.linspace(0,1,5))     
    plot_parameters.lw               = 2                              # line_width                
    plot_parameters.legend_font      = 20                             # legend_font_size          
    
    flyover_noise(plot_parameters,save_figure)
    return 


def flyover_noise(PP,save_figure): 
    idx            = 30
    aircraft_tags  = ['SR','TR','TW','HC'] 
    altitudes      = np.array([200,500,1000,1500,2000])
    PP.colors      = cm.inferno(np.linspace(0.2,0.8,5))[::-1]     
    num_aircrafts  = len(aircraft_tags)
    
    
    # vehicle noise (taken from interstate in colorado )
    # sources https://www.codot.gov/programs/research/pdfs/2005/tnm.pdf
    # page 29 
      
    car_distances       = np.array([0, 50, 71, 111,130,180,230,500,750,1000])  
    
    for aircraft in range(num_aircrafts):
        
        # initialize plots
        fig_name  = aircraft_tags[aircraft]  + '_Flyover' 
        fig = plt.figure(fig_name)    
        fig.set_size_inches(PP.figure_width,PP.figure_height)  
        axes = fig.add_subplot(1,1,1)      
        axes.set_ylabel('SPL')
        axes.set_xlabel('range (nmi)')   
        
        for i in range(len(altitudes)): 
        
            # load results 
            filename     = 'Results/' + aircraft_tags[aircraft]  + '_Flyover_' + str(altitudes[i]) + 'ft' 
            results      = load_results(filename) 
            
            # format results 
            noise_data   = post_process_noise_data(results) 
            N_gm_y       = noise_data.N_gm_y
            SPL          = noise_data.SPL_dBA_ground_mic      
            gm           = noise_data.SPL_dBA_ground_mic_loc  
            gm_y         = gm[:,:,1]  
            max_SPL      = np.max(SPL,axis=0)          
            semi_pts    = int(np.floor(len(gm_y[idx,:])/2))
            
            # plot results  
            axes.plot(gm_y[idx,semi_pts:]/Units.feet, max_SPL[idx,semi_pts:],markersize = PP.marker_size,  marker = PP.markers[i], linewidth = PP.line_width , color = PP.colors[i], label= str(altitudes[i]) + ' ft' ) 
            #axes.fill_between(gm_y[idx,semi_pts:]/Units.feet, max_SPL[idx,semi_pts:], color = PP.colors[i] , alpha=1.0, label= str(altitudes[i]) + ' ft' )
            
            axes.set_ylabel(r'SPL$_{max}$ (dBA)')
            axes.set_xlabel('Distance from Roadway (ft)')    
            axes.set_xlim([0,500/Units.feet])  
            axes.set_ylim([32,85])  
    
        car_noise = 1E-04*car_distances**2 - 0.15*car_distances + 67.656 
        #axes.plot(car_distances/Units.feet, np.maximum(car_noise,35), marker = PP.markers[i],markersize = PP.marker_size  , linewidth = PP.line_width , color = 'black', label= 'highway traffic' ) 
        axes.plot(car_distances/Units.feet, np.maximum(car_noise,35) , linewidth = PP.line_width , color = 'black', label= 'highway' )  
        axes.legend(loc='upper right', prop={'size': PP.legend_font_size})  
        plt.tight_layout()
        # save results 
        if save_figure:
            plt.savefig(fig_name + '.png')   
         
    
    return  
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
 
