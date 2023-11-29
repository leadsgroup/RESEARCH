# MARC Imports 
from  MARC.Core import Data, Units  

# Python Imports  
import numpy as np   
import pickle 
import matplotlib.pyplot as plt  
import matplotlib.cm as cm 
from MARC.Visualization.Performance.Common.post_process_noise_data import post_process_noise_data
 
# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():  

    save_figure = True 
    filetype    = '.png'
    
    # Universal Plot Settings 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 28,
                  'xtick.labelsize': 24,
                  'ytick.labelsize': 24,
                  'axes.titlesize': 28}
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
    
  
    aircraft_tags  = ['TR','HC','SR'] 
    angles         = [0,45,90,135,180]      
    num_aircrafts  = len(aircraft_tags) 
    
    
    for a_i in range(num_aircrafts):
        
        for i in range(len(angles)):  
            
            # initialize plots 
            filename     = 'Raw_Data_' + aircraft_tags[a_i]  + '/'  + aircraft_tags[a_i]  + '_TG_Noise_Angle_' + str(angles[i]) 
            results      = load_results(filename)  

            fig = plt.figure(filename)    
            fig.set_size_inches(plot_parameters.figure_width,plot_parameters.figure_height)  
            axis = fig.add_subplot(1,1,1)       
             
            # format results 
            noise_data   = post_process_noise_data(results)  
            SPL          = noise_data.SPL_dBA    
            max_SPL      = np.flip(np.max(SPL,axis=0),axis = 1)      
               
        
            Range_x    = noise_data.receptor_locations[:,0,0]/Units.nmi + 0.5
            Range_y    = noise_data.receptor_locations[0,:,1]/Units.nmi + 1.0
            Y, X       = np.meshgrid(Range_y, Range_x)     
            

            #  quantify nodes in bands at airports   
            aircraft_band_45_65   = percent_exposure(max_SPL, 45, 65)     
            aircraft_band_65_100   = percent_exposure(max_SPL, 65, 100)   
            
            print(filename)
            print('% 45_65 ', aircraft_band_45_65   ) 
            print('% 65_90 ', aircraft_band_65_100   )      
            
                        
            levs                = np.linspace(45,90,10)  # np.array([45,50,55,60,70,80,90])  
            CS                  = axis.contourf(X , Y, max_SPL, levels  = levs, cmap=plt.cm.jet, extend='both')   
            cbar                = fig.colorbar(CS)
            cbar.ax.set_ylabel('LA$_{max}$ [dBA]', rotation =  90)      
            axis.set_xlabel('Longitudinal Distance [nmi]')
            axis.set_ylabel('Latitudinal Distance [nmi]',labelpad = 15) 
            
            axis.grid(False)  
            axis.minorticks_on()    
            
            plt.tight_layout()
            # save results 
            if save_figure:
                safe_filename = '../Papers_and_Presentation/Images/' + aircraft_tags[a_i]  + '_TG_Noise_Angle_' + str(angles[i]) 
                plt.savefig(safe_filename + filetype)   
         
    
    return  

def percent_exposure(X, lower_bound, upper_bound):
    count = 0
    x = X.flatten()
    for i in range(len(x)):
        if x[i] >= lower_bound and x[i] < upper_bound: 
            count += 1  
    return 100*count/len(x)  
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
 
