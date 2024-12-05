#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
from RCAIDE.Framework.Core import Units, Data   
from RCAIDE.Library.Plots import *       
from RCAIDE import  load 
from RCAIDE import  save
from matplotlib import  pyplot as  plt
import  pickle

# python imports 
import os 
import pickle
import sys 
import pandas as pd
import numpy as  np 

local_path_1 =  os.path.split(sys.path[0])[0] 
local_path_2 =  os.path.split(os.path.split(sys.path[0])[0])[0]  

sys.path.append( os.path.join(local_path_1, 'Flight_Path_Functions'))
sys.path.append( os.path.join(local_path_1, 'Post_Processing_Functions'))
from Aircraft_Noise_Emissions         import post_process_approach_departure_noise_data
# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def main():            
     
    # ----------------------------------------------------------------------------------------------------------------------
    #  SIMULATION SETTINGS 
    # ----------------------------------------------------------------------------------------------------------------------  
    aircraft_codes                   = ['HC','TSR','TR']  

    number_of_microphone_in_stencil  = 90000
    noise_evaluation_pitch           = 100 * Units.feet  
    aircraft_origin_location         = np.array([0, 0])
    N_gm_x                           = 300 # mics every 50 feet
    N_gm_y                           = 300 # mics every 50 feet
    Mic_min_x                        = -7500 *  Units.feet
    Mic_min_y                        = -7500*  Units.feet
    Mic_max_x                        =  7500*  Units.feet
    Mic_max_y                        =  7500*  Units.feet
    
    
    filename_list = ['HC_Hemisphere_Analyzed_Noise', 'TR_Hemisphere_Analyzed_Noise', 'TSR_Hemisphere_Analyzed_Noise']
    processed_filenames =  []
    for i,  analyzed_filename in enumerate(filename_list):        
        analyzed_results   = load(analyzed_filename + '.res')
        processed_results  = post_process_approach_departure_noise_data(analyzed_results, number_of_microphone_in_stencil, noise_evaluation_pitch, aircraft_origin_location,  N_gm_x, N_gm_y, Mic_min_x, Mic_min_y, Mic_max_x, Mic_max_y)
        processed_filename = aircraft_codes[i]  +  '_Hemisphere_Processed_Noise'
        processed_filenames.append(processed_filename)
        save(processed_results, processed_filename + '.res') 
        
        F =  Data(filename_list=processed_filenames)
        postprocessed_filename_name = '_Hemisphere_Postprocessed_Noise' 
        save(F, postprocessed_filename_name + '.res')                          
              
    return

 
if __name__ == '__main__': 
    main()     