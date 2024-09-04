import numpy as np
import matplotlib.pyplot as plt  
import pandas as pd
import pickle


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():  
     
    # tag for city  
    city                        = 'LA' 
    aircraft_models             = ['HC','SR','TR'] 
    altitudes                   = ['1000']
    flight_intervals            = ['60','30','10'] 
      
    for m_i in range(len(aircraft_models)): 
        for alt_i in range(len(altitudes)):  
            for f_i in range(len(flight_intervals)):
                
                aircraft  = aircraft_models[m_i] 
                altitude  = altitudes[alt_i]
                frequency = flight_intervals[f_i]
            
                # load data  
                processed_results_filename = 'Raw_Data_' + aircraft + '/'+ aircraft + '_' + altitude + 'ft_' + city + '_' + frequency +  'min_All'   
                Results = load_results(processed_results_filename) 
                
                
                # create excel files of data below 
                elevation      = (Results.elevation).flatten() # this is in 2D 
                LAT            = (Results.lat_deg).flatten()
                LONG           = (Results.long_deg).flatten()    
                L_Amax         = (Results.L_Amax).flatten() 
                L_AeqT         = (Results.L_AeqT).flatten()  
                L_AeqT_24hr    = (Results.L_AeqT_24hr).flatten()
                SEL            = (Results.SEL).flatten()
                L_dn           = (Results.L_dn).flatten()
                L_Aeq_jetliner = (Results.L_AeqT_24hr_total).flatten()   
                
                M                   = {}
                M['elevation']      = elevation   
                M['Latitude']       = LAT           
                M['Longitude']      = LONG          
                M['L_Amax']         = L_Amax        
                M['L_AeqT']         = L_AeqT        
                M['L_AeqT_24hr']    = L_AeqT_24hr   
                M['SEL']            = SEL           
                M['L_dn']           = L_dn          
                M['L_Aeq_jetliner'] = L_Aeq_jetliner 
                
                Mat =  pd.DataFrame(M)
                
                csv_filename = aircraft + '_' + altitude + 'ft_' + city + '_' + frequency +  'min_All.csv' 
                Mat.to_csv(csv_filename)
                

                    
                    
                    
                    
                    
                    
                # save data in excel file 
                
                

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
