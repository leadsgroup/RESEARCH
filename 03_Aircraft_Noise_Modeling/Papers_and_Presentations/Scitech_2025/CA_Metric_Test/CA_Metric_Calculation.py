import pandas as pd
import os 
import  numpy as np 

def main():

    ospath          = os.path.abspath(__file__)
    separator       = os.path.sep
    relative_path   = os.path.dirname(ospath) + separator 
    routes_filepath = relative_path  + '_Raw_Data_Tract_Revised'
    
    gdf_census                  = pd.read_csv( routes_filepath + '.csv')  
    sensitivity_levels          = sensitivity_levels = {"Churches": .5,"Schools": .83,"Hospitals/Medical Centers": .5}
    demographic_col_headers     = ['Asian', 'Black', 'Native American' , 'Pacific Islander' , 'Hispanic', 'White'] 
    new_demographic_headers     = ['Asian', 'Black', 'Native_American' , 'Pacific_Islander' , 'Hispanic', 'White'] 
    income_columns              = ['<50K','50k to 100k','100k to 150k','150k to 200k','200k+'] 
    noise_sensitive_structures  = ['Churches','Hospitals/Medical Centers','Schools']
    noise_columns               = ['L_dn']  
    
    for i in range(len(demographic_col_headers)):
        tract_race =  gdf_census[demographic_col_headers[i]]
        PR         =  tract_race / np.sum(tract_race)
        
        Noise      =  gdf_census[noise_columns[0]]
         
         
        summation_N_j  = np.array(gdf_census[noise_sensitive_structures[0]] + gdf_census[noise_sensitive_structures[1]] +   gdf_census[noise_sensitive_structures[2]])
        S_net = 1 + ( np.array(gdf_census[noise_sensitive_structures[0]] * sensitivity_levels[noise_sensitive_structures[0]])+ \
                       np.array(gdf_census[noise_sensitive_structures[1]] * sensitivity_levels[noise_sensitive_structures[1]]) + \
                       np.array(gdf_census[noise_sensitive_structures[2]] * sensitivity_levels[noise_sensitive_structures[2]]) ) /summation_N_j
          
        S_net[np.isnan(S_net)] = 1 
         
        demograhpic_CA =  Noise * PR *  S_net *  (np.log(summation_N_j+1) +  1)
         

        gdf_census[f'{new_demographic_headers[i]}_L_nd_PR'] =  Noise * PR     
        gdf_census[f'{new_demographic_headers[i]}_CA']      =  demograhpic_CA
        
    
    file_name =  'HC_Raw_Data_Tract_Updated' 
    gdf_census.to_csv( file_name+ '.csv', index=False)
    return  


# Run the main function
if __name__ == '__main__':
    main()
    plt.show()