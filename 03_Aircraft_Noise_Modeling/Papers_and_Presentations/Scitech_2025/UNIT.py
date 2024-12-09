from DataPreProcessing import*
from NoiseImpactMapping import*
from CommunityImpactMapping import*

import geopandas as gpd
import pandas as pd
import os

import time 

race_dictionary = {
    "Total": "B03002001",
    "White": "B03002003",
    "Black or African American": "B03002004",
    "Native American": "B03002005",
    "Asian": "B03002006",
    "Native Hawaiian and Other Pacific Islander": "B03002007",
    "Hispanic or Latino": "B03002012"
}

income_dictionary = {
    "Total": "B19001001",
    "Less than $10,000": "B19001002",
    "$10,000 to $14,999": "B19001003",
    "$15,000 to $19,999": "B19001004",
    "$20,000 to $24,999": "B19001005",
    "$25,000 to $29,999": "B19001006",
    "$30,000 to $34,999": "B19001007",
    "$35,000 to $39,999": "B19001008",
    "$40,000 to $44,999": "B19001009",
    "$45,000 to $49,999": "B19001010",
    "$50,000 to $59,999": "B19001011",
    "$60,000 to $74,999": "B19001012",
    "$75,000 to $99,999": "B19001013",
    "$100,000 to $124,999": "B19001014",
    "$125,000 to $149,999": "B19001015",
    "$150,000 to $199,999": "B19001016",
    "$200,000 or More": "B19001017"
}

ti=time.time()

# Define the separator for file paths
separator = os.path.sep

# Update the filenames for raw census data (race and income) in the Processed_Data directory
base_geoJson = 'LA_Area_Tract_TR.geojson'
noise_filename = 'Noise_Data' + separator + 'Cumulative_TR_LA_1000ft.csv'

# base_geoJson = 'LA_Area_Tract_HC.geojson'
# noise_filename = 'Noise_Data' + separator + 'Cumulative_HC_LA_1000ft.csv'

struct_filename1 = 'Raw_Data' + separator + 'Churches.csv'
struct_filename2 = 'Raw_Data' + separator + 'Schools_Colleges_and_Universities.csv'
struct_filename3 = 'Raw_Data' + separator + 'Hospitals_and_Medical_Centers.csv'


#load in files
geoJson_file = gpd.read_file(base_geoJson)
noise_file = pd.read_csv(noise_filename)
struct_file1 = pd.read_csv(struct_filename1)
struct_file2 = pd.read_csv(struct_filename2)
struct_file3 = pd.read_csv(struct_filename3)


struct_list=[struct_file1,struct_file2,struct_file3]


#-------------------UNIT---------------#

#Define Variables Based on Analysis
sensitivity_levels = {
"Churches": .5,
"Schools": .83,
"Hospitals/Medical Centers": .5
}

#Array of values you wish to analyze for noise or census
census_columns=[race_dictionary['White'],race_dictionary['Black or African American'],race_dictionary['Hispanic or Latino'],
               race_dictionary['Asian'], race_dictionary['Native Hawaiian and Other Pacific Islander']]
noise_columns=['L_dn']

#-------------Preprocessing Data---------#

#Merging Data into one geoJSON
noise__struct_file = sensitive_area_preprocess(struct_list,geoJson_file,list(sensitivity_levels.keys()))
# save_file(noise__struct_file,'LA_Area_Tract.geojson')

# noise_census_file = noise_preprocess(noise_file,geoJson_file)
# convert_csv(noise_census_file,'TR_Data_Tract.csv')

# save_file(noise_census_file,'LA_Area_Tract_TR.geojson')


c_a_gdf = community_annoyance(geoJson_file,list(sensitivity_levels.keys()),sensitivity_levels)


#Save merged data file (optional)
# save_file(noise_census_file,'LA_Area_Tract_HC.geojson')

# convert_csv(noise_census_file,'HC_Data_Tract.csv')
# convert_csv(c_a_gdf,'HC_Raw_Data_Tract.csv')
convert_csv(c_a_gdf,'TR_Raw_Data_Tract.csv')

# convert_csv(c_a_gdf,'TR_Data_Tract_LogPOP_Only.csv')



#-----------------Plotting--------------#

#Plotting function
# impact_plotting(noise_census_file,noise_columns,min_lon= -119,min_lat=33,max_lon= -117.25,max_lat=35)
# census_plotting(noise_census_file,census_columns,min_lon= -119,min_lat=33,max_lon= -117.25,max_lat=35)



to = time.time()

print((to-ti)/60)