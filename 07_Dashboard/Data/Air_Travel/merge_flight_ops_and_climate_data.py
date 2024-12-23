import numpy as np           
import pandas as pd 

def main():
    '''This scripts appends the temperature data of the counties within the US to data frame with airports recursively '''
    
    # load data
    airline          = 'American_Airlines'
    airline_file     = 'American_Airlines_Flight_Operations_2019.csv'
    temperature_file = '../US_Climate/Monthly_US_County_Temperature_2019.xlsx'
    Airline_data     = pd.read_csv(airline_file)     # pd.read_excel(airline_file,sheet_name=[airline])[airline] 
    Temperature_data = pd.read_excel(temperature_file,sheet_name=['US_County_Temperature_F'])['US_County_Temperature_F']  
    
    # county lats and longs 
    county_lats   = np.array(Temperature_data['Latitude'])
    county_longs  = np.array(Temperature_data['Longitude'])
    month_headers = list(Temperature_data.columns.values)[6:18]

    # loop through latitude and longitude points 
    temperature = np.empty((24 , len(Airline_data)))

    for i in range(len(Airline_data)):
        # find index of nearest county 
        origin_lat        = Airline_data['Origin Latitude (Deg.)'][i]
        origin_long       = Airline_data['Origin Longitude (Deg.)'][i]
        destination_lat   = Airline_data['Destination Latitude (Deg.)'][i]
        destination_long  = Airline_data['Destination Longitude (Deg.)'][i]
        
        # minimim degree
        origin_min_deg              = ((county_lats-origin_lat)**2 + (county_longs-origin_long)**2)**(0.5)
        origin_nearest_county       = np.where(origin_min_deg == origin_min_deg.min())[0][0] 
        destination_min_deg         = ((county_lats-destination_lat)**2 + (county_longs-destination_long)**2)**(0.5)
        destination_nearest_county  = np.where(destination_min_deg == destination_min_deg.min())[0][0]
        
        # save data 
        temperature[:12,i] = np.array(Temperature_data[month_headers].loc[origin_nearest_county])
        temperature[12:,i] = np.array(Temperature_data[month_headers].loc[destination_nearest_county])
    
    Aggredated_data         = Airline_data[['Passengers', 'Distance (miles)', 'Carrier Code', 'Airline',
       'Origin Airport', 'Origin City Name',
       'Origin Latitude (Deg.)', 'Origin Longitude (Deg.)',
       'Destination Airport', 'Destination City Name',
       'Destination Latitude (Deg.)', 'Destination Longitude (Deg.)', 'Month',
       'Estimated Aircraft Capacity', 'SFC (L/100km per Pax)',
       'No of Flights Per Month', 'Climb-Descent Penalty (mi)',
       'Fuel Consumed Per Flight (Liters)', 'Total Fuel Per Route (Gal)',
       'Fuel Cost Per Gal', 'Fuel Cost']]
        
    # save new data file
    headers = ['Origin January', 'Origin February', 'Origin March', 'Origin April', 'Origin May', 'Origin June',
               'Origin July', 'Origin August', 'Origin September', 'Origin October', 'Origin November', 'Origin December' ,
               'Destination January', 'Destination February', 'Destination March', 'Destination April', 'Destination May', 'Destination June',
               'Destination July', 'Destination August', 'Destination September', 'Destination October', 'Destination November', 'Destination December']
    Aggredated_data[headers] = temperature.T
    Aggredated_data.to_csv( airline + "_Flight_Ops_and_Climate.csv") 
    
    return 
     
if __name__ == '__main__':
    main()