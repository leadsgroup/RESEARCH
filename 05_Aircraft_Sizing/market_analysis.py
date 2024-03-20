# imports  
import numpy as np           
import pandas as pd  
import pandas as pd
import matplotlib.pyplot as plt
import os


def main():  
    
    # Read the data from the Excel file
    ospath = os.path.abspath(__file__)
    separator = os.path.sep
    rel_path = os.path.dirname(ospath)
    file_path = rel_path + separator + 'T_T100D_MARKET_ALL_CARRIER_data.xlsx'
    df = pd.read_excel(file_path)
    
    # Extract passenger and distance data
    no_passengers = df['PASSENGERS']
    distance = df['DISTANCE']
    
    # Create a dictionary to store the sum of passengers for each distance
    passenger_sum_per_distance = {}
    
    # Sum up the number of passengers for each distance
    for i in range(len(distance)):
        dist = distance[i]
        passengers = no_passengers[i]
        if dist in passenger_sum_per_distance:
            passenger_sum_per_distance[dist] += passengers
        else:
            passenger_sum_per_distance[dist] = passengers
    
    # Plot histogram with customizations
    plt.figure(figsize=(12, 8))
    plt.hist(
        passenger_sum_per_distance.keys(),
        weights=passenger_sum_per_distance.values(),
        bins=range(0, max(passenger_sum_per_distance.keys()) + 100, 100),
        edgecolor='black',
        color='skyblue',
        alpha=0.7  # Transparency
    )
    plt.xlabel('Distance')
    plt.ylabel('Number of Passengers')
    plt.title('Histogram of Number of Passengers for Each Distance')
    plt.grid(True, linestyle='--', alpha=0.5)  # Add grid lines
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    plt.tight_layout()  # Adjust layout to prevent clipping of labels
    plt.show()
    
    
    return 
     
if __name__ == '__main__':
    main()