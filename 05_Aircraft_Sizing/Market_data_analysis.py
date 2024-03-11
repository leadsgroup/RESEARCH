import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
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
    
    # annual US economy 
    font_size           = 12
    fig               = go.Figure()
    sector_colors       = px.colors.qualitative.Pastel 
    fig.add_trace(go.Histogram(histfunc="sum",
                               x=  df['DISTANCE'],
                               y =  df['PASSENGERS'], 
                               xbins=dict(start=0, end=4000, size=500),
                               marker_color=sector_colors[0],
                               )
                  ) 
    
    # The two histograms are drawn on top of another
    fig.update_layout(barmode='stack', 
                      xaxis_title_text='Distance (miles)', # xaxis label
                      yaxis_title_text='Passengers', # yaxis label 
                      height        = 300, 
                      width         = 600, 
                      margin        = {'t':0,'l':0,'b':0,'r':0},  
                      bargap        = 0.1,
                      font=dict(  size=font_size ))     
    
    
    plt.figure(figsize=(12, 8))
    plt.hist(
        passenger_sum_per_distance.keys(),
        weights=passenger_sum_per_distance.values(),
        bins=range(0, max(passenger_sum_per_distance.keys()) + 100, 100),
        edgecolor='black',
        color='skyblue',
        alpha=0.7  # Transparency
    )
    
    fig.show()
     
    # estimated daily operations  
     
    
    #plt.xlabel('Distance')
    #plt.ylabel('Number of Passengers')
    #plt.title('Histogram of Number of Passengers for Each Distance')
    #plt.grid(True, linestyle='--', alpha=0.5)  # Add grid lines
    #plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
    #plt.tight_layout()  # Adjust layout to prevent clipping of labels
    
if __name__ == '__main__':
    main() 