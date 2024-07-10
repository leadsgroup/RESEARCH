'''
# Simulation_Repeated_Flight_Operations.py
#
# Created: May 2024, S S. Shekar

'''

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units, Data   
from RCAIDE.Library.Plots                                           import *  

import time  
import numpy as np
import pylab as plt
import matplotlib.pyplot as plt

import folium
from selenium import webdriver



def main():

    create_map_png()
    convert_html()
   
    
    return



def create_map_png():
    
    # Coordinates for each location
    locations = {
        "Memphis": [35.0443, -89.9766],
        "Louisville": [38.1707, -85.7308],
        "Nashville": [36.1245, -86.6781],
        "St. Louis": [38.7499, -90.3748],
        "Cincinnati": [39.0514, -84.6671],
        "Columbus": [39.9999, -82.8872],
        "Champaign": [40.0365, -88.2640]
    }
    
    # Create a map centered at the average location of all coordinates
    average_lat = sum([coord[0] for coord in locations.values()]) / len(locations)
    average_lon = sum([coord[1] for coord in locations.values()]) / len(locations)
    map_plot = folium.Map(location=[average_lat, average_lon], zoom_start=6)
    
    # Add markers for each location
    for name, coordinates in locations.items():
        folium.Marker(
            location=coordinates,
            popup=name,
        ).add_to(map_plot)
   
   
    
    # Create a map centered at O'Hare International Airport
    #map_plot = folium.Map(location=ord_coordinates, zoom_start=8)
    
    # Add a marker for O'Hare International Airport
    #folium.Marker(ord_coordinates, popup="MEM").add_to(map_plot)
    folium.TileLayer('openstreetmap').add_to(map_plot)
    
    #folium.TileLayer('cartodb positron').add_to(chicago_map)
    #folium.TileLayer('cartodb dark_matter').add_to(chicago_map)
    #folium.TileLayer('Stamen Toner').add_to(chicago_map)
  
    
    
    #tiles='Stamen Terrain'
    
    # Radius in miles (e.g., 10 miles)
    radius_miles_1 = 110
       
    
    # Convert miles to meters (1 mile = 1609.34 meters)
    radius_meters_1 = radius_miles_1 * 1609.34

    
    #for name, coordinates in locations.items():
        #folium.Circle(
              #location=coordinates,
              #radius=radius_meters_1,
              #color='blue',
              #fill=True,
              #fill_color='blue',
              #fill_opacity=0.2,
          #).add_to(map_plot)
       
    
    # Save the map as an HTML file
    map_plot.save("folium_maps/midwest_map.html")
    
    map_plot
    
    return

def convert_html():
    # Path to the HTML file
    html_path = "file:///Users/sai/Documents/Research/RESEARCH/14_BTMS_Liqiud_Cooling/AIAA_Aviation_2024/Map_Plots/folium_maps/midwest_map.html"
    
    ## Set up the Selenium WebDriver (ensure you have the correct path to your webdriver)
    #options = webdriver.ChromeOptions()
    #options.add_argument("--headless")
    #options.add_argument("--disable-gpu")
    #options.add_argument("--no-sandbox")
    driver = webdriver.Chrome()
    
    # Open the HTML file in the browser
    driver.get(html_path)
    
    # Set the window size to the desired dimensions
    driver.set_window_size(1200, 1200)
    
    # Give it a few seconds to load
    time.sleep(2)
    
    # Save the screenshot
    screenshot_path = "folium_maps/midwest_map.png"
    driver.save_screenshot(screenshot_path)
    
    # Close the browser
    driver.quit()
    
      
    return


if __name__ == '__main__': 
    main()    
    plt.show()
