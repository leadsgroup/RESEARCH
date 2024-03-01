# imports 
import plotly.express as px
import plotly.graph_objects as go
from plotly.graph_objs import *
import numpy as np           
import pandas as pd 
from urllib.request import urlopen
import json 

def main(): 
 
    with urlopen('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json') as response:
        counties = json.load(response)    
        
    filename = '../Data/LEADS_SAT_Dashboard_Data.xlsx'
    SAT_data = pd.read_excel(filename,sheet_name=['Commercial_Batteries','Battery_Research', 'Air_Travel','US_Temperature_F']) 
     
    US_Temperature_F   = SAT_data['US_Temperature_F']   
     
     
    names = list(US_Temperature_F.columns.values)[4:18]
          
    fips = list(US_Temperature_F['FIPS'])  
    US_Temperature_F['FIPS'] = ["%05d" % i for i in fips] 
    us_temperature_map= px.choropleth(US_Temperature_F, geojson=counties, locations='FIPS', color = month,
                           color_continuous_scale="RdYlBu_r",
                           hover_data=["Name","State", month],
                           scope='usa',
                           range_color=(20, 90),  
                          )  
     
    us_temperature_map.update_layout(coloraxis_colorbar=dict(title=" "),
                                     coloraxis_colorbar_x=0.85, 
                                     height    = 400, 
                                     margin={'t':0,'l':0,'b':0,'r':0},                              
                              )  
    us_temperature_map.update_coloraxes(colorbar_orientation= "v",
                                        colorbar_tickvals= np.linspace(0,90,11),
                                        colorbar_tickfont_size=20)
    us_temperature_map.show()
    
    return 
     
if __name__ == '__main__':
    main()