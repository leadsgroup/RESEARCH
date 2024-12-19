from dash import html, dcc  
import dash_bootstrap_components as dbc 
from dash import Dash, html, dcc, Input, Output, callback
import plotly.express as px
import plotly.graph_objects as go
import json
import numpy as np
import pandas as pd

def main():
 
    routes_filename            = 'Data/American_Airlines_Monthly_Temp.xlsx'
    corn_feedstock_filename    = 'Data/corn_2022_production.csv'
    soybean_feedstock_filename = 'Data/soybean_2022_production.csv' 
    SAF_filename               = 'Data/Technology_Data.xlsx'
    
    Commercial_SAF                  = pd.read_excel(SAF_filename,sheet_name=['Commercial_SAF']) 
    Routes                          = pd.read_csv(routes_filename,sheet_name=['Sheet1']) 
    corn_feedstock                  = pd.read_csv(corn_feedstock_filename,sheet_name=['corn_2022_production']) 
    soybean_feedstock                = pd.read_excel(soybean_feedstock_filename,sheet_name=['soybean_2022_production']) 
    selected_fuels                 = ["Neste Jet-A","Neste SAF HEFA form Corn","World Energy SAF HEFA form Barley"]
    blend_ratios                    = [100,50,25]
    percent_fuel_use                = [75,95] 
    flight_range                    = [40000]
    SAF_production_states_1         = ["California"]  
    SAF_production_states_2         = ["Illinois"]  
    aircraft                        = ["Boeing 787-8"]
    selected_airpots                = ["Top 10 Airports","Top 20 Airports", "Top 50 Airports",  "All Airports"]
    percent_adoption                = [50]
    selected_feedstock              = ["Corn"]
    switch_off                      = False 
    
    generate_saf_production_map(Routes,Commercial_SAF,corn_feedstock,soybean_feedstock,selected_fuels,blend_ratios,percent_fuel_use,SAF_production_states_1,SAF_production_states_2,aircraft,flight_range,selected_airpots,percent_adoption,selected_feedstock,switch_off)
    #generate_saf_flight_ops_map(Routes,Commercial_SAF,corn_feedstock,soybean_feedstock,selected_fuels,blend_ratios,percent_fuel_use,SAF_production_states_1,SAF_production_states_2,aircraft,selected_airpots,percent_adoption,selected_feedstock,switch_off):  
    
    return 
 
def generate_saf_production_map(Routes,Commercial_SAF,corn_feedstock,soybean_feedstock,selected_fuels,blend_ratios,percent_fuel_use,SAF_production_states_1,SAF_production_states_2,aircraft,flight_range,selected_airpots,percent_adoption,selected_feedstock,switch_off):  
    template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]    
    font_size            = 16  
    
    with urlopen('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json') as response:
        counties = json.load(response)  
        
         
    # Step 1: compute the percentages of diffent types of fuels used 
    fuel_percentages_list = [0]
    fuel_percentages_list += percent_fuel_use
    fuel_percentages_list += [100]
    fuels_percentages     = np.diff(np.array(fuel_percentages_list))/100  
    
    
    # Step 2: determine the percentage of neat (pure) saf and Jet-A1 using blending ratios   
    num_fuels = len(selected_fuels)
    cumulative_fuel_use = np.zeros(num_fuels)  
    for i in range(1,num_fuels):
        # loop through fuels and get percentage of fuel used by each type
        cumulative_fuel_use[i] = fuels_percentages[i]* blend_ratios[i]
    cumulative_fuel_use[0]  = 1 - np.sum(cumulative_fuel_use[1:])
    
    # Step 3: Based on aircraft's fuel consumption, determine the total volume of fuel used oer flight 
    if aircraft == "Boeing 787-8": 
        C_T             = 0.000017678 
        max_fuel_volume = 126206  # liters 
        fuel_consumtion = 2.68 # L/100km per Pax 
        passengers      = 238
        W_0             = 227900   
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 377
    if aircraft == "Boeing 737 MAX-8":
        C_T             = 0.000017678  
        max_fuel_volume = 26024 # liters 
        fuel_consumtion = 2.04 # L/100km per Pax 
        passengers      = 162
        W_0             = 79015.8  
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 127 
    elif aircraft == 'Airbus A350-1000': 
        C_T             = 0.000017678   
        max_fuel_volume = 166488 
        fuel_consumtion = 2.85 # L/100km per Pax 
        passengers      = 327  
        W_0             = 322050  
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 464.3
    elif aircraft == 'Airbus A320 neo': 
        C_T             = 0.000017678  
        max_fuel_volume = 27200  # liters 
        fuel_consumtion = 1.94   # L/100km per Pax 
        passengers      = 180
        W_0             = 79015
        CL_cruise       = 0.055 
        CD_cruise       = 0.035 
        S_ref           = 122.6   
    
    Routes       =  Routes[ Routes['Range'] < flight_range[0]] # filter list of routes by range 
    miles_to_km  = 1.60934 
    flight_range = Routes['Range']*miles_to_km 
    passengers   = Routes['Passengers']
    fuel_volumes = fuel_consumtion*(flight_range/100)*passengers # Compute fuel volume used for each flight 
    Routes['Fuel_Volume'] = fuel_volumes
        
    # Step 4: Filter flight data based on option selected: i.e. top 10, top 20, top 50, all airpots 
    Airport_Routes     = Routes[['Passengers','Origin Airport','Destination City']]
    Cumulative_Flights = Airport_Routes.groupby(['Origin Airport']).sum()    
    if  selected_airpots == "Top 10 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(10) 
    elif  selected_airpots =="Top 20 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(20) 
    elif  selected_airpots == "Top 50 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(50) 
    elif  selected_airpots == "All Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False) 
    Airport_List = Busiest_Airports['Origin Airport']
    
    # Step 5: Filter airports that will support SAF and those that wont support SAF 
    mask                  = Routes['Origin Airport'].isin(Airport_List)
    SAF_Possible_Airports = Routes[mask]
    Non_SAF_Airports      = Routes[~mask]  
    
    # Step 6: Out of SAF supporting airports, use the percent adoption to detemrine how many flights at that airport will use SAF
    mask             = np.len(0,len(SAF_Possible_Airports), int(len(SAF_Possible_Airports) *percent_adoption ))
    SAF_Flights      = SAF_Possible_Airports[mask]  # flights on SAF 
    SAF_Flights_N    = SAF_Possible_Airports[~mask] 
    Non_SAF_Flights  = pd.concat([Non_SAF_Airports, SAF_Flights_N] )   # add list of flights from non supporting airports to non-SAF flights  
     
    # Step 7: Get total volume of eacb SAF required at the airports 
    total_fuel_volume_required = np.sum(np.float(SAF_Flights['Fuel_Volume'][1:]))
    unique_fuel_volumes        = cumulative_fuel_use*total_fuel_volume_required  
    
    # Step 8: Sort SAF's by feedstock and sum all fuel volumes based on fuel (FUTURE WORK: determine how much corn is needed for EACH type of process)
    mask                       = Commercial_SAF['Feedstock'].isin(selected_fuels)
    fuels_used                 = Commercial_SAF[mask]    
    fuels_used['Fuels Used']   = unique_fuel_volumes 
    relative_crop_volume       = fuels_used['Fuels Used']*fuels_used['Fuel-to-Crop-factor']
    fuels_used['Crop-Volume']  = relative_crop_volume
    df2                        = fuels_used.groupby(['Feedstock']).sum()  
    required_crop_vol          = df2[selected_feedstock]['Crop-Volume'] 
    
    # Step 9: Determine how many states will source the feedstock under consideration 
    if selected_feedstock == "Corn": 
        crop_data = corn_feedstock
    elif selected_feedstock == "Barley": 
        crop_data = soybean_feedstock 
    fips_data  = "{:02d}".format(crop_data['State ANSI']) +  "{:03d}".format(int(crop_data['County ANSI']))
    crop_data['FIPS'] = fips_data 
    
    # Step 10: Filter out states have not been selected and get number of states  
    crop_data_states      = SAF_production_states_1 + SAF_production_states_2 
    mask                  = crop_data['State'].isin(crop_data_states)
    feedstock_states      = crop_data[mask]    
    non_feedstock_states  = crop_data[~mask]     

    # Step 11: Randomize tracts in terms of crop area/usage then Recursively add rows until requied volume is met 
    Used_Feedstock        = feedstock_states.sample(frac = 1)
    tract_use_flag        = np.zeros(len(feedstock_states))
    Used_Feedstock["Feedstock_Usage"] = tract_use_flag
      
    add_row   = True
    i         = 1
    total_vol = 0
    while add_row == True:
        total_vol += np.float(Used_Feedstock.iloc[i]['fuel vols'])
        Used_Feedstock["Feedstock_Usage"][i] = 1
        i += 1
        if total_vol>required_crop_vol:
            add_row = False  
     
    # Step 11: Plot Land Use 
    Land_Use  = pd.concat([Used_Feedstock, non_feedstock_states] )   
    fips      = list(Land_Use['FIPS'])  
    Land_Use['FIPS'] = ["%05d" % i for i in fips] 
    fig= px.choropleth(Land_Use, geojson=counties, locations='FIPS', color = Feedstock_Usage,
                           color_continuous_scale="RdYlBu_r", 
                           hover_data=["Name","State", Feedstock_Usage],
                           scope='usa',
                           range_color=(0, 90),  
                          )   
    fig.update_layout(coloraxis_colorbar=dict(title=" "),
                                     coloraxis_colorbar_x=0.85, 
                                     height    = 400, 
                                     margin={'t':0,'l':0,'b':0,'r':0},                              
                              )  
    fig.update_coloraxes( colorbar_tickvals= np.linspace(0,90,11),
                                        colorbar_tickfont_size=font_size)  
    fig.show()    
  
    return  

def generate_saf_flight_ops_map(Routes,Commercial_SAF,corn_feedstock,soybean_feedstock,selected_fuels,blend_ratios,percent_fuel_use,SAF_production_states_1,SAF_production_states_2,aircraft,selected_airpots,percent_adoption,selected_feedstock,switch_off):  
    mapbox_access_token  = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ"     
    map_style            = None if switch_off else 'dark'  
    template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]    
    
    
    if aircraft == "Boeing 787-8": 
        C_T             = 0.000017678 
        max_fuel_volume = 126206  # liters 
        fuel_consumtion = 2.68 # L/100km per Pax 
        passengers      = 238
        W_0             = 227900   
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 377
    if aircraft == "Boeing 737 MAX-8":
        C_T             = 0.000017678  
        max_fuel_volume = 26024 # liters 
        fuel_consumtion = 2.04 # L/100km per Pax 
        passengers      = 162
        W_0             = 79015.8  
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 127 
    elif aircraft == 'Airbus A350-1000': 
        C_T             = 0.000017678   
        max_fuel_volume = 166488 
        fuel_consumtion = 2.85 # L/100km per Pax 
        passengers      = 327  
        W_0             = 322050  
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 464.3
    elif aircraft == 'Airbus A320 neo': 
        C_T             = 0.000017678  
        max_fuel_volume = 27200  # liters 
        fuel_consumtion = 1.94   # L/100km per Pax 
        passengers      = 180
        W_0             = 79015
        CL_cruise       = 0.055 
        CD_cruise       = 0.035 
        S_ref           = 122.6 
             
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # Compute distance between departure and destimation points
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    Routes_Mo                    = Routes[Routes['Month'] == month_no+1 ]  
    des_lon                      = np.array(Routes_Mo['Destination Longitude (Deg.)'])
    des_lat                      = np.array(Routes_Mo['Destination Latitude (Deg.)'])
    org_lon                      = np.array(Routes_Mo['Origin Longitude (Deg.)'])
    org_lat                      = np.array(Routes_Mo['Origin Latitude (Deg.)']) 
    origin_coordinates           = np.stack((des_lat,des_lon))
    destination_coordinates      = np.stack((org_lat, org_lon)) 
    R                            = 6371.0088*1000 
    coord0_rad                   = origin_coordinates*0.017453292519943295
    coord1_rad                   = destination_coordinates*0.017453292519943295
    angles                       = np.arccos(np.sin(coord0_rad[0,:])*np.sin(coord1_rad[0,:]) + 
                                            np.cos(coord0_rad[0,:])*np.cos(coord1_rad[0,:])*np.cos(coord0_rad[1,:] - coord1_rad[1,:]))
    distances                    = R*angles  
    Routes_Mo['Range']  = distances.tolist() 

    
    fuel_reserve        = 1800 # 30 mins 
    reserve_speed       = 250 # m/s
    reserve_distance    = reserve_speed*fuel_reserve
    air_density         = 8.21392E-3
    liters_to_kg        = 1/1.250
    W_f                 = liters_to_kg*fuel_consumtion*passengers*(distances+ reserve_distance)/(1000*100)
    W_1                 = W_0*percent_fuel_use - W_f
    Range               = (1/C_T)*(np.sqrt(CL_cruise)/CD_cruise)*np.sqrt(2/(air_density*S_ref))*(np.sqrt(W_0) -np.sqrt(W_1) ) 
    Infeasible_Routes   = Routes_Mo[Routes_Mo['Range'] > Range ]  
    Feasible_Routes     = Routes_Mo[Routes_Mo['Range'] < Range ]  
     
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # Routes 
    # ---------------------------------------------------------------------------------------------------------------------------------------------- 
    mapbox_access_token  = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ"     
    fig                  = go.Figure()
    airport_marker_size  = 5
    airport_marker_color = "white"

    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # Flight Paths
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    lons       = np.empty(3 * len(Infeasible_Routes))
    lons[::3]  = Infeasible_Routes['Origin Longitude (Deg.)']
    lons[1::3] = Infeasible_Routes['Destination Longitude (Deg.)']
    lons[2::3] = None
    lats       = np.empty(3 * len(Infeasible_Routes))
    lats[::3]  = Infeasible_Routes['Origin Latitude (Deg.)']
    lats[1::3] = Infeasible_Routes['Destination Latitude (Deg.)']
    lats[2::3] = None    
  
    fig.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            line = dict(width = 0.1,color = 'grey'),  
        )
    )
    
    lons       = np.empty(3 * len(Feasible_Routes))
    lons[::3]  = Feasible_Routes['Origin Longitude (Deg.)']
    lons[1::3] = Feasible_Routes['Destination Longitude (Deg.)']
    lons[2::3] = None
    lats       = np.empty(3 * len(Feasible_Routes))
    lats[::3]  = Feasible_Routes['Origin Latitude (Deg.)']
    lats[1::3] = Feasible_Routes['Destination Latitude (Deg.)']
    lats[2::3] = None    
  
    fig.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            line = dict(width = 2,color = "aquamarine"), 
        )
    )
 

    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # Airports 
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    fig.add_trace(go.Scattergeo( 
        lon = Routes['Destination Longitude (Deg.)'],
        lat = Routes['Destination Latitude (Deg.)'], 
        text = Routes['Destination City'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color, 
        ))) 

    fig.add_trace(go.Scattergeo( 
        lon = Routes['Origin Longitude (Deg.)'],
        lat = Routes['Origin Latitude (Deg.)'], 
        text = Routes['Origin City'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color, 
        ))) 
    
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # Flight Paths
    # ---------------------------------------------------------------------------------------------------------------------------------------------- 
    fig.update_layout(mapbox_style  = "open-street-map",      
                      showlegend    = False, 
                      height        = 400, 
                      geo_scope     ='usa',
                      margin        = {'t':0,'l':0,'b':0,'r':0},  
                      mapbox        = dict( accesstoken=mapbox_access_token,style=map_style,
                                            center=go.layout.mapbox.Center( lat=30, lon= 230 ))  )   
     

    fig["layout"]["template"] = template
    
    

    fig.show()    
        
    return fig


    
if __name__ == '__main__':
    main()