import plotly.express as px
import plotly.graph_objects as go
from plotly.graph_objs import *
import numpy as np           
import pandas as pd
from scipy.io import netcdf 
from  mpl_toolkits.basemap import Basemap
from urllib.request import urlopen
import json
import plotly.io as pio
import plotly.figure_factory as ff 
import os
 

def main():
   
    separator                       = os.path.sep
    SAF_filename                    = 'Data' + separator + 'Technology' + separator +  'Technology_Data.xlsx'
    routes_filename                 = 'Data' + separator +  'Air_Travel' + separator + 'American_Airlines_Flight_Ops_and_Climate.xlsx'
    crops_filename                  = 'Data' + separator +  'Crops' + 'All_Crops_2017.csv' 
    
    SAT_data                       = pd.read_excel(SAF_filename,sheet_name=['Commercial_SAF','SAF_Research']) 
    Commercial_SAF                 = SAT_data['Commercial_SAF'] 
    SAF_Research                   = SAT_data['SAF_Research']  
    a                              = Commercial_SAF['Brand']  
    b                              = Commercial_SAF['Fuel Type']
    c                              = Commercial_SAF['Process']
    d                              = Commercial_SAF['Source']
    e                              = Commercial_SAF['Feedstock']
    Commercial_SAF["SAF Name"]     = a + ' ' + b + 'from' + c  + 'using ' + d + ' ' + e 
    Flight_Ops                     = pd.read_excel(routes_filename,sheet_name=['Sheet1']) 
    Flight_Ops                     = Flight_Ops['Sheet1']     
    feedstocks                     = pd.read_csv(crops_filename,sheet_name=['Corn','Soyabean','Barley','Canola'])  
  

    selected_process   = 'All'
    selected_feedstock = 'All'
    selected_x_axis    = list(Commercial_SAF.columns.values)[4:18][7] # NEED TO CORRECT 
    selected_y_axis    = list(Commercial_SAF.columns.values)[4:18][9] # NEED TO CORRECT 
    switch_off         = False     
    fig_1 = generate_saf_scatter_plot(Commercial_SAF,selected_process,selected_feedstock,selected_x_axis,selected_y_axis,switch_off)
    

    saf_1           = Commercial_SAF['SAF Name'][4] 
    saf_2           = Commercial_SAF['SAF Name'][1] 
    saf_3           = Commercial_SAF['SAF Name'][3] 
    fig_2 = generate_saf_spider_plot(Commercial_SAF,saf_1,saf_2,saf_3,switch_off)

    selected_sector = 'All'
    selected_type   = 'All'
    fig_3 = generate_saf_dev_map(SAF_Research,selected_sector,selected_type,switch_off)
    
    
    selected_fuels                  = [saf_1,saf_2,saf_3]
    blend_ratios                    = [100,50,25]
    percent_fuel_use                = [75,95] 
    flight_range                    = [40000] # given in miles 
    SAF_production_states_1         = ["California"]  
    SAF_production_states_2         = ["Illinois"]  
    aircraft                        = ["Boeing 787-8"]
    selected_airpots                = ["Top 10 Airports","Top 20 Airports", "Top 50 Airports",  "All Airports"]
    percent_adoption                = [50] 
    month_no                        = 1
    switch_off                      = False 
    
    fig_4 = generate_feedstock_usage_plots(Flight_Ops,Commercial_SAF,feedstocks,selected_fuels,blend_ratios,percent_fuel_use,SAF_production_states_1,SAF_production_states_2,aircraft,flight_range,selected_airpots,percent_adoption,selected_feedstock,month_no,switch_off)
      
    fig_5 = generate_saf_flight_operations_plots(Flight_Ops,Commercial_SAF,feedstocks,selected_fuels,blend_ratios,percent_fuel_use,SAF_production_states_1,SAF_production_states_2,aircraft,flight_range,selected_airpots,percent_adoption,selected_feedstock,month_no,switch_off)
    
    return 
     
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# SAF Plots 
# ---------------------------------------------------------------------------------------------------------------------------------------------------
def generate_saf_scatter_plot(Commercial_SAF,selected_process,selected_feedstock,selected_x_axis,selected_y_axis,switch_off): 
    template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]  
    unique_processes     = list(Commercial_SAF['Process'][1:].unique())
    unique_feedstocks    = list(Commercial_SAF['Feedstock'][1:].unique())
    marker_size          = 15
    opacity_ratio        = 0.8 if switch_off else 1.0
    font_size            = 16 

    # Process Colors: greenyellow, aquamarine, paleturquoise, lightcoral, yellow, lavender ,thistle ,orangered   
    Process_Colors      = px.colors.qualitative.Pastel 
    Feedstock_Markers  = ['square','x','circle','cross','diamond','triangle-up','triangle-down','star','hourglass'] 
    fig = go.Figure()
    if selected_process == 'All' and  selected_feedstock == 'All':
        # for each Process 
        for i in range(len(unique_processes)): 
            # for each chemsitry 
            for j in range(len(unique_feedstocks)):
                data_1 = Commercial_SAF.loc[ Commercial_SAF['Process'] == unique_processes[i]] 
                data_2 = data_1.loc[Commercial_SAF['Feedstock'] == unique_feedstocks[j]]   
                fig.add_trace(go.Scatter( x        = np.array(data_2[selected_x_axis]), 
                                     y             = np.array(data_2[selected_y_axis]),  
                                     mode          = 'markers', 
                                     name          ="",                                    
                                     marker        = dict(size=marker_size,color=Process_Colors[i],opacity=opacity_ratio,symbol = Feedstock_Markers[j]),
                                     hovertemplate = 'Process: ' + unique_processes[i] + '<br>' + 'Feedstock: ' + unique_feedstocks[j] + '<br>' + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                                     )) 
    elif selected_process != 'All' and  selected_feedstock == 'All':
        # for each Process  
        for j in range(len(unique_feedstocks)):
            data_1   = Commercial_SAF.loc[ Commercial_SAF['Process'] == selected_process] 
            data_2   = data_1.loc[Commercial_SAF['Feedstock'] == unique_feedstocks[j]]  
            i_index  = unique_processes.index(selected_process) 
            fig.add_trace(go.Scatter( x             = np.array(data_2[selected_x_axis]), 
                                 y             = np.array(data_2[selected_y_axis]),  
                                 mode          = 'markers',
                                 name          = "",            
                                 marker        = dict(size=marker_size,color=Process_Colors[i_index ],opacity=opacity_ratio,symbol = Feedstock_Markers[j]),
                                 hovertemplate = 'Process: ' + selected_process + '<br>' + 'Feedstock: ' + unique_feedstocks[j] + '<br>' + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                                 ))
    elif selected_process == 'All' and selected_feedstock != 'All':
        # for each Process 
        for i in range(len(unique_processes)):  
            data_1   = Commercial_SAF.loc[ Commercial_SAF['Process'] == unique_processes[i]] 
            data_2   = data_1.loc[Commercial_SAF['Feedstock'] == selected_feedstock] 
            j_index  = unique_feedstocks.index(selected_feedstock)
            models   = data_2["Model"]
            config   = data_2["Configuration"]
            fig.add_trace(go.Scatter( x             = np.array(data_2[selected_x_axis]), 
                                 y             = np.array(data_2[selected_y_axis]),  
                                 mode          = 'markers',
                                 name          = "",      
                                 marker        = dict(size=marker_size,color=Process_Colors[i],opacity=opacity_ratio,symbol = Feedstock_Markers[j_index]),
                                 hovertemplate = 'Process: ' + unique_processes[i] + '<br>' + 'Feedstock: ' + selected_feedstock + '<br>' + 'Configuration: ' + config + '<br>' + 'Model: ' + models + '<br>' + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                                 ))
    else:
        data_1  = Commercial_SAF.loc[ Commercial_SAF['Process'] == selected_process ] 
        data_2  = data_1.loc[Commercial_SAF['Feedstock'] == selected_feedstock] 
        i_index = unique_processes.index(selected_process)
        j_index = unique_feedstocks.index(selected_feedstock) 
        fig.add_trace(go.Scatter( x        = np.array(data_2[selected_x_axis]), 
                             y             = np.array(data_2[selected_y_axis]),  
                             mode          = 'markers', 
                             name          = "",       
                             marker        = dict(size=marker_size,color=Process_Colors[i_index],opacity=opacity_ratio,symbol = Feedstock_Markers[j_index]),
                             hovertemplate = 'Process: ' + selected_process  + '<br>' + 'Feedstock: ' + selected_feedstock + '<br>' + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                             ))
 
    fig.update_layout(xaxis_title = selected_x_axis,
                       yaxis_title = selected_y_axis,
                       showlegend  = False, 
                       height      = 400,
                       margin      ={'t':0,'l':0,'b':0,'r':0},
                       font=dict(  size=font_size ),  
                       )     
    fig["layout"]["template"] = template 
    return fig 

def generate_saf_spider_plot(Commercial_SAF,saf_1,saf_2,saf_3,switch_off): 
    template      = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]   
    metrics_list  = ['Density @15°C, (g/cm3)' ,
                     'Boiling Point (°C)',
                     'Freezing Point (°C)',
                     'Autoignition Temperature (°C)',
                     'Vapor Pressure @ 20°C',
                     'Upper Explosive Limit',
                     'Lower Explosive Limit']    
    font_size     = 16     
    data_1        = Commercial_SAF.loc[Commercial_SAF['Fuel Name'] == saf_1] 
    vals_1        = np.zeros((1,8))
    vals_1[:,:7]  = np.array(data_1[metrics_list])
    vals_1[0,7]   = vals_1[0,0]    
        
    data_2       = Commercial_SAF.loc[ Commercial_SAF['Fuel Name'] == saf_2]  
    vals_2       = np.zeros((1,8))
    vals_2[:,:7] = np.array(data_2[metrics_list]) 
    vals_2[0,7]  = vals_2[0,0]   

    data_3       = Commercial_SAF.loc[ Commercial_SAF['Fuel Name'] == saf_3]  
    vals_3       = np.zeros((1,8))
    vals_3[:,:7] = np.array(data_3[metrics_list]) 
    vals_3[0,7]  = vals_3[0,0]  
    
    # Stack data
    saf_data = np.vstack((vals_1,vals_2 ))    
    saf_data = np.vstack((saf_data,vals_3))    
     
    scaling = np.array([np.maximum(1000,np.max(saf_data[:,0])*1.05),
                        np.maximum(5,np.max(saf_data[:,1])*1.05) ,
                        np.maximum(5000,np.max(saf_data[:,2])*1.05),
                        np.maximum(20,np.max(saf_data[:,3])*1.05), 
                        np.maximum(50,np.max(saf_data[:,4])*1.05),
                        np.maximum(500,np.max(saf_data[:,5])*1.05),
                        np.maximum(500,np.max(saf_data[:,6])*1.05),
                        np.maximum(1000,np.max(saf_data[:,7])*1.05)])
    scaling = np.tile(scaling[:,None],(1,3)).T
    fig = go.Figure()  
    
    scaled_data = np.zeros((2,len(scaling))) 
    scaled_data = saf_data/scaling 
    fig.add_trace(
                go.Scatterpolar(
                                r    = scaled_data[0],
                                theta=['Density @15°C  <br>  (g/cm3)' ,
                                       'Boiling Point  <br>  (°C)',
                                       'Freezing Point  <br>  (°C)',
                                       'Autoignition Temperature  <br>  (°C)',
                                       'Vapor Pressure  <br>  @ 20°C',
                                       'Upper Explosive <br>  Limit',
                                       'Lower Explosive  <br>  Limit',
                                       'Density @15°C  <br>  (g/cm3)' ],  
                                fill      ='toself', 
                                name      = 'SAF 1', 
                                marker    = None,
                                line=dict(color=px.colors.qualitative.Pastel[0],width = 4),  
                                showlegend=True, 
                                )
                )
    

    fig.add_trace(
                go.Scatterpolar(
                                r    = scaled_data[1],
                                theta=['Density @15°C  <br>  (g/cm3)' ,
                                       'Boiling Point  <br>  (°C)',
                                       'Freezing Point  <br>  (°C)',
                                       'Autoignition Temperature  <br>  (°C)',
                                       'Vapor Pressure  <br>  @ 20°C',
                                       'Upper Explosive <br>  Limit',
                                       'Lower Explosive  <br>  Limit',
                                       'Density @15°C  <br>  (g/cm3)' ],  
                                fill      ='toself', 
                                name      = 'SAF 2', 
                                marker    = None,
                                line=dict(color=px.colors.qualitative.Pastel[1],width = 4),  
                                showlegend=True, 
                                )
                )        


    fig.add_trace(
                go.Scatterpolar(
                                r    = scaled_data[2],
                                theta=['Density @15°C  <br>  (g/cm3)' ,
                                       'Boiling Point  <br>  (°C)',
                                       'Freezing Point  <br>  (°C)',
                                       'Autoignition Temperature  <br>  (°C)',
                                       'Vapor Pressure  <br>  @ 20°C',
                                       'Upper Explosive <br>  Limit',
                                       'Lower Explosive  <br>  Limit',
                                       'Density @15°C  <br>  (g/cm3)' ],  
                                fill      ='toself', 
                                name      = 'SAF 3', 
                                marker    = None,
                                line=dict(color=px.colors.qualitative.Pastel[4],width = 4),  
                                showlegend=True, 
                                )
                )        

    fig.update_layout(height      = 400, 
                       margin={'t':50},
                       font=dict(  size=font_size ),  
                           )            
    fig.update_polars(radialaxis=dict(visible=False,range=[0, 1]))  
    fig["layout"]["template"] = template 
    return fig 


def generate_saf_dev_map(SAF_Development,selected_sector,selected_type,switch_off): 
    template            = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]    
    map_style           = None if switch_off else 'dark'
    map_style           = None if switch_off else 'dark'  
    unique_sectors      = ['Industry', 'Academia', 'Government']
    unique_process      = ['HEFA','Alcohol-to-Jet','Gasification/FT','Catalytic Hydrothermolysis']  # change to process 
    sector_colors       = px.colors.qualitative.Pastel 
    mapbox_access_token = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ" 
    
    fig = go.Figure()
    if selected_sector == 'All' and  selected_type == 'All': 
        for i in range(len(unique_sectors)):
            for j in range(len(unique_process)):
                data_1 = SAF_Development.loc[SAF_Development['Sector'] == unique_sectors[i]] 
                data_2 = data_1.loc[SAF_Development[unique_process[j]] == 1]   
                fig2   = px.scatter_mapbox(data_2, lat="Latitude", lon="Longitude",
                                          hover_name="Entity",
                                          hover_data=["City"],
                                         color_discrete_sequence=[sector_colors[i]], zoom=1 ,)
                fig.add_trace(fig2.data[0])
    
    elif selected_sector == 'All' and  selected_type != 'All': 
        for i in range(len(unique_sectors)): 
            data_1 = SAF_Development.loc[SAF_Development['Sector'] == unique_sectors[i]] 
            data_2 = data_1.loc[SAF_Development[selected_type] == 1]   
            fig2   = px.scatter_mapbox(data_2, lat="Latitude", lon="Longitude",
                                      hover_name="Entity",
                                      hover_data=["City"],
                                     color_discrete_sequence=[sector_colors[i]], zoom=1 ,)
            fig.add_trace(fig2.data[0])
            
    
    elif selected_sector != 'All' and  selected_type == 'All': 
        color_idx = unique_sectors.index(selected_sector)
        for j in range(len(unique_sectors)): 
            data_1 = SAF_Development.loc[SAF_Development['Sector'] == selected_sector] 
            data_2 = data_1.loc[SAF_Development[unique_process[j]] == 1]  
            fig2   = px.scatter_mapbox(data_2, lat="Latitude", lon="Longitude",
                                      hover_name="Entity",
                                      hover_data=["City"],
                                     color_discrete_sequence=[sector_colors[color_idx]], zoom=1 ,)
            fig.add_trace(fig2.data[0])      
    
    else:  
        color_idx = unique_sectors.index(selected_sector)
        data_1    = SAF_Development.loc[SAF_Development['Sector'] == selected_sector] 
        data_2    = data_1.loc[SAF_Development[selected_type] == 1]
        fig2      = px.scatter_mapbox(data_2, lat="Latitude", lon="Longitude",
                                  hover_name="Entity",
                                  hover_data=["City"],
                                 color_discrete_sequence=[sector_colors[color_idx]], zoom=1 ,)
        fig.add_trace(fig2.data[0])          

    
    fig.update_traces(marker={"size": 10})
    fig.update_layout(mapbox_style  = "open-street-map",      
                      showlegend    = False, 
                      height        = 300, 
                      margin        = {'t':0,'l':0,'b':0,'r':0}, 
                      mapbox        = dict( accesstoken=mapbox_access_token,style=map_style,
                                          center=go.layout.mapbox.Center( lat=20, lon= 200 ))  )     

    fig["layout"]["template"] = template 
    return fig 




 
def generate_saf_flight_operations_plots(Flight_Ops,Commercial_SAF,corn_feedstock,soybean_feedstock,selected_fuels,blend_ratios,percent_fuel_use,SAF_production_states_1,SAF_production_states_2,aircraft,flight_range,selected_airpots,percent_adoption,selected_feedstock,month_no,switch_off):  
    template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]    
    font_size            = 16  
    with urlopen('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json') as response:
        counties = json.load(response)  
         
    # Step 1: compute the percentages of different types of fuels used 
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
    if aircraft == 'ATR 72-600': 
        P_max           = 1568083 
        W_0             = 23000   
        thrust_coef     = 0.000017678 # incorrect  
        fuel_volume     = 6410.2
        fuel_economy    = 3.27 
        passengers      = 72 
        CL_cruise       = 0.55   # incorrect  
        CD_cruise       = 0.035  # incorrect  
        S_ref           = 61.0   
        L_div_D         = CL_cruise/CD_cruise 
        
    #elif aircraft == 'Embraer 190': 
        #P_max           =  
        #W_0             = 52290 
        #thrust_coef     = 0.000017678 # incorrect    
        #fuel_volume     = 16629 
        #fuel_economy    = 3.54 
        #passengers      = 100 
        #CL_cruise       = 0.55   # incorrect 
        #CD_cruise       = 0.035  # incorrect 
        #S_ref           = 92.53 
        #L_div_D         = CL_cruise/CD_cruise  
    elif aircraft == "Boeing 737 MAX-8":
        P_max           = 15000000 
        W_0             = 79015.8   
        thrust_coef     = 0.000017678  
        fuel_volume     = 26024 
        fuel_economy    = 2.04 
        passengers      = 162  
        CL_cruise       = 0.55 
        CD_cruise       = 0.035
        S_ref           = 127 
        L_div_D         = CL_cruise/CD_cruise
        
    elif aircraft == 'Airbus A320 neo': 
        P_max           = 15000000 
        W_0             = 78000
        thrust_coef     = 0.000017678  
        fuel_volume     = 26024  
        fuel_economy    = 2.04  
        passengers      = 162
        W_0             = 79015.8  
        CL_cruise       = 0.55 
        CD_cruise       = 0.035
        S_ref           = 127   
        L_div_D         = CL_cruise/CD_cruise
        
    #elif aircraft == "Boeing 787-8": 
        #P_max           =  
        #thrust_coef     = 0.000017678 
        #fuel_volume     = 126206  
        #fuel_economy    = 2.68    
        #passengers      = 238
        #W_0             = 227900   
        #CL_cruise       = 0.055 
        #CD_cruise       = 0.035
        #S_ref           = 377
        #L_div_D         = CL_cruise/CD_cruise
        
    #elif aircraft == 'Airbus A350-1000':
        #P_max          =  
        #thrust_coef    = 0.000017678   
        #fuel_volume    = 166488        # liters 
        #fuel_economy   = 2.85          # L/100km per Pax 
        #passengers     = 327  
        #W_0            = 322050  
        #CL_cruise      = 0.055 
        #CD_cruise      = 0.035
        #S_ref          = 464.3
        #L_div_D        = CL_cruise/CD_cruise      
        
    # Step 4: filter list of routes by range 
    Flight_Ops   =  Flight_Ops[Flight_Ops['Distance (miles)'] < flight_range[0]] 
    miles_to_km  = 1.60934 
    flight_range = Flight_Ops['Distance (miles)']*miles_to_km 
    passengers   = Flight_Ops['Passengers'] 
        
    # Step 6: Filter flight data based on option selected: i.e. top 10, top 20, top 50, all airpots 
    Airport_Routes     = Flight_Ops[['Passengers','Origin Airport','Destination City']]
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
    
    # Step 7: Filter airports that will support SAF and those that wont support SAF 
    mask_1                = Flight_Ops['Origin Airport'].isin(Airport_List)
    SAF_Possible_Airports = Flight_Ops[mask_1]
    Non_SAF_Airports      = Flight_Ops[~mask_1]  
    
    # Step 8: Out of SAF supporting airports, use the percent adoption to detemrine how many flights at that airport will use SAF
    mask_2                      = np.len(0,len(SAF_Possible_Airports), int(len(SAF_Possible_Airports)*percent_adoption ))
    SAF_Airports_Using_SAF      = SAF_Possible_Airports[mask_2]  # flights on SAF 
    SAF_Airports_Using_Jet_A    = SAF_Possible_Airports[~mask_2] 
    Non_SAF_Flights             = pd.concat([Non_SAF_Airports, SAF_Airports_Using_Jet_A] )   # add list of flights from non supporting airports to non-SAF flights  
     
    # Step 9: Get total volume of each SAF required at the airports 
    total_fuel_volume_required = np.sum(np.float(SAF_Airports_Using_SAF['Total Fuel Per Route (Gal)'][1:]))
    fuel_volumes               = cumulative_fuel_use*total_fuel_volume_required  
    
    # Step 10: Sort SAF's by feedstock and sum all fuel volumes based on fuel (FUTURE WORK: determine how much corn is needed for EACH type of process)
    mask                       = Commercial_SAF['Feedstock'].isin(selected_fuels)
    fuels_used                 = Commercial_SAF[mask]    
    fuels_used['Fuels Used']   = fuel_volumes 
    relative_crop_volume       = fuels_used['Fuels Used']*fuels_used['Fuel-to-Crop-factor']
    fuels_used['Crop-Volume']  = relative_crop_volume
    df2                        = fuels_used.groupby(['Feedstock']).sum()  
    required_crop_vol          = df2[selected_feedstock]['Crop-Volume'] 
    
    # Step 11: Determine how many states will source the feedstock under consideration 
    if selected_feedstock == "Corn": 
        crop_data = corn_feedstock
    elif selected_feedstock == "Barley": 
        crop_data = soybean_feedstock 
    fips_data  = "{:02d}".format(crop_data['State ANSI']) +  "{:03d}".format(int(crop_data['County ANSI']))
    crop_data['FIPS'] = fips_data 
    
    # Step 12: Filter out states have not been selected and get number of states  
    crop_data_states      = SAF_production_states_1 + SAF_production_states_2 
    mask                  = crop_data['State'].isin(crop_data_states)
    feedstock_states      = crop_data[mask]    
    non_feedstock_states  = crop_data[~mask]     

    # Step 13: Randomize tracts in terms of crop area/usage then Recursively add rows until requied volume is met 
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
     
    # Step 14: Plot Land Use 
    Land_Use  = pd.concat([Used_Feedstock, non_feedstock_states] )   
    fips      = list(Land_Use['FIPS'])  
    Land_Use['FIPS'] = ["%05d" % i for i in fips] 
    fig_4 = px.choropleth(Land_Use, geojson=counties, locations='FIPS',color = Feedstock_Usage,
                           color_continuous_scale="RdYlBu_r", 
                           hover_data=["Name","State", Feedstock_Usage],
                           scope='usa',
                           range_color=(0, 90),)   
    fig_4.update_layout(coloraxis_colorbar=dict(title=" "),
                                     coloraxis_colorbar_x=0.85, 
                                     height    = 400, 
                                     margin={'t':0,'l':0,'b':0,'r':0}, )  
    fig_4.update_coloraxes( colorbar_tickvals= np.linspace(0,90,11),
                                        colorbar_tickfont_size=font_size)   
    

    fig_4["layout"]["template"] = template
     
    return fig_4

def generate_saf_flight_ops_map(Flight_Ops,Commercial_SAF,corn_feedstock,soybean_feedstock,selected_fuels,blend_ratios,percent_fuel_use,SAF_production_states_1,SAF_production_states_2,aircraft,selected_airpots,percent_adoption,selected_feedstock,month_no,switch_off):  
    template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]  
    mapbox_access_token  = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ"
    map_style            = None if switch_off else 'dark'   
    

    fuel_reserve        = 1800 # 30 mins 
    reserve_speed       = 250 # m/s
    air_density         = 8.21392E-3
    Jet_A_density       = 775.0 #  kg/m3 
    
    # Step 3: Based on aircraft's fuel consumption, determine the total volume of fuel used oer flight 
    if aircraft == "Boeing 787-8": 
        C_T             = 0.000017678 
        max_fuel_volume = 126206  # liters 
        fuel_economy = 2.68 # L/100km per Pax 
        passengers      = 238
        W_0             = 227900   
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 377
    if aircraft == "Boeing 737 MAX-8":
        C_T             = 0.000017678  
        max_fuel_volume = 26024 # liters 
        fuel_economy = 2.04 # L/100km per Pax 
        passengers      = 162
        W_0             = 79015.8  
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 127 
    elif aircraft == 'Airbus A350-1000': 
        C_T             = 0.000017678   
        max_fuel_volume = 166488 
        fuel_economy = 2.85 # L/100km per Pax 
        passengers      = 327  
        W_0             = 322050  
        CL_cruise       = 0.055 
        CD_cruise       = 0.035
        S_ref           = 464.3
    elif aircraft == 'Airbus A320 neo': 
        C_T             = 0.000017678  
        max_fuel_volume = 27200  # liters 
        fuel_economy = 1.94   # L/100km per Pax 
        passengers      = 180
        W_0             = 79015
        CL_cruise       = 0.055 
        CD_cruise       = 0.035 
        S_ref           = 122.6   
        
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # Compute distance between departure and destimation points
    # ----------------------------------------------------------------------------------------------------------------------------------------------
    Routes_Mo                    = Flight_Ops[Flight_Ops['Month'] == month_no+1 ]  
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
    Routes_Mo['Range']           = distances.tolist() 

    
    reserve_distance    = reserve_speed*fuel_reserve
    fuel_volume_L       = fuel_economy*passengers*(distances+ reserve_distance)/(1000*100)
    fuel_volume_m_3     = fuel_volume_L*0.001  
    W_f                 = Jet_A_density*fuel_volume_m_3
    W_1                 = W_0*percent_fuel_use - W_f
    Range               = (1/C_T)*(np.sqrt(CL_cruise)/CD_cruise)*np.sqrt(2/(air_density*S_ref))*(np.sqrt(W_0) -np.sqrt(W_1) ) 
    Infeasible_Routes   = Routes_Mo[Routes_Mo['Range'] > Range ]  
    Feasible_Routes     = Routes_Mo[Routes_Mo['Range'] < Range ]  
    
    #================================================================================================================================================  
    # Plot Flight_Ops 
    #================================================================================================================================================     
    # Flight_Ops   
    fig_5                = go.Figure()
    airport_marker_size  = 5
    airport_marker_color = "white"
 
    # Flight Paths
    lons       = np.empty(3 * len(Infeasible_Routes))
    lons[::3]  = Infeasible_Routes['Origin Longitude (Deg.)']
    lons[1::3] = Infeasible_Routes['Destination Longitude (Deg.)']
    lons[2::3] = None
    lats       = np.empty(3 * len(Infeasible_Routes))
    lats[::3]  = Infeasible_Routes['Origin Latitude (Deg.)']
    lats[1::3] = Infeasible_Routes['Destination Latitude (Deg.)']
    lats[2::3] = None    
  
    fig_5.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            line = dict(width = 0.1,color = 'coral'), ))
    
    lons       = np.empty(3 * len(Feasible_Routes))
    lons[::3]  = Feasible_Routes['Origin Longitude (Deg.)']
    lons[1::3] = Feasible_Routes['Destination Longitude (Deg.)']
    lons[2::3] = None
    lats       = np.empty(3 * len(Feasible_Routes))
    lats[::3]  = Feasible_Routes['Origin Latitude (Deg.)']
    lats[1::3] = Feasible_Routes['Destination Latitude (Deg.)']
    lats[2::3] = None    
  
    fig_5.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            line = dict(width = 2,color = "aquamarine"), ))
  
    # Airports  
    fig_5.add_trace(go.Scattergeo( 
        lon = Flight_Ops['Destination Longitude (Deg.)'],
        lat = Flight_Ops['Destination Latitude (Deg.)'], 
        text = Flight_Ops['Destination City'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color, ))) 

    fig_5.add_trace(go.Scattergeo( 
        lon = Flight_Ops['Origin Longitude (Deg.)'],
        lat = Flight_Ops['Origin Latitude (Deg.)'], 
        text = Flight_Ops['Origin City'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color,))) 
     
    # Flight Paths 
    fig_5.update_layout(mapbox_style  = "open-street-map",      
                      showlegend    = False, 
                      height        = 400, 
                      geo_scope     ='usa',
                      margin        = {'t':0,'l':0,'b':0,'r':0},  
                      mapbox        = dict( accesstoken=mapbox_access_token,style=map_style,
                                            center=go.layout.mapbox.Center( lat=30, lon= 230 )))   
     
    #================================================================================================================================================  
    #  
    #================================================================================================================================================     


    #================================================================================================================================================  
    #  
    #================================================================================================================================================     



    #================================================================================================================================================  
    #  
    #================================================================================================================================================     


    #================================================================================================================================================  
    #  
    #================================================================================================================================================     


    fig_5["layout"]["template"] = template 
             
   
        
    return fig_5

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# SAF Bar
# ---------------------------------------------------------------------------------------------------------------------------------------------------
def generate_saf_slider_bar(Commercial_SAF,fuel_values): 
    
    v    = [0] 
    v    = v + fuel_values
    v    = v + [100]
    vals = np.array(v)
    fuel_quantity = np.diff(vals)
    df = pd.DataFrame({'Quantity': list(fuel_quantity),
                       'Source'  : ['Fuel','Fuel','Fuel','Fuel'] ,
                       'Fuel'    : list(Commercial_SAF['Fuel Name'][1:].unique()) }) 
    fig = px.histogram(df, y="Source",
                   x="Quantity", color="Fuel",
                   barnorm='percent')
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                      xaxis_title=None,
                      yaxis_title=None, 
                      xaxis     =  {'title': 'x-label','visible': False,'showticklabels': True},
                      yaxis     =  {'title': 'y-label','visible': False,'showticklabels': False},
                              )
    
    return fig 



 
if __name__ == '__main__':
    main()