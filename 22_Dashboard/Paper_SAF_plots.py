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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.colors as mcolors
 

def main():
   
    separator                      = os.path.sep
    SAF_filename                   = 'Data' + separator + 'Technology' + separator +  'Technology_Data.xlsx'
    routes_filename                = 'Data' + separator +  'Air_Travel' + separator + 'American_Airlines_Flight_Ops_and_Climate.xlsx'
    crops_filename                 = 'Data' + separator +  'Crops'     + separator + 'All_Crops_2017.xlsx' 
    
    SAT_data                       = pd.read_excel(SAF_filename,sheet_name=['Commercial_SAF']) 
    Commercial_SAF                 = SAT_data['Commercial_SAF']  
    a                              = Commercial_SAF['Brand']  
    b                              = Commercial_SAF['Fuel Type']
    c                              = Commercial_SAF['Process']
    d                              = Commercial_SAF['Source']
    e                              = Commercial_SAF['Feedstock']
    Commercial_SAF["Fuel Name"]    = a + ' ' + b + ' from ' + c  + ' using ' + d + '-' + e 
    Flight_Ops                     = pd.read_excel(routes_filename,sheet_name=['Sheet1']) 
    Flight_Ops                     = Flight_Ops['Sheet1']     
    feedstocks                     = pd.read_excel(crops_filename,sheet_name=['Corn','Soybean','Canola','Sunflower','Sorghum','Wheat'])  

    fuel_values             = [40,60,75,95]   
    saf_1                   = Commercial_SAF['Fuel Name'][0] 
    saf_2                   = Commercial_SAF['Fuel Name'][12] 
    saf_3                   = Commercial_SAF['Fuel Name'][3]  
    saf_5                   = Commercial_SAF['Fuel Name'][11]  
    saf_4                   = Commercial_SAF['Fuel Name'][8]     
    fuels_1                   = [saf_1,saf_2]  
    fuels_2                  = [saf_3,saf_4,saf_5]
    fuels = fuels_1 + fuels_2
    #fig_x = generate_saf_slider_bar(Commercial_SAF,fuel_values,fuels )
    #fig_x.show()
    
    selected_process   = 'All'
    selected_feedstock = 'All'
    selected_x_axis    = list(Commercial_SAF.columns.values)[10:39][7]  
    selected_y_axis    = list(Commercial_SAF.columns.values)[10:39][9]  
    switch_off         = False     
    #fig_1 = generate_saf_scatter_plot(Commercial_SAF,selected_process,selected_feedstock,selected_x_axis,selected_y_axis,switch_off)
    #fig_1.show()

    saf_1           = Commercial_SAF['Fuel Name'][4] 
    saf_2           = Commercial_SAF['Fuel Name'][1] 
    saf_3           = Commercial_SAF['Fuel Name'][3] 
    switch_off         = False
    #fig_2 = generate_saf_spider_plot(Commercial_SAF,saf_1,saf_2,saf_3,switch_off)
    #fig_2.show()

    selected_feedstock = 'UCO'
    selected_process   = 'HEFA'
    switch_off      = False
    #fig_3 = generate_saf_dev_map(Commercial_SAF,selected_feedstock,selected_process,switch_off)
    #fig_3.show()

    saf_1                           = Commercial_SAF['Fuel Name'][0] 
    saf_2                           = Commercial_SAF['Fuel Name'][1] 
    #saf_3                           = Commercial_SAF['Fuel Name'][2] 
    #saf_4                           = Commercial_SAF['Fuel Name'][3]     
    #selected_fuels                  = [saf_1,saf_2,saf_3,saf_4]
    selected_fuels                  = [saf_1,saf_2]
    #blend_ratios                    = [100,40,50,25]
    #percent_fuel_use                = [75,80,95] 
    blend_ratios                    = [100,40]
    State_List_1                    = ["California","Colorado"]  
    State_List_2                    = ["Georgia","Idaho","Illinois", "Indiana","Iowa","Kansas","Kentucky","Louisiana"] 
    State_List_3                    = ["Missouri", "Nebraska"] 
    State_List_4                    = ["Ohio","Oklahoma"] 
    State_List_5                    = [""]         
    aircraft                        = "Boeing 787-8"
    selected_airpots                = "All Airports" # ["Top 5 Airports","Top 10 Airports","Top 20 Airports", "Top 50 Airports",  "All Airports"]

    month_no                        = 1
    selected_feedstock              = 'Canola'
    switch_off                      = False  
    feedstock_producing_states     = State_List_1+State_List_2+State_List_3+State_List_4+State_List_5
    vector, w_SAF, w_o_SAF, CASM_SAF, CASM_jet_A = generate_saf_flight_operations_plots(Flight_Ops,Commercial_SAF,feedstocks,selected_fuels,selected_feedstock,blend_ratios,feedstock_producing_states,aircraft, selected_airpots,month_no,switch_off)
    generate_plot(vector, w_SAF, w_o_SAF, CASM_SAF, CASM_jet_A)
    return 
     
## ---------------------------------------------------------------------------------------------------------------------------------------------------
## SAF Plots 
## ---------------------------------------------------------------------------------------------------------------------------------------------------
#def generate_saf_scatter_plot(Commercial_SAF,selected_process,selected_feedstock,selected_x_axis,selected_y_axis,switch_off): 
    ##template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]  
    #unique_processes     = list(Commercial_SAF['Process'][1:].unique())
    #unique_feedstocks    = list(Commercial_SAF['Feedstock'][1:].unique())
    #marker_size          = 15
    #opacity_ratio        = 0.8 if switch_off else 1.0
    #font_size            = 16 

    ## Process Colors: greenyellow, aquamarine, paleturquoise, lightcoral, yellow, lavender ,thistle ,orangered   
    #Process_Colors      = px.colors.qualitative.Pastel 
    #Feedstock_Markers  = ['square','x','circle','cross','diamond','triangle-up','triangle-down','star','hourglass'] 
    #fig = go.Figure()
    #if selected_process == 'All' and  selected_feedstock == 'All':
        ## for each Process 
        #for i in range(len(unique_processes)): 
            ## for each chemsitry 
            #for j in range(len(unique_feedstocks)):
                #data_1 = Commercial_SAF.loc[ Commercial_SAF['Process'] == unique_processes[i]] 
                #data_2 = data_1.loc[Commercial_SAF['Feedstock'] == unique_feedstocks[j]]   
                #fig.add_trace(go.Scatter( x        = np.array(data_2[selected_x_axis]), 
                                     #y             = np.array(data_2[selected_y_axis]),  
                                     #mode          = 'markers', 
                                     #name          ="",                                    
                                     #marker        = dict(size=marker_size,color=Process_Colors[i],opacity=opacity_ratio,symbol = Feedstock_Markers[j]),
                                     #hovertemplate = 'Process: ' + unique_processes[i] + '<br>' + 'Feedstock: ' + unique_feedstocks[j] + '<br>' + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                                     #)) 
    #elif selected_process != 'All' and  selected_feedstock == 'All':
        ## for each Process  
        #for j in range(len(unique_feedstocks)):
            #data_1   = Commercial_SAF.loc[ Commercial_SAF['Process'] == selected_process] 
            #data_2   = data_1.loc[Commercial_SAF['Feedstock'] == unique_feedstocks[j]]  
            #i_index  = unique_processes.index(selected_process) 
            #fig.add_trace(go.Scatter( x             = np.array(data_2[selected_x_axis]), 
                                 #y             = np.array(data_2[selected_y_axis]),  
                                 #mode          = 'markers',
                                 #name          = "",            
                                 #marker        = dict(size=marker_size,color=Process_Colors[i_index ],opacity=opacity_ratio,symbol = Feedstock_Markers[j]),
                                 #hovertemplate = 'Process: ' + selected_process + '<br>' + 'Feedstock: ' + unique_feedstocks[j] + '<br>' + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                                 #))
    #elif selected_process == 'All' and selected_feedstock != 'All':
        ## for each Process 
        #for i in range(len(unique_processes)):  
            #data_1   = Commercial_SAF.loc[ Commercial_SAF['Process'] == unique_processes[i]] 
            #data_2   = data_1.loc[Commercial_SAF['Feedstock'] == selected_feedstock] 
            #j_index  = unique_feedstocks.index(selected_feedstock)
            #models   = data_2["Model"]
            #config   = data_2["Configuration"]
            #fig.add_trace(go.Scatter( x             = np.array(data_2[selected_x_axis]), 
                                 #y             = np.array(data_2[selected_y_axis]),  
                                 #mode          = 'markers',
                                 #name          = "",      
                                 #marker        = dict(size=marker_size,color=Process_Colors[i],opacity=opacity_ratio,symbol = Feedstock_Markers[j_index]),
                                 #hovertemplate = 'Process: ' + unique_processes[i] + '<br>' + 'Feedstock: ' + selected_feedstock + '<br>' + 'Configuration: ' + config + '<br>' + 'Model: ' + models + '<br>' + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                                 #))
    #else:
        #data_1  = Commercial_SAF.loc[ Commercial_SAF['Process'] == selected_process ] 
        #data_2  = data_1.loc[Commercial_SAF['Feedstock'] == selected_feedstock] 
        #i_index = unique_processes.index(selected_process)
        #j_index = unique_feedstocks.index(selected_feedstock) 
        #fig.add_trace(go.Scatter( x        = np.array(data_2[selected_x_axis]), 
                             #y             = np.array(data_2[selected_y_axis]),  
                             #mode          = 'markers', 
                             #name          = "",       
                             #marker        = dict(size=marker_size,color=Process_Colors[i_index],opacity=opacity_ratio,symbol = Feedstock_Markers[j_index]),
                             #hovertemplate = 'Process: ' + selected_process  + '<br>' + 'Feedstock: ' + selected_feedstock + '<br>' + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                             #))
 
    #fig.update_layout(xaxis_title = selected_x_axis,
                       #yaxis_title = selected_y_axis,
                       #showlegend  = False, 
                       #height      = 400,
                       #margin      ={'t':0,'l':0,'b':0,'r':0},
                       #font=dict(  size=font_size ),  
                       #)     
    ##fig["layout"]["template"] = template 
    #return fig 

#def generate_saf_spider_plot(Commercial_SAF,saf_1,saf_2,saf_3,switch_off): 
    ##template      = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]   
    #metrics_list  = ['Density @15°C, (g/cm3)' ,
                     #'Boiling Point (°C)',
                     #'Freezing Point (°C)',
                     #'Autoignition Temperature (°C)',
                     #'Vapor Pressure @ 20°C (kPa)',
                     #'Upper Explosive Limit (%)',
                     #'Lower Explosive Limit (%)']    
    #font_size     = 16     
    #data_1        = Commercial_SAF.loc[Commercial_SAF['Fuel Name'] == saf_1] 
    #vals_1        = np.zeros((1,8))
    #vals_1[:,:7]  = np.array(data_1[metrics_list])
    #vals_1[0,7]   = vals_1[0,0]    
        
    #data_2       = Commercial_SAF.loc[ Commercial_SAF['Fuel Name'] == saf_2]  
    #vals_2       = np.zeros((1,8))
    #vals_2[:,:7] = np.array(data_2[metrics_list]) 
    #vals_2[0,7]  = vals_2[0,0]   

    #data_3       = Commercial_SAF.loc[ Commercial_SAF['Fuel Name'] == saf_3]  
    #vals_3       = np.zeros((1,8))
    #vals_3[:,:7] = np.array(data_3[metrics_list]) 
    #vals_3[0,7]  = vals_3[0,0]  
    
    ## Stack data
    #saf_data = np.vstack((vals_1,vals_2 ))    
    #saf_data = np.vstack((saf_data,vals_3))    
     
    #scaling = np.array([np.maximum(1,np.max(saf_data[:,0])*1.05),
                        #np.maximum(400,np.max(saf_data[:,1])*1.05) ,
                        #np.minimum(-100,np.max(saf_data[:,2])*1.05),
                        #np.maximum(250,np.max(saf_data[:,3])*1.05), 
                        #np.maximum(4,np.max(saf_data[:,4])*1.05),
                        #np.maximum(15,np.max(saf_data[:,5])*1.05),
                        #np.maximum(1,np.max(saf_data[:,6])*1.05),
                        #np.maximum(1,np.max(saf_data[:,7])*1.05)])
    #scaling = np.tile(scaling[:,None],(1,3)).T
    #fig = go.Figure()  
    
    #scaled_data = np.zeros((2,len(scaling))) 
    #scaled_data = saf_data/scaling 
    #fig.add_trace(
                #go.Scatterpolar(
                                #r    = scaled_data[0],
                                #theta=['Density @15°C  <br>  (g/cm3)' ,
                                       #'Boiling Point  <br>  (°C)',
                                       #'Freezing Point  <br>  (°C)',
                                       #'Autoignition  <br> Temperature (°C)',
                                       #'Vapor Pressure  <br>  @ 20°C (kPa)',
                                       #'Upper Explosive <br>  Limit (%)',
                                       #'Lower Explosive  <br>  Limit (%)',
                                       #'Density @15°C  <br>  (g/cm3)' ],  
                                #fill      ='toself', 
                                #name      = 'SAF 1', 
                                #marker    = None,
                                #line=dict(color=px.colors.qualitative.Pastel[0],width = 4),  
                                #showlegend=True, 
                                #)
                #)
    

    #fig.add_trace(
                #go.Scatterpolar(
                                #r    = scaled_data[1],
                                #theta=['Density @15°C  <br>  (g/cm3)' ,
                                       #'Boiling Point  <br>  (°C)',
                                       #'Freezing Point  <br>  (°C)',
                                       #'Autoignition  <br> Temperature (°C)',
                                       #'Vapor Pressure  <br>  @ 20°C (kPa)',
                                       #'Upper Explosive <br>  Limit (%)',
                                       #'Lower Explosive  <br>  Limit (%)',
                                       #'Density @15°C  <br>  (g/cm3)' ], 
                                #fill      ='toself', 
                                #name      = 'SAF 2', 
                                #marker    = None,
                                #line=dict(color=px.colors.qualitative.Pastel[1],width = 4),  
                                #showlegend=True, 
                                #)
                #)        


    #fig.add_trace(
                #go.Scatterpolar(
                                #r    = scaled_data[2],
                                #theta=['Density @15°C  <br>  (g/cm3)' ,
                                       #'Boiling Point  <br>  (°C)',
                                       #'Freezing Point  <br>  (°C)',
                                       #'Autoignition  <br> Temperature (°C)',
                                       #'Vapor Pressure  <br>  @ 20°C (kPa)',
                                       #'Upper Explosive <br>  Limit (%)',
                                       #'Lower Explosive  <br>  Limit (%)',
                                       #'Density @15°C  <br>  (g/cm3)' ],
                                #fill      ='toself', 
                                #name      = 'SAF 3', 
                                #marker    = None,
                                #line=dict(color=px.colors.qualitative.Pastel[4],width = 4),  
                                #showlegend=True, 
                                #)
                #)        

    #fig.update_layout(height      = 400, 
                       #margin={'t':50},
                       #font=dict(  size=font_size ),  
                           #)            
    #fig.update_polars(radialaxis=dict(visible=False,range=[0, 1]))  
    ##fig["layout"]["template"] = template 
    #return fig 


#def generate_saf_dev_map(Commercial_SAF,selected_feedstock,selected_process,switch_off): 
    ##template            = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]    
    #map_style           = None if switch_off else 'dark'
    #map_style           = None if switch_off else 'dark'  
    #unique_feedstocks   = list(Commercial_SAF['Feedstock'][1:].unique())
    #unique_process      = list(Commercial_SAF['Process'][1:].unique())   
    #sector_colors       = px.colors.qualitative.Pastel 
    #mapbox_access_token = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ" 
    
    #fig = go.Figure()
    #if selected_feedstock == 'All' and  selected_process == 'All': 
        #for i in range(len(unique_feedstocks)):
            #for j in range(len(unique_process)):
                #data_1 = Commercial_SAF.loc[Commercial_SAF['Feedstock'] == unique_feedstocks[i]] 
                #data_2 = data_1.loc[Commercial_SAF['Process'] == unique_process[j]]   
                #fig2   = px.scatter_mapbox(data_2, lat="Latitude", lon="Longitude",
                                          #hover_name="Fuel Name",
                                          #hover_data=['Feedstock','Process','Source','Maximum Blend Ratio','LCA Value'],
                                         #color_discrete_sequence=[sector_colors[i]], zoom=1 ,)
                #fig.add_trace(fig2.data[0])
    
    #elif selected_feedstock == 'All' and  selected_process != 'All': 
        #for i in range(len(unique_feedstocks)): 
            #data_1 = Commercial_SAF.loc[Commercial_SAF['Feedstock'] == unique_feedstocks[i]] 
            #data_2 = data_1.loc[Commercial_SAF['Process'] == selected_process]   
            #fig2   = px.scatter_mapbox(data_2, lat="Latitude", lon="Longitude",
                                      #hover_name="Fuel Name",
                                      #hover_data=['Feedstock','Process','Source','Maximum Blend Ratio','LCA Value'],
                                     #color_discrete_sequence=[sector_colors[i]], zoom=1 ,)
            #fig.add_trace(fig2.data[0])
            
    
    #elif selected_feedstock != 'All' and  selected_process == 'All': 
        #color_idx = unique_feedstocks.index(selected_feedstock)
        #for j in range(len(unique_process)): 
            #data_1 = Commercial_SAF.loc[Commercial_SAF['Feedstock'] == selected_feedstock] 
            #data_2 = data_1.loc[Commercial_SAF['Process'] == unique_process[j]]  
            #fig2   = px.scatter_mapbox(data_2, lat="Latitude", lon="Longitude",
                                      #hover_name="Fuel Name",
                                      #hover_data=['Feedstock','Process','Source','Maximum Blend Ratio','LCA Value'],
                                     #color_discrete_sequence=[sector_colors[color_idx]], zoom=1 ,)
            #fig.add_trace(fig2.data[0])      
    
    #else:  
        #color_idx = unique_feedstocks.index(selected_feedstock)
        #data_1    = Commercial_SAF.loc[Commercial_SAF['Feedstock'] == selected_feedstock] 
        #data_2    = data_1.loc[Commercial_SAF['Process'] == selected_process]
        #fig2      = px.scatter_mapbox(data_2, lat="Latitude", lon="Longitude",
                                  #hover_name="Fuel Name",
                                  #hover_data=["City"],
                                 #color_discrete_sequence=[sector_colors[color_idx]], zoom=1 ,)
        #fig.add_trace(fig2.data[0])          

    
    #fig.update_traces(marker={"size": 10})
    #fig.update_layout(mapbox_style  = "open-street-map",      
                      #showlegend    = False, 
                      #height        = 300, 
                      #margin        = {'t':0,'l':0,'b':0,'r':0}, 
                      #mapbox        = dict( accesstoken=mapbox_access_token,style=map_style,
                                          #center=go.layout.mapbox.Center( lat=20, lon= 200 ))  )     

    ##fig["layout"]["template"] = template 
    #return fig 

 
def generate_saf_flight_operations_plots(Flight_Ops,Commercial_SAF,feedstocks,selected_fuels,selected_feedstock,blend_ratios,feedstock_producing_states,aircraft,selected_airpots,month_no,switch_off):  
    
    vector      = np.arange(0,101,1)
    CASM_jet_A  = np.zeros(len(vector))
    CASM_SAF    = np.zeros(len(vector))
    w_SAF       = np.zeros(len(vector))
    w_o_SAF     = np.zeros(len(vector))
    
    for w in range(len(vector)):
        #percent_fuel_use     = [vector[w]]  
        percent_fuel_use      = [75]
        SAF_dollars_per_gal   = 10
        percent_adoption      = vector[w]   
        #template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]    
        #map_style            = None if switch_off else 'dark'     
        #font_size            = 16  
        #mapbox_access_token  = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ"
        #with urlopen('https://raw.githubusercontent.com/plotly/datasets/master/geojson-counties-fips.json') as response:
            #counties = json.load(response)  
             
        # Step 1: compute the percentages of different types of fuels used 
        fuel_percentages_list = [0]
        fuel_percentages_list += percent_fuel_use
        fuel_percentages_list += [100]
        fuels_percentages     = np.diff(np.array(fuel_percentages_list))/100  
        
        
        # Step 2: determine the percentage of neat (pure) saf and Jet-A1 using blending ratios
        if Commercial_SAF['Fuel Name'][0] not in selected_fuels: # 
            selected_fuels    = [Commercial_SAF['Fuel Name'][0]] + selected_fuels
            fuels_percentages = np.hstack((np.array([0]),fuels_percentages))
            blend_ratios      = [0] + blend_ratios
        mask                = Commercial_SAF['Fuel Name'].isin(selected_fuels)
        fuels_used          = Commercial_SAF[mask] 
        
        num_fuels = len(selected_fuels)
        cumulative_fuel_use   = np.zeros(num_fuels) 
        SAF_LCA_val           = np.zeros(num_fuels) 
        Jet_A_LCA_val         = np.zeros(num_fuels)  
        
        for i in range(1,num_fuels):
            # loop through fuels and get percentage of fuel used by each type
            cumulative_fuel_use[i]  = fuels_percentages[i]* blend_ratios[i]/100
            SAF_LCA_val[i]          = Commercial_SAF[Commercial_SAF['Fuel Name'] == selected_fuels[i]]['LCA Value']
        cumulative_fuel_use[0]      = 1 - np.sum(cumulative_fuel_use[1:])
        SAF_LCA_val[0]   = 89
        Jet_A_LCA_val[0] = 89
        
        # Step 6: Filter flight data based on option selected: i.e. top 10, top 20, top 50, all airpots 
        Airport_Routes     = Flight_Ops[['Passengers','Origin Airport','Destination City']]
        Cumulative_Flights = Airport_Routes.groupby('Origin Airport', as_index=False).sum()[Airport_Routes.columns]
        if  selected_airpots == "Top 5 Airports":
            Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(5) 
        if  selected_airpots == "Top 10 Airports":
            Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(10) 
        elif  selected_airpots =="Top 20 Airports":
            Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(20) 
        elif  selected_airpots == "Top 50 Airports":
            Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(50) 
        elif  selected_airpots == "All Airports":
            Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False) 
        Airport_List = list(Busiest_Airports['Origin Airport'])
        
        # Step 7: Filter airports that will support SAF and those that wont support SAF 
        mask_1                = Flight_Ops['Origin Airport'].isin(Airport_List)
        SAF_Airports          = Flight_Ops[mask_1]
        Non_SAF_Airports      = Flight_Ops[~mask_1]  
        
        # Step 8: Out of SAF supporting airports, use the percent adoption to determine how many flights at that airport will use SAF 
        Flight_at_SAF_Airports_Using_SAF      = SAF_Airports.sample(frac=(percent_adoption/100))
        Flight_at_SAF_Airports_Using_Jet_A    = SAF_Airports[~SAF_Airports.index.isin(Flight_at_SAF_Airports_Using_SAF.index)]
        Non_SAF_Flights                       = pd.concat([Non_SAF_Airports, Flight_at_SAF_Airports_Using_Jet_A] )   # add list of flights from non supporting airports to non-SAF flights  
         
        # Step 9: Get total volume of each SAF required at the airports 
        total_fuel_volume_required = np.sum(np.array(Flight_at_SAF_Airports_Using_SAF['Total Fuel Per Route (Gal)']))
        fuel_volumes               = cumulative_fuel_use*total_fuel_volume_required  
        
        # Step 10: Sort SAF's by feedstock and sum all fuel volumes based on fuel     
        fuels_used['Total Fuel Volume']   = fuel_volumes 
        relative_crop_area                = fuels_used['Total Fuel Volume']/fuels_used['SAF Gallons per Acre']
        fuels_used['Requires Acres']      = relative_crop_area
        select_fuels                      = fuels_used[['Source','Total Fuel Volume','Requires Acres']].groupby('Source', as_index=False).sum()[fuels_used[['Source','Total Fuel Volume','Requires Acres']].columns]
        required_crop_area                = np.array(select_fuels.loc[select_fuels['Source'] == selected_feedstock]['Requires Acres'])
        
        # Step 11: Determine how many states will source the feedstock under consideration  
        crop_data  = feedstocks[selected_feedstock]   
        crop_data['FIPS'] = crop_data['FIPS'].apply('{:0>5}'.format)
        
        # Step 12: Filter out states have not been selected and get number of states   
        mask                  = crop_data['State'].isin(feedstock_producing_states)
        feedstock_states      = crop_data[mask]    
        non_feedstock_states  = crop_data[~mask]  
        non_feedstock_states["Feedstock_Usage"] = np.ones(len(non_feedstock_states))
    
        # Step 13: Randomize tracts in terms of crop area/usage then Recursively add rows until requied volume is met 
        Used_Feedstock        = feedstock_states.sample(frac = 1)
        tract_use_flag        =  np.ones(len(feedstock_states))
        Used_Feedstock["Feedstock_Usage"] = tract_use_flag 
        idx       = 0
        total_vol = 0
        if len(required_crop_area) == 0:
            RCA  = 0
        else:
            RCA  = required_crop_area[0]  
        available_tracts = len(Used_Feedstock)
        while total_vol<RCA: 
            total_vol += Used_Feedstock.loc[Used_Feedstock.index[idx]]['Acres Harvested']
            Used_Feedstock["Feedstock_Usage"][Used_Feedstock.index[idx]] = 0.3
            idx += 1      
            if available_tracts == idx:
                total_vol = 1E9  
            
        # Steps 14 and 15 Determine Cost per Seat Mile and Emissions  
        #CASM_jet_A  = np.zeros(12) 
        #CASM_SAF    = np.zeros(12) 
        #w_SAF       = np.zeros(12) 
        #w_o_SAF     = np.zeros(12) 
        gallons_to_Liters  = 3.78541
        Non_SAF_Flights_Mo                  = Non_SAF_Flights.loc[Non_SAF_Flights['Month'] == 5 ]   
        Flight_at_SAF_Airports_Using_SAF_Mo = Flight_at_SAF_Airports_Using_SAF.loc[Flight_at_SAF_Airports_Using_SAF['Month'] == 5 ]   
    
        total_SAF_fuel_volume_required_mo   = np.sum(np.array(Flight_at_SAF_Airports_Using_SAF_Mo['Total Fuel Per Route (Gal)']))  
        SAF_volumes_mo                      = cumulative_fuel_use*total_SAF_fuel_volume_required_mo 
        Jet_A_fuel_volume_required_mo       = np.sum(np.array(Non_SAF_Flights_Mo['Total Fuel Per Route (Gal)']))    
        
        if len(Non_SAF_Flights_Mo ) == 0:
            pass
        else:
            ASM_jet_A                           = np.sum(Non_SAF_Flights_Mo['Distance (miles)'] * Non_SAF_Flights_Mo['Passengers'])
            Total_Fuel_Cost_jet_A               = np.sum(Non_SAF_Flights_Mo['Fuel Cost']) 
            CASM_jet_A[w]                     = 100*Total_Fuel_Cost_jet_A/ASM_jet_A     
       
        if len(Flight_at_SAF_Airports_Using_SAF_Mo)  == 0:
            pass
        else:    
            ASM_SAF                             = np.sum(Flight_at_SAF_Airports_Using_SAF_Mo['Distance (miles)'] * Flight_at_SAF_Airports_Using_SAF_Mo['Passengers']) 
            Total_Fuel_Cost_SAF                 = np.sum(Flight_at_SAF_Airports_Using_SAF_Mo['Total Fuel Per Route (Gal)'] ) * SAF_dollars_per_gal 
            CASM_SAF[w]                       = 100*Total_Fuel_Cost_SAF/ASM_SAF   
             
        w_SAF[w]   = np.sum(fuels_used['LCEF (gCO2e/MJ)']*fuels_used['Volumetric Energy Density (MJ/L)']*SAF_volumes_mo*gallons_to_Liters) +\
            Commercial_SAF.loc[0]['LCEF (gCO2e/MJ)']*Commercial_SAF.loc[0]['Volumetric Energy Density (MJ/L)']*np.sum(Jet_A_fuel_volume_required_mo)*gallons_to_Liters
        
        w_o_SAF[w] = Commercial_SAF.loc[0]['LCEF (gCO2e/MJ)']*Commercial_SAF.loc[0]['Volumetric Energy Density (MJ/L)']*(np.sum(SAF_volumes_mo) + Jet_A_fuel_volume_required_mo) *gallons_to_Liters
    

           
    
    ##================================================================================================================================================  
    ## Plot Flight_Ops 
    ##================================================================================================================================================     
    ## Flight_Ops   
    #colors                = px.colors.qualitative.Pastel 
    #fig_1                = go.Figure()
    #airport_marker_size  = 5
    #airport_marker_color = "white"
 
    ## Flight Paths
    #lons       = np.empty(3 * len(Non_SAF_Flights))
    #lons[::3]  = Non_SAF_Flights['Origin Longitude (Deg.)']
    #lons[1::3] = Non_SAF_Flights['Destination Longitude (Deg.)']
    #lons[2::3] = None
    #lats       = np.empty(3 * len(Non_SAF_Flights))
    #lats[::3]  = Non_SAF_Flights['Origin Latitude (Deg.)']
    #lats[1::3] = Non_SAF_Flights['Destination Latitude (Deg.)']
    #lats[2::3] = None    
  
    #fig_1.add_trace(
        #go.Scattergeo( 
            #lon = lons,
            #lat = lats,
            #mode = 'lines',
            #opacity= 0.5,
            #line = dict(width = 1,color = colors[2]), ))
    
    #lons       = np.empty(3 * len(Flight_at_SAF_Airports_Using_SAF))
    #lons[::3]  = Flight_at_SAF_Airports_Using_SAF['Origin Longitude (Deg.)']
    #lons[1::3] = Flight_at_SAF_Airports_Using_SAF['Destination Longitude (Deg.)']
    #lons[2::3] = None
    #lats       = np.empty(3 * len(Flight_at_SAF_Airports_Using_SAF))
    #lats[::3]  = Flight_at_SAF_Airports_Using_SAF['Origin Latitude (Deg.)']
    #lats[1::3] = Flight_at_SAF_Airports_Using_SAF['Destination Latitude (Deg.)']
    #lats[2::3] = None    
  
    #fig_1.add_trace(
        #go.Scattergeo( 
            #lon = lons,
            #lat = lats,
            #mode = 'lines',
            #opacity= 0.5,
            #line = dict(width = 1,color = colors[4]), ))
  
    ## Airports  
    #fig_1.add_trace(go.Scattergeo( 
        #lon = Flight_Ops['Destination Longitude (Deg.)'],
        #lat = Flight_Ops['Destination Latitude (Deg.)'], 
        #text = Flight_Ops['Destination City'],
        #mode = 'markers',
        #marker = dict(
            #size = airport_marker_size,
            #color = airport_marker_color, ))) 

    #fig_1.add_trace(go.Scattergeo( 
        #lon = Flight_Ops['Origin Longitude (Deg.)'],
        #lat = Flight_Ops['Origin Latitude (Deg.)'], 
        #text = Flight_Ops['Origin City'],
        #mode = 'markers',
        #marker = dict(
            #size = airport_marker_size,
            #color = airport_marker_color,))) 
     
    ## Flight Paths 
    #fig_1.update_layout(mapbox_style  = "open-street-map",      
                      #showlegend    = False, 
                      #height        = 400, 
                      #geo_scope     ='usa',
                      #margin        = {'t':0,'l':0,'b':0,'r':0},  
                      #mapbox        = dict( accesstoken=mapbox_access_token,style=map_style,
                                            #center=go.layout.mapbox.Center( lat=30, lon= 230 )))   

    ##================================================================================================================================================      
    ## Passenger vs Distance Traveled 
    ##================================================================================================================================================     
    #fig_2               = go.Figure()  
    #fig_2.add_trace(go.Histogram(histfunc="sum",
                               #x= Flight_at_SAF_Airports_Using_SAF['Distance (miles)'],
                               #y = Flight_at_SAF_Airports_Using_SAF['Passengers'],
                               #name='SAF', 
                               #xbins=dict(start=0, end=4000, size=500),
                               #marker_color=colors[4],))
    #fig_2.add_trace(go.Histogram(histfunc="sum",
                               #x= Non_SAF_Flights['Distance (miles)'],
                               #y = Non_SAF_Flights['Passengers'],
                               #name='Jet-A',
                               #xbins=dict(start=0, end=4000, size=500),
                               #marker_color=colors[2],)) 
    
    ## The two histograms are drawn on top of another
    #fig_2.update_layout(barmode='stack', 
                      #xaxis_title_text='Distance (miles)', 
                      #yaxis_title_text='Passengers',
                      #height        = 300, 
                      #width         = 600, 
                      #margin        = {'t':0,'l':0,'b':0,'r':0},  
                      #bargap        = 0.1,
                      #font=dict(  size=font_size ))   
    
    ##================================================================================================================================================      
    ## Busiest Airports 
    ##================================================================================================================================================    
    #fig_3 = go.Figure()
    #Airport_Routes     = Flight_at_SAF_Airports_Using_SAF[['Passengers','Origin Airport','Destination City']]
    #Cumulative_Flights = Airport_Routes.groupby(['Origin Airport']).sum()
    #Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(10) 
    #Alphabetical_List  = Busiest_Airports.sort_values(by=['Origin Airport'])  
    #fig_3.add_trace(go.Bar( x=list(Alphabetical_List['Passengers'].index),
                       #y=np.array(Alphabetical_List['Passengers']),
                       #marker_color=colors[4])) 
    #fig_3.update_layout(xaxis_title_text='Airport', 
                      #yaxis_title_text='Passengers', 
                      #height        = 300, 
                      #width         = 600, 
                      #margin        = {'t':0,'l':0,'b':0,'r':0},  
                      #bargap        = 0.1,
                      #font=dict(  size=font_size ))  
    
    ##================================================================================================================================================      
    ## Determine Ratio of SAF to Jet-A Routes
    ##================================================================================================================================================    
    #fig_4                       = go.Figure() 
    #sector_colors               = [colors[4],colors[2]]
    #Feasible_Passenger_Miles    = np.sum(np.array(Flight_at_SAF_Airports_Using_SAF['Passengers'])* np.array(Flight_at_SAF_Airports_Using_SAF['Distance (miles)']))
    #Infeasible_Passenger_Miles  = np.sum(np.array(Non_SAF_Flights[['Passengers']])* np.array(Non_SAF_Flights[['Distance (miles)']]))
    #labels                      = ["SAF", "Jet-A"] 
    #fig_4.add_trace(go.Pie(labels=labels,
                         #values=[Feasible_Passenger_Miles, Infeasible_Passenger_Miles],
                         #marker_colors=sector_colors)) 
    #fig_4.update_traces(hole=.4, hoverinfo="label+percent+name") 
    #fig_4.update_layout( height     = 400, 
                      #width         = 600, 
                      #margin        = {'t':50,'l':0,'b':0,'r':0},  
                      #font=dict(  size=font_size ))  
    
    
    ##================================================================================================================================================  
    ## Plot Land Use  
    ##================================================================================================================================================     
 
    #Land_Use  = pd.concat([Used_Feedstock, non_feedstock_states] )    
    #fig_5 = px.choropleth(Land_Use, geojson=counties, locations='FIPS',color = 'Feedstock_Usage',
                           #color_continuous_scale="algae", 
                           #hover_data=["County","State","Acres Harvested"],
                           #scope='usa',
                           #range_color=(0,1)
                           #) 
    #fig_5.update(layout_coloraxis_showscale=False)
    
    ##================================================================================================================================================      
    ## Cost Per Seat Mile
    ##================================================================================================================================================   
    #fig_6 = go.Figure()       
    #month_names         = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']    
    #fig_6.add_trace(go.Scatter(x=month_names, y=CASM_SAF, name = 'Electric',
                             #line=dict(color= colors[4], width=4)))  
    #fig_6.add_trace(go.Scatter(x=month_names, y=CASM_jet_A, name='Jet-A',
                             #line=dict(color= colors[2], width=4)))  
    #fig_6.update_layout( 
                      #height           = 400, 
                      #width            = 600, 
                      #margin           = {'t':50,'l':0,'b':0,'r':0},
                      #yaxis_title_text ='Cost Per Seat Mile* (cents)', 
                      #font=dict(  size=font_size ),
                      #legend=dict(
                          #yanchor="top",
                          #y=0.99,
                          #xanchor="center",
                          #x=0.4 )) 

    ##================================================================================================================================================      
    ## Emissions Comparison 
    ##================================================================================================================================================   
    #month_names         = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']      
    #fig_8               = go.Figure() 
    #fig_8.add_trace(go.Scatter(x=month_names, y=w_SAF*1e-6, name = 'SAF Aircraft in Fleet',
                             #line=dict(color=colors[1], width=4)))  
    #fig_8.add_trace(go.Scatter(x=month_names, y=w_o_SAF*1e-6, name='No SAF Aircraft in Fleet',
                             #line=dict(color=colors[10], width=4)))   
    #fig_8.update_layout( 
                      #height           = 400, 
                      #width            = 600, 
                      #margin           = {'t':50,'l':0,'b':0,'r':0},
                      #yaxis_title_text ='CO2e (Ton)', # yaxis label
                      #font=dict(  size=font_size ),
                      #legend=dict(
                          #yanchor="top",
                          #y=0.99,
                          #xanchor="center",
                          #x=0.4 )) 
    
    ##================================================================================================================================================      
    ##  Life Cycle Analysis
    ##================================================================================================================================================   
    #Jet_A_name = ["Jet-A"]*num_fuels
    #Jet_A_data = {'Cumulative Fuel': Jet_A_name,
            #'Fuels' : selected_fuels,
            #'Cumulative LCA Value'        : Jet_A_LCA_val,
            #} 
    #SAF_name = ["SAF"]*num_fuels
    #SAF_data = {'Cumulative Fuel': SAF_name,
            #'Fuels' : selected_fuels,
            #'Cumulative LCA Value'        : SAF_LCA_val*cumulative_fuel_use,
            #} 
    #Jet_A_Emissions = pd.DataFrame(Jet_A_data)
    #SAF_Emissions   = pd.DataFrame(SAF_data)
    #frames = [Jet_A_Emissions,SAF_Emissions] 
    #Emissions = pd.concat(frames)
    #fig_7               =   px.bar(Emissions,
                                   #x="Cumulative Fuel",
                                   #y="Cumulative LCA Value",
                                   #color= "Fuels",
                                   #color_discrete_sequence=px.colors.qualitative.Pastel)  
    #fig_7.update_layout(xaxis_title=None)
    
    ##fig_1["layout"]["template"] = template  
    ##fig_2["layout"]["template"] = template 
        
    return vector, w_SAF, w_o_SAF, CASM_SAF, CASM_jet_A



#def save_current_csv(res_No_Hex,
                     #save_csv = True,
                        #show_legend = True,
                        #save_filename = "Voltage_variation",
                        #file_type = ".csv"):
    #current = []
    #t    = []
    #for network in res_No_Hex.segments[0].analyses.energy.networks: 
        #busses  = network.busses
        #for bus in busses: 
            #for battery in bus.batteries:         
    
                #for i in range(len(res_No_Hex.segments)): 
        
                    #time     = res_No_Hex.segments[i].conditions.frames.inertial.time[:,0] / Units.s
                    #battery_conditions_No_Hex  = res_No_Hex.segments[i].conditions.energy[bus.tag][battery.tag]         
                    #pack_current        = battery_conditions_No_Hex.cell.current[:,0]
                    #reservoir_temperature = battery_conditions_No_Hex.thermal_management_system.RES.coolant_temperature[:,0]
                    #voltage = battery_conditions_No_Hex.cell.voltage_under_load[:,0]
                    #t_ambient           = res_No_Hex.segments[i].conditions.freestream.temperature[:,0]
                    
                    
                    #current.append(voltage)
                    #t.append(time)
          
        #current = np.hstack(current)
        #t       = np.hstack(t)
        #df = pd.DataFrame({'Time': t, 'Voltage': current})
        
        #if save_csv:
            #df.to_csv('{save_filename}{file_type}', index=False)
        #return  
 
# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS COMMON AXIS 
# ----------------------------------------------------------------------------------------------------------------------    
## @ingroup Visualization-Performance-Common
def set_axes(axes):
    """This sets the axis parameters for all plots

    Assumptions:
    None

    Source:
    None

    Inputs
    axes

    Outputs:
    axes

    Properties Used:
    N/A 
    """

    axes.minorticks_on()
    #axes.grid(which='major', linestyle='-', linewidth=0.5, color='grey')
    #axes.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')
    axes.grid(True)
    #axes.get_yaxis().get_major_formatter().set_scientific(False)
    #axes.get_yaxis().get_major_formatter().set_useOffset(False)

    return


# ----------------------------------------------------------------------------------------------------------------------
#  PLOTS STYLE
# ----------------------------------------------------------------------------------------------------------------------   
## @ingroup Visualization-Performance-Common
def plot_style(): 
    plt.rcParams['axes.linewidth'] = 1.
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 20,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14,
                  'axes.titlesize': 18,
                  'figure.dpi': 1200
                  }


    # Universal Plot Settings  
    plt.rcParams.update(parameters)
    plot_parameters                        = Data()
    plot_parameters.line_width             = 1.5  
    plot_parameters.line_style             = ['-','-']
    plot_parameters.marker_size            = 4
    plot_parameters.legend_fontsize        = '12'
    plot_parameters.legend_title_font_size = 14
    plot_parameters.axis_font_size         = 16
    plot_parameters.title_font_size        = 16   
    plot_parameters.markers                =  ['o','x','o','v','P','p','^','D','*']
    plot_parameters.color                  = 'black'


    return plot_parameters


def generate_plot(vector, w_SAF, w_o_SAF, CASM_SAF, CASM_jet_A,
                  save_figure = True,
                  show_legend = True,
                  save_filename = "Flight_Profile_aircraft_speed",
                  file_type = ".png"):
    # get plotting style 
    ps      = plot_style()   

    # get line colors for plots 
    fig = plt.figure(save_filename)
    fig.set_size_inches(16,8)   
    axis_1 = plt.subplot(1,2,1)
    axis_2 = plt.subplot(1,2,2)

    axis_1.grid(False)
   
    axis_1.plot(vector, w_SAF, linestyle = ps.line_style[0], markersize= ps.marker_size, marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width)
    axis_1.plot(vector, w_o_SAF, linestyle = ps.line_style[1], markersize= ps.marker_size, marker = ps.markers[1],  fillstyle = 'none', linewidth = ps.line_width)
    #axis_1.plot(vector, CASM,_SAF)
    axis_1.set_xlabel('Adoption [%]',fontsize = ps.axis_font_size)
    axis_1.set_ylabel('CO2e [-]',fontsize = ps.axis_font_size)
    #axis_1.set_title('CO2e vs. Adoption',fontsize = ps.axis_font_size)    
    set_axes(axis_1)
    
   
    axis_2.plot(vector[1:-1], CASM_SAF[1:-1], linestyle = ps.line_style[0], markersize= ps.marker_size, marker = ps.markers[0],  fillstyle = 'none', linewidth = ps.line_width,label= 'SAF')
    axis_2.plot(vector[1:-1], CASM_jet_A[1:-1], linestyle = ps.line_style[1], markersize= ps.marker_size, marker = ps.markers[1],  fillstyle = 'none', linewidth = ps.line_width,label= 'Jet A')
    axis_2.set_xlabel('Adoption [%]',fontsize = ps.axis_font_size)
    axis_2.set_ylabel('CASM [-]',fontsize = ps.axis_font_size)
    #axis_2.set_title('CO2e vs. Adoption',fontsize = ps.axis_font_size)    
    set_axes(axis_2)    
    
    
    if show_legend:        
        #leg =  axis_1.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.70, 0.9),loc='lower left')   
        leg =  fig.legend(fontsize = ps.legend_fontsize,bbox_to_anchor=(0.9, 0.85),loc='lower left')   
    
    #fig_1 = plt.figure()
 
    fig.tight_layout()

    if save_figure:
        plt.savefig(f'Plots_for_paper/{save_filename}{file_type}')    
    plt.show()
    return fig
 
 
 
 
if __name__ == '__main__':
    main()