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
import json
 
 
def main():
    
    # ---------------------------------------------------------------------------------------------------------------------------------------------------
    # Data
    # ---------------------------------------------------------------------------------------------------------------------------------------------------
    
   
    separator                       = os.path.sep
    H2_filename                     = 'Data' + separator + 'Technology' + separator +  'Technology_Data.xlsx'
    routes_filename                 = 'Data' + separator +  'Air_Travel' + separator + 'American_Airlines_Flight_Ops_and_Climate.xlsx' 
    
    SAT_data                       = pd.read_excel(H2_filename,sheet_name=['Commercial_H2']) 
    Commercial_H2                  = SAT_data['Commercial_H2']   
    a                              = Commercial_H2['Process'] 
    b                              = Commercial_H2['Feedstock']
    Commercial_H2["H2 Fuel Name"]  =  a  + ' from ' + b 
    Flight_Ops                     = pd.read_excel(routes_filename,sheet_name=['Sheet1']) 
    Flight_Ops                     = Flight_Ops['Sheet1']
    
    aircraft              = 'Boeing 737 MAX-8' 
    volume_fraction       = 35   
    selected_airpots      = "All Airports" # ["Top 5 Airports","Top 10 Airports","Top 20 Airports", "Top 50 Airports",  "All Airports"]
    percent_adoption      = 100 
    month_no              = 1 
    switch_off            = False  
    H2_dollars_per_gal    =  0 
    fig_5, fig_6, fig_7, fig_8,fig_9,fig_10 = generate_electric_flight_operations_plots(Flight_Ops,Commercial_H2,aircraft,selected_airpots,percent_H2_process , volume_fraction,percent_adoption,month_no,H2_dollars_per_gal ,switch_off)
                
      
    return 
    
     
def generate_electric_flight_operations_plots(Flight_Ops,Commercial_H2,aircraft,selected_airpots,percent_H2_process , volume_fraction,percent_adoption,month_no,H2_dollars_per_gal ,switch_off): 
    mapbox_access_token  = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ"     
    map_style            = None if switch_off else 'dark'  
    #template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"] 
    font_size            = 16     
    
    if aircraft == 'ATR 72-600': 
        P_max           = 1568083 
        W_0             = 23000   
        thrust_coef     = 0.000017678 # incorrect  
        fuel_volume     = 6410.2
        fuel_economy    = 3.27 
        passengers      = 72 
        CL_cruise       = 0.55   # incorrect  
        fuselage_volume = 100 # CHANGE 
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
        #fuselage_volume = 100 # CHANGE 
        #SFC             = 100 # CHANGE  
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
        #SFC             = 100 # CHANGE  
        CL_cruise       = 0.55
        fuselage_volume = 100 # CHANGE 
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
        SFC             = 100 # CHANGE  
        fuselage_volume = 100 # CHANGE 
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
        #SFC             = 100 # CHANGE  
        #W_0             = 227900   
        #CL_cruise       = 0.055 
        #fuselage_volume = 100 # CHANGE 
        #CD_cruise       = 0.035
        #S_ref           = 377
        #L_div_D         = CL_cruise/CD_cruise
        
    #elif aircraft == 'Airbus A350-1000':
        #P_max          =  
        #thrust_coef    = 0.000017678   
        #fuel_volume    = 166488        # liters 
        #fuel_economy   = 2.85          # L/100km per Pax
        #SFC             = 100 # CHANGE  
        #passengers     = 327  
        #W_0            = 322050 
        #fuselage_volume = 100 # CHANGE  
        #CL_cruise      = 0.055 
        #CD_cruise      = 0.035
        #S_ref          = 464.3
        #L_div_D        = CL_cruise/CD_cruise      
              
    
    density_JetA  = 820.0  
    density_H2    = 70 # kg /m3  
    gallons_to_Liters  = 3.78541
    airspeed   =  0.75 * 343
    SFC_H2    =  SFC * 0.35 

    
    H2_vol   =  volume_fraction * fuselage_volume
    W_final  =  W_0  -  (H2_vol * density_H2) *0.9 # only 90% of fuel is used up 
    
    #================================================================================================================================================  
    # Compute Feasible Routes 
    #================================================================================================================================================        
    # Step 1: Compute Range   
    H2_volume       = fuselage_volume * volume_fraction   
    Range           = (airspeed / 9.81) * (1 / SFC_H2) *  (L_div_D) *  (W_0 / W_final)
    Range_mi        = Range * 0.000621371 
    
    #  Step 2: Compute modified passenger capacity  
    Routes_and_Temp_Mo                   = Flight_Ops[Flight_Ops['Month'] == month_no+1 ]    
    passenger_reductions                 = H2_volume * (fuselage_volume /passengers)
    original_Pax_volume                  = np.array(Routes_and_Temp_Mo['Passengers'])
    remaining_Pax                        = original_Pax_volume - passenger_reductions*np.array(Routes_and_Temp_Mo['No of Flights Per Month'])
    remaining_Pax[remaining_Pax<0]       = 0
    Routes_and_Temp_Mo['H2_Passengers']  = remaining_Pax
    
       
    # Step 1: compute the percentages of different types of fuels used 
    H2_process_percentages_list  = [0]
    H2_process_percentages_list  += percent_H2_process
    H2_process_percentages_list  += [100]
    H2_process_percentages       = np.diff(np.array(H2_process_percentages_list))/100  
    
    
    # Step 2: determine the percentage of different H2 production processes 
    if Commercial_H2['H2 Fuel Name'][0] not in selected_fuels:  
        selected_fuels         = [Commercial_H2['H2 Fuel Name'][0]] + selected_fuels
        H2_process_percentages = np.hstack((np.array([0]),H2_process_percentages)) 
    mask                  = Commercial_H2['H2 Fuel Name'].isin(selected_fuels)
    fuels_used            = Commercial_H2[mask] 
    
    num_fuels             = len(selected_fuels)
    cumulative_fuel_use   = np.zeros(num_fuels) 
    H2_LCA_val            = np.zeros(num_fuels) 
    Jet_A_LCA_val         = np.zeros(num_fuels)   
    for i in range(1,num_fuels):
        # loop through fuels and get percentage of fuel used by each type
        cumulative_fuel_use[i]  = H2_process_percentages[i]
        H2_LCA_val[i]           = Commercial_H2[Commercial_H2['H2 Fuel Name'] == selected_fuels[i]]['LCA Value'] 
    H2_LCA_val[0]     = 89
    Jet_A_LCA_val[0]  = 89 
      
    # Step 3: Filter flight data based on option selected: i.e. top 10, top 20, top 50, all airpots 
    Airport_Routes     = Flight_Ops[['Passengers','Origin Airport','Destination City']]
    Cumulative_Flights = Airport_Routes.groupby('Origin Airport', as_index=False).sum()[Airport_Routes.columns]
    if  selected_airpots   == "Top 5 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(5) 
    if  selected_airpots   == "Top 10 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(10) 
    elif  selected_airpots =="Top 20 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(20) 
    elif  selected_airpots == "Top 50 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(50) 
    elif  selected_airpots == "All Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False) 
    Airport_List = list(Busiest_Airports['Origin Airport'])
    
    # Step 4: Filter airports that will support H2 and those that wont support H2 
    mask_1               = Flight_Ops['Origin Airport'].isin(Airport_List)
    H2_Airports          = Flight_Ops[mask_1]
    Non_H2_Airports      = Flight_Ops[~mask_1] 
    Feasible_Routes_0    = H2_Airports[H2_Airports['H2_Passengers'] > 0 ] 
    Infeasible_Routes_0  = H2_Airports[H2_Airports['H2_Passengers'] < 0 ]   
    Feasible_Routes_1    = Feasible_Routes_0[Feasible_Routes_0['Distance (miles)'] < Range_mi ] 
    Infeasible_Routes_1  = Feasible_Routes_0[Feasible_Routes_0['Distance (miles)'] > Range_mi ] 
    
    # Step 8: Out of H2 supporting airports, use the percent adoption to determine how many flights at that airport will use H2     
    Flight_at_H2_Airports_Using_H2  = Feasible_Routes_1.sample(frac=(percent_adoption/100))
    Infeasible_Routes_2             = Feasible_Routes_1[~Feasible_Routes_1.index.isin(Flight_at_H2_Airports_Using_H2.index)]
    Non_H2_Flights                  = pd.concat([Non_H2_Airports,Infeasible_Routes_0,Infeasible_Routes_1,Infeasible_Routes_2])# add list of flights from non supporting airports to non-H2 flights  
        
    # Step 9: Get total volume of each H2 required at the airports 
    equivalent_JetA_consumed = np.array(Flight_at_H2_Airports_Using_H2['Total Fuel Per Route (Gal)'])
    airport_H2_volumes       = (density_JetA / density_H2) *  equivalent_JetA_consumed  
        
    # Steps 14 and 15 Determine Cost per Seat Mile and Emissions  
    CASM_jet_A  = np.zeros(12) 
    CASM_H2     = np.zeros(12) 
    Emissions_w_H2       = np.zeros(12) 
    Emissions_w_o_H2     = np.zeros(12) 
    for m_i in range(12): 
        Flight_at_H2_Airports_Using_H2_Mo  = Flight_at_H2_Airports_Using_H2.loc[Flight_at_H2_Airports_Using_H2['Month'] == m_i+1 ]    
        H2_volumes_mo                      = np.sum(np.array(Flight_at_H2_Airports_Using_H2_Mo['Total Fuel Per Route (Gal)']))
        
        Non_H2_Flights_Mo                  = Non_H2_Flights.loc[Non_H2_Flights['Month'] == m_i+1 ]   
        Jet_A_fuel_volume_required_mo      = np.sum(np.array(Non_H2_Flights_Mo['Total Fuel Per Route (Gal)']))    
        
        
        if len(Non_H2_Flights_Mo ) == 0:
            CASM_jet_A[m_i] = 0
        else:
            ASM_jet_A               = np.sum(Non_H2_Flights_Mo['Distance (miles)'] * Non_H2_Flights_Mo['Passengers'])
            Total_Fuel_Cost_jet_A   = np.sum(Non_H2_Flights_Mo['Fuel Cost']) 
            CASM_jet_A[m_i]         = 100*Total_Fuel_Cost_jet_A/ASM_jet_A     
       
        if len(Flight_at_H2_Airports_Using_H2_Mo)  == 0:
            CASM_H2[m_i]   
        else:    
            ASM_H2               = np.sum(Flight_at_H2_Airports_Using_H2_Mo['Distance (miles)'] * Flight_at_H2_Airports_Using_H2_Mo['Passengers']) 
            Total_Fuel_Cost_H2   = np.sum(Flight_at_H2_Airports_Using_H2_Mo['Total Fuel Per Route (Gal)'] ) * H2_dollars_per_gal 
            CASM_H2[m_i]         = 100*Total_Fuel_Cost_H2/ASM_H2   
             
        Emissions_w_H2[m_i]   = np.sum(fuels_used['LCEF (gCO2e/MJ)']*fuels_used['Volumetric Energy Density (MJ/L)']*H2_volumes_mo*gallons_to_Liters) +\
                                Commercial_H2.loc[0]['LCEF (gCO2e/MJ)']*Commercial_H2.loc[0]['Volumetric Energy Density (MJ/L)']*np.sum(Jet_A_fuel_volume_required_mo)*gallons_to_Liters 
        Emissions_w_o_H2[m_i] = Commercial_H2.loc[0]['LCEF (gCO2e/MJ)']*Commercial_H2.loc[0]['Volumetric Energy Density (MJ/L)']*(np.sum(H2_volumes_mo) + Jet_A_fuel_volume_required_mo) *gallons_to_Liters
           
    
    #================================================================================================================================================  
    # Plot Flight_Ops 
    #================================================================================================================================================     
    # Flight_Ops   
    colors               = px.colors.qualitative.Pastel 
    fig_1                = go.Figure()
    airport_marker_size  = 5
    airport_marker_color = "white"
 
    # Flight Paths
    lons       = np.empty(3 * len(Non_H2_Flights))
    lons[::3]  = Non_H2_Flights['Origin Longitude (Deg.)']
    lons[1::3] = Non_H2_Flights['Destination Longitude (Deg.)']
    lons[2::3] = None
    lats       = np.empty(3 * len(Non_H2_Flights))
    lats[::3]  = Non_H2_Flights['Origin Latitude (Deg.)']
    lats[1::3] = Non_H2_Flights['Destination Latitude (Deg.)']
    lats[2::3] = None    
  
    fig_1.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            opacity= 0.5,
            line = dict(width = 1,color = colors[2]), ))
    
    lons       = np.empty(3 * len(Flight_at_H2_Airports_Using_H2))
    lons[::3]  = Flight_at_H2_Airports_Using_H2['Origin Longitude (Deg.)']
    lons[1::3] = Flight_at_H2_Airports_Using_H2['Destination Longitude (Deg.)']
    lons[2::3] = None
    lats       = np.empty(3 * len(Flight_at_H2_Airports_Using_H2))
    lats[::3]  = Flight_at_H2_Airports_Using_H2['Origin Latitude (Deg.)']
    lats[1::3] = Flight_at_H2_Airports_Using_H2['Destination Latitude (Deg.)']
    lats[2::3] = None    
  
    fig_1.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            opacity= 0.5,
            line = dict(width = 1,color = colors[4]), ))
  
    # Airports  
    fig_1.add_trace(go.Scattergeo( 
        lon = Flight_Ops['Destination Longitude (Deg.)'],
        lat = Flight_Ops['Destination Latitude (Deg.)'], 
        text = Flight_Ops['Destination City'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color, ))) 

    fig_1.add_trace(go.Scattergeo( 
        lon = Flight_Ops['Origin Longitude (Deg.)'],
        lat = Flight_Ops['Origin Latitude (Deg.)'], 
        text = Flight_Ops['Origin City'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color,))) 
     
    # Flight Paths 
    fig_1.update_layout(mapbox_style  = "open-street-map",      
                      showlegend    = False, 
                      height        = 400, 
                      geo_scope     ='usa',
                      margin        = {'t':0,'l':0,'b':0,'r':0},  
                      mapbox        = dict( accesstoken=mapbox_access_token,style=map_style,
                                            center=go.layout.mapbox.Center( lat=30, lon= 230 )))   

    #================================================================================================================================================      
    # Passenger vs Distance Traveled 
    #================================================================================================================================================     
    fig_2               = go.Figure()  
    fig_2.add_trace(go.Histogram(histfunc="sum",
                               x= Flight_at_H2_Airports_Using_H2['Distance (miles)'],
                               y = Flight_at_H2_Airports_Using_H2['Passengers'],
                               name='H2', 
                               xbins=dict(start=0, end=4000, size=500),
                               marker_color=colors[4],))
    fig_2.add_trace(go.Histogram(histfunc="sum",
                               x= Non_H2_Flights['Distance (miles)'],
                               y = Non_H2_Flights['Passengers'],
                               name='Jet-A',
                               xbins=dict(start=0, end=4000, size=500),
                               marker_color=colors[2],)) 
    
    # The two histograms are drawn on top of another
    fig_2.update_layout(barmode='stack', 
                      xaxis_title_text='Distance (miles)', 
                      yaxis_title_text='Passengers',
                      height        = 300, 
                      width         = 600, 
                      margin        = {'t':0,'l':0,'b':0,'r':0},  
                      bargap        = 0.1,
                      font=dict(  size=font_size ))   
    
    #================================================================================================================================================      
    # Busiest Airports 
    #================================================================================================================================================    
    fig_3 = go.Figure()
    Airport_Routes     = Flight_at_H2_Airports_Using_H2[['Passengers','Origin Airport','Destination City']]
    Cumulative_Flights = Airport_Routes.groupby(['Origin Airport']).sum()
    Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(10) 
    Alphabetical_List  = Busiest_Airports.sort_values(by=['Origin Airport'])  
    fig_3.add_trace(go.Bar( x=list(Alphabetical_List['Passengers'].index),
                       y=np.array(Alphabetical_List['Passengers']),
                       marker_color=colors[4])) 
    fig_3.update_layout(xaxis_title_text='Airport', 
                      yaxis_title_text='Passengers', 
                      height        = 300, 
                      width         = 600, 
                      margin        = {'t':0,'l':0,'b':0,'r':0},  
                      bargap        = 0.1,
                      font=dict(  size=font_size ))  
    
    #================================================================================================================================================      
    # Determine Ratio of H2 to Jet-A Routes
    #================================================================================================================================================    
    fig_4                       = go.Figure() 
    sector_colors               = [colors[4],colors[2]]
    Feasible_Passenger_Miles    = np.sum(np.array(Flight_at_H2_Airports_Using_H2['Passengers'])* np.array(Flight_at_H2_Airports_Using_H2['Distance (miles)']))
    Infeasible_Passenger_Miles  = np.sum(np.array(Non_H2_Flights[['Passengers']])* np.array(Non_H2_Flights[['Distance (miles)']]))
    labels                      = ["H2", "Jet-A"] 
    fig_4.add_trace(go.Pie(labels=labels,
                         values=[Feasible_Passenger_Miles, Infeasible_Passenger_Miles],
                         marker_colors=sector_colors)) 
    fig_4.update_traces(hole=.4, hoverinfo="label+percent+name") 
    fig_4.update_layout( height     = 400, 
                      width         = 600, 
                      margin        = {'t':50,'l':0,'b':0,'r':0},  
                      font=dict(  size=font_size ))  
     
    #================================================================================================================================================      
    # Cost Per Seat Mile
    #================================================================================================================================================   
    fig_6 = go.Figure()       
    month_names         = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']    
    fig_6.add_trace(go.Scatter(x=month_names, y=CASM_H2, name = 'Electric',
                             line=dict(color= colors[4], width=4)))  
    fig_6.add_trace(go.Scatter(x=month_names, y=CASM_jet_A, name='Jet-A',
                             line=dict(color= colors[2], width=4)))  
    fig_6.update_layout( 
                      height           = 400, 
                      width            = 600, 
                      margin           = {'t':50,'l':0,'b':0,'r':0},
                      yaxis_title_text ='Cost Per Seat Mile* (cents)', 
                      font=dict(  size=font_size ),
                      legend=dict(
                          yanchor="top",
                          y=0.99,
                          xanchor="center",
                          x=0.4 )) 

    #================================================================================================================================================      
    # Emissions Comparison 
    #================================================================================================================================================   
    month_names         = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']      
    fig_8               = go.Figure() 
    fig_8.add_trace(go.Scatter(x=month_names, y=Emissions_w_H2*1e-6, name = 'H2 Aircraft in Fleet',
                             line=dict(color=colors[1], width=4)))  
    fig_8.add_trace(go.Scatter(x=month_names, y=Emissions_w_o_H2*1e-6, name='No H2 Aircraft in Fleet',
                             line=dict(color=colors[10], width=4)))   
    fig_8.update_layout( 
                      height           = 400, 
                      width            = 600, 
                      margin           = {'t':50,'l':0,'b':0,'r':0},
                      yaxis_title_text ='CO2e (Ton)', # yaxis label
                      font=dict(  size=font_size ),
                      legend=dict(
                          yanchor="top",
                          y=0.99,
                          xanchor="center",
                          x=0.4 )) 
    
    #================================================================================================================================================      
    #  Life Cycle Analysis
    #================================================================================================================================================   
    Jet_A_name = ["Jet-A"]*num_fuels
    Jet_A_data = {'Cumulative Fuel': Jet_A_name,
            'Fuels' : selected_fuels,
            'Cumulative LCA Value'        : Jet_A_LCA_val,
            } 
    H2_name = ["H2"]*num_fuels
    H2_data = {'Cumulative Fuel': H2_name,
            'Fuels' : selected_fuels,
            'Cumulative LCA Value'        : H2_LCA_val*cumulative_fuel_use,
            } 
    Jet_A_Emissions = pd.DataFrame(Jet_A_data)
    H2_Emissions   = pd.DataFrame(H2_data)
    frames = [Jet_A_Emissions,H2_Emissions] 
    Emissions = pd.concat(frames)
    fig_7               =   px.bar(Emissions,
                                   x="Cumulative Fuel",
                                   y="Cumulative LCA Value",
                                   color= "Fuels",
                                   color_discrete_sequence=px.colors.qualitative.Pastel)  
    fig_7.update_layout(xaxis_title=None)
    
    #fig_1["layout"]["template"] = template  
    #fig_2["layout"]["template"] = template
 
    
    return fig_4,fig_5,fig_6,fig_7,fig_8,fig_9


# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Motor Plots
# ---------------------------------------------------------------------------------------------------------------------------------------------------
def generate_motor_scatter_plot(Electric_Motor_Development,selected_x_axis,selected_y_axis,switch_off): 
    #template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"]  
    unique_manufacturers = list(Electric_Motor_Development['Manufacturer'][1:].unique())
    unique_motor_types   = list(Electric_Motor_Development['Motor Type'][1:].unique())
    marker_size          = 15
    opacity_ratio        = 0.8 if switch_off else 1.0
    font_size            = 16 

    # Brand Colors: greenyellow, aquamarine, paleturquoise, lightcoral, yellow, lavender ,thistle ,orangered   
    Brand_Colors      = px.colors.qualitative.Pastel 
    Motor_Type_Markers = ['square','x','circle','cross','diamond','triangle-up','triangle-down','star','hourglass'] 
    fig = go.Figure()  
    for i in range(len(unique_manufacturers)):  
        for j in range(len(unique_motor_types)):
            data_1    = Electric_Motor_Development.loc[Electric_Motor_Development['Manufacturer'] == unique_manufacturers[i]] 
            data_2    = data_1.loc[Electric_Motor_Development['Motor Type'] == unique_motor_types[j]]  
            developer = data_2["Manufacturer"]                 
            coolant   = data_2["Coolant"]                 
            fig.add_trace(go.Scatter( x        = np.array(data_2[selected_x_axis]), 
                                 y             = np.array(data_2[selected_y_axis]),  
                                 mode          = 'markers', 
                                 name          ="",                                    
                                 marker        = dict(size=marker_size,color=Brand_Colors[i],opacity=opacity_ratio,symbol = Motor_Type_Markers[j]),
                                 hovertemplate = 'Manufacturer: ' + developer  + '<br>' + 'Motor Type: ' + unique_motor_types[j] + '<br>' + 'Coolant: ' + coolant   + '<br>'  + selected_x_axis + ': %{x} <br>' + selected_y_axis + ': %{y}',
                                 ))  
 
    fig.update_layout(xaxis_title = selected_x_axis,
                       yaxis_title = selected_y_axis,
                       showlegend  = False, 
                       height      = 400,
                       margin      ={'t':0,'l':0,'b':0,'r':0},
                       font=dict(  size=font_size ),  
                       )     
    #fig["layout"]["template"] = template 
    return fig 


     
if __name__ == '__main__':
    main()
