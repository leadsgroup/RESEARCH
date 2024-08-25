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
    
    separator            = os.path.sep
    airline              =  'Top_10_Major_US_Airlines'  
    technology_filename  = 'Data' + separator + 'Technology' + separator + 'Technology_Data.xlsx'
    SAT_data             = pd.read_excel(technology_filename,sheet_name=['Commercial_Batteries','Battery_Development','Electric_Motor_Development']) 
    Commercial_Batteries = SAT_data['Commercial_Batteries'] 
    Electric_Motor_Development = SAT_data['Electric_Motor_Development']
    a                    = Commercial_Batteries['Brand']  
    b                    = Commercial_Batteries['Chemistry']
    c                    = Commercial_Batteries['Model']
    Commercial_Batteries["Battery Name"] = a + ': ' + b + '-' + c  
    Battery_Development  = SAT_data['Battery_Development']   
    route_temp_filename  = 'Data' + separator +  'Air_Travel' + separator + airline +  '_Flight_Ops_and_Climate.csv'  
    Routes_and_Temp      = pd.read_csv(route_temp_filename) 
    filename             = 'Data' + separator  +  'US_Climate' + separator +  'Monthly_US_County_Temperature_2019.csv'
    US_Temperature_F     =  pd.read_csv(filename)   
 
  
    cell_V                = 4
    capacity              = 4
    cell_C_max            = 5
    cell_e0               = 1200
    cell_Temp             = [10, 80]
    weight_fraction       = 30
    system_voltage_in_kV  = 1
    propulsive_efficiency = 95
    percent_adoption      = 100 
    month_no              = 1 
    switch_off            = False 
    cost_of_electricity   = 0.15 # $ per Kwh
    fig_5, fig_6, fig_7, fig_8,fig_9,fig_10 = generate_EX_aircraft_flight_ops(Routes_and_Temp,system_voltage_in_kV,weight_fraction,propulsive_efficiency,cell_V,
                                              capacity,cell_C_max,cell_e0,cell_Temp,percent_adoption,month_no,cost_of_electricity,switch_off)
                
    fig_5.show()
    fig_6.show()
    fig_7.show()
    fig_8.show()
    fig_9.show()
    fig_10.show()
  
    return  
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Flight Operations Plots 
# ---------------------------------------------------------------------------------------------------------------------------------------------------
def generate_EX_aircraft_flight_ops(Routes_and_Temp,system_voltage,weight_fraction,propulsive_efficiency,cell_V,
                                              capacity,cell_C_max,cell_e0,cell_Temp,percent_adoption,month_no,cost_of_electricity,switch_off):  
    mapbox_access_token  = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ"     
    map_style            = None if switch_off else 'dark'  
    #template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"] 
    font_size            = 16
    
    
    """
    Units in metric system if not stated otherwise 
    P_max          : maxiumum power          [Watts]
    TOGW           : takeoff gross weight    [kg]
    fuel_volume    : fuel volume             [L]
    passengers     : passenger capacity      [-]
    fuel_economy   : fuel per seat           [L/100 km]
    CL_cruise      : cruise lift coefficient [-]
    CD_cruise      : cruise drag coefficient [-]
    L_div_D        : lift to drag ratio      [-] 
    """
        
    #================================================================================================================================================  
    # Unit Conversions 
    #================================================================================================================================================    
   
    CO2e_per_mile          = 9.0736                 
    Wh_per_kg_to_J         = 3600.0
    Ah_to_C                = 3600.0
    meters_to_miles        = 0.000621371  
    V_bat                  = system_voltage*1000
    joule_to_kWh           = 2.77778e-7
    JetA_GHG               = 4.36466 # CO2e/kg fuel 
    gallons_to_Liters      = 3.78541
    liters_to_cubic_meters = 0.001
    Jet_A_density          = 800.0  
    density_JetA           = 820.0  # kg /m3  
    kg_to_Megaton          = 1E-9
    eta_0                  = propulsive_efficiency/100 
    
    #================================================================================================================================================  
    # Compute Feasible Routes 
    #================================================================================================================================================    
    # Compute Range  
    months              = ['January', 'February', 'March', 'April', 'May', 'June', 'July','August', 'September', 'October', 'November', 'December']      
    month               = months[month_no] 
    Routes_and_Temp_Mo  = Routes_and_Temp[Routes_and_Temp['Month'] == month_no+1 ]    
    V_cell              = cell_V
    e_cell              = cell_e0 *Wh_per_kg_to_J
    q_cell              = capacity * Ah_to_C  # conversion to coulombs
    i_max               = (capacity*cell_C_max) # amps   
    Min_Temp            = cell_Temp[0]
    Max_Temp            = cell_Temp[1]  
    
    # Determine Max power of aircraft fleet
    original_Pax_volume  = np.array(Routes_and_Temp_Mo['Passengers'])
    est_pax_capacity     = np.array(Routes_and_Temp_Mo['Estimated Aircraft Capacity'])

    Max_Power                                 = np.zeros_like(original_Pax_volume)
    Lift_to_Drag_ratio                        = np.zeros_like(original_Pax_volume)
    Takeoff_Weight                            = np.zeros_like(original_Pax_volume)
    Max_Power[est_pax_capacity==19]           = 1158000
    Max_Power[est_pax_capacity==88]           = 2050000  
    Max_Power[est_pax_capacity==120]          = 13567500
    Max_Power[est_pax_capacity==189]          = 15270000
    Max_Power[est_pax_capacity==368]          = 68000000
    Lift_to_Drag_ratio[est_pax_capacity==19]  = 14 
    Lift_to_Drag_ratio[est_pax_capacity==88]  = 15
    Lift_to_Drag_ratio[est_pax_capacity==120] = 16 
    Lift_to_Drag_ratio[est_pax_capacity==189] = 18
    Lift_to_Drag_ratio[est_pax_capacity==368] = 19
    Takeoff_Weight[est_pax_capacity==19]      = 6575     
    Takeoff_Weight[est_pax_capacity==88]      = 23000  
    Takeoff_Weight[est_pax_capacity==120]     = 63100      
    Takeoff_Weight[est_pax_capacity==189]     = 79015.8     
    Takeoff_Weight[est_pax_capacity==368]     = 227900
    
    # Determine limiting condition on cell configuration 
    I_bat           = Max_Power/V_bat
    n_series        = V_bat/V_cell
    W_bat           = (weight_fraction/100) * Takeoff_Weight
    E_bat           = W_bat * e_cell  
    Q_bat           = E_bat/V_bat
    n_parallel      = Q_bat/q_cell 
    n_parallel_min  = I_bat/i_max 
     
    # Compute distances between departure and destimation points  
    fuel_volume_L                        = Routes_and_Temp_Mo['Fuel Consumed Per Flight (Liters)']  
    fuel_volume_m_3                      = fuel_volume_L*0.001  
    W_f                                  = Jet_A_density*fuel_volume_m_3
    W_residual                           = W_bat-W_f  
    weight_per_pass                      = 158.757 # in kg  (250 lb for person, 100 lb for luggage) 
    passenger_reductions                 = np.ceil(np.array(W_residual)/weight_per_pass)   
    passenger_reductions[passenger_reductions < 0] = 0
    remaining_Pax                        = original_Pax_volume - passenger_reductions*np.array(Routes_and_Temp_Mo['No of Flights Per Month'])
    remaining_Pax[remaining_Pax<0]       = 0
    Routes_and_Temp_Mo['E_Passengers']   = remaining_Pax 
    
    # compute range 
    Range_mi = meters_to_miles * (e_cell/9.81) * Lift_to_Drag_ratio * (weight_fraction/100)* eta_0
    Range_mi[n_parallel_min  > n_parallel] = 0 
    
    # filter for flights that do not meet range between two airports 
    Feasible_Routes_1    = Routes_and_Temp_Mo[Routes_and_Temp_Mo['Distance (miles)'] < Range_mi ] 
    Infeasible_Routes_1  = Routes_and_Temp_Mo[Routes_and_Temp_Mo['Distance (miles)'] > Range_mi ]
    
    # filter for flights based on airline  
    Feasible_Routes_2    = Feasible_Routes_1  
    Infeasible_Routes_2  = Feasible_Routes_1.head(0)  

    # filter for flights based on aircraft type 
    Feasible_Routes_3    = Feasible_Routes_2  
    Infeasible_Routes_3  = Feasible_Routes_2.head(0)  
    
    # filter for flights that have no passengers due to too much battery weight 
    Feasible_Routes_4    = Feasible_Routes_3[Feasible_Routes_3['E_Passengers'] > 0 ] 
    Infeasible_Routes_4  = Feasible_Routes_3[Feasible_Routes_3['E_Passengers'] < 0 ]
    
    # filter for flights where batteries cannot operate 
    Feasible_Routes_5    = Feasible_Routes_4[Feasible_Routes_4['Origin ' + month] > Min_Temp] 
    Infeasible_Routes_5  = Feasible_Routes_4[Feasible_Routes_4['Origin ' + month] < Min_Temp] 
    Feasible_Routes_6    = Feasible_Routes_5[Feasible_Routes_5['Origin ' + month] < Max_Temp] 
    Infeasible_Routes_6  = Feasible_Routes_5[Feasible_Routes_5['Origin ' + month] > Max_Temp]  
    Feasible_Routes_7    = Feasible_Routes_6[Feasible_Routes_6['Destination ' + month] > Min_Temp] 
    Infeasible_Routes_7  = Feasible_Routes_6[Feasible_Routes_6['Destination ' + month] < Min_Temp] 
    Feasible_Routes_8    = Feasible_Routes_7[Feasible_Routes_7['Destination ' + month] < Max_Temp] 
    Infeasible_Routes_8  = Feasible_Routes_7[Feasible_Routes_7['Destination ' + month] > Max_Temp] 
    Feasible_Routes      = Feasible_Routes_8.head(int(len(Feasible_Routes_8)*percent_adoption/100 )) 
    Infeasible_Routes_9  = Feasible_Routes_8.tail(int(len(Feasible_Routes_8)*(100 - percent_adoption)/100 ))
    Infeasible_Routes    = pd.concat([Infeasible_Routes_1,Infeasible_Routes_2,Infeasible_Routes_3,
                                      Infeasible_Routes_4,Infeasible_Routes_5,Infeasible_Routes_6,
                                      Infeasible_Routes_7,Infeasible_Routes_8,Infeasible_Routes_9])    
    #================================================================================================================================================      
    # Monthly Analysis 
    #================================================================================================================================================     
    
    Emissions_wo_Electric_Aircraft    = np.zeros(12)
    Emissions_w_Electric_Aircraft     = np.zeros(12) 
    CASM_wo_E_Aircraft                = np.zeros(12) 
    CASM_w_E_Aircraft                 = np.zeros(12) 

    for m_i in range(12):  
        Routes_and_Temp_Mo  = Routes_and_Temp.loc[Routes_and_Temp['Month'] == m_i+1 ] 
        V_cell                  = cell_V
        e_cell                  = cell_e0 *Wh_per_kg_to_J
        q_cell                  = capacity * Ah_to_C  # conversion to coulombs
        i_max                   = (capacity*cell_C_max) # amps   
        Min_Temp                = cell_Temp[0]
        Max_Temp                = cell_Temp[1]
        
        # Determine Max power of aircraft fleet
        original_Pax_volume                       = np.array(Routes_and_Temp_Mo['Passengers'])    
        est_pax_capacity                          = np.array(Routes_and_Temp_Mo['Estimated Aircraft Capacity'])
        Max_Power                                 = np.zeros_like(original_Pax_volume)
        Lift_to_Drag_ratio                        = np.zeros_like(original_Pax_volume)
        Takeoff_Weight                            = np.zeros_like(original_Pax_volume)
        Max_Power[est_pax_capacity==19]           = 1158000
        Max_Power[est_pax_capacity==88]           = 2050000  
        Max_Power[est_pax_capacity==120]          = 13567500
        Max_Power[est_pax_capacity==189]          = 15270000
        Max_Power[est_pax_capacity==368]          = 68000000
        Lift_to_Drag_ratio[est_pax_capacity==19]  = 14 
        Lift_to_Drag_ratio[est_pax_capacity==88]  = 15
        Lift_to_Drag_ratio[est_pax_capacity==120] = 16 
        Lift_to_Drag_ratio[est_pax_capacity==189] = 18
        Lift_to_Drag_ratio[est_pax_capacity==368] = 19
        Takeoff_Weight[est_pax_capacity==19]      = 6575     
        Takeoff_Weight[est_pax_capacity==88]      = 23000  
        Takeoff_Weight[est_pax_capacity==120]     = 63100      
        Takeoff_Weight[est_pax_capacity==189]     = 79015.8     
        Takeoff_Weight[est_pax_capacity==368]     = 227900   
        
        # Determine limiting condition on cell configuration 
        I_bat           = Max_Power/V_bat
        n_series        = V_bat/V_cell
        W_bat           = (weight_fraction/100) *Takeoff_Weight
        E_bat           = W_bat * e_cell  
        Q_bat           = E_bat/V_bat
        n_parallel      = Q_bat/q_cell 
        n_parallel_min  = I_bat/i_max 
         
        # Compute distances between departure and destimation points  
        fuel_volume_L                                  = Routes_and_Temp_Mo['Fuel Consumed Per Flight (Liters)']  
        fuel_volume_m_3                                = fuel_volume_L*0.001  
        W_f                                            = Jet_A_density*fuel_volume_m_3
        W_residual                                     = W_bat-W_f  
        weight_per_pass                                = 158.757 # in kg  (250 lb for person, 100 lb for luggage) 
        passenger_reductions                           = np.ceil(np.array(W_residual)/weight_per_pass)   
        passenger_reductions[passenger_reductions < 0] = 0
        remaining_Pax                                  = original_Pax_volume - passenger_reductions*np.array(Routes_and_Temp_Mo['No of Flights Per Month'])
        remaining_Pax[remaining_Pax<0]                 = 0
        Routes_and_Temp_Mo['E_Passengers']             = remaining_Pax 
        Routes_and_Temp_Mo['Aircraft_Battery_Energy']  = E_bat
        
        # compute range 
        Range_mi = meters_to_miles * (e_cell/9.81) * Lift_to_Drag_ratio * (weight_fraction/100)* eta_0
        Range_mi[n_parallel_min  > n_parallel] = 0  
        
        # filter for flights that do not meet range between two airports 
        Feasible_Routes_1    = Routes_and_Temp_Mo[Routes_and_Temp_Mo['Distance (miles)'] < Range_mi ] 
        Infeasible_Routes_1  = Routes_and_Temp_Mo[Routes_and_Temp_Mo['Distance (miles)'] > Range_mi ]
        
        # filter for flights based on airline   
        Feasible_Routes_2    = Feasible_Routes_1  
        Infeasible_Routes_2  = Feasible_Routes_1.head(0)  
    
        # filter for flights based on aircraft type 
        Feasible_Routes_3    = Feasible_Routes_2  
        Infeasible_Routes_3  = Feasible_Routes_2.head(0) 
        
        # filter for flights that have no passengers due to too much battery weight 
        Feasible_Routes_4    = Feasible_Routes_3[Feasible_Routes_3['E_Passengers'] > 0 ] 
        Infeasible_Routes_4  = Feasible_Routes_3[Feasible_Routes_3['E_Passengers'] < 0 ]
        
        # filter for flights where batteries cannot operate 
        Feasible_Routes_5    = Feasible_Routes_4[Feasible_Routes_4['Origin ' + month] > Min_Temp] 
        Infeasible_Routes_5  = Feasible_Routes_4[Feasible_Routes_4['Origin ' + month] < Min_Temp] 
        Feasible_Routes_6    = Feasible_Routes_5[Feasible_Routes_5['Origin ' + month] < Max_Temp] 
        Infeasible_Routes_6  = Feasible_Routes_5[Feasible_Routes_5['Origin ' + month] > Max_Temp]  
        Feasible_Routes_7    = Feasible_Routes_6[Feasible_Routes_6['Destination ' + month] > Min_Temp] 
        Infeasible_Routes_7  = Feasible_Routes_6[Feasible_Routes_6['Destination ' + month] < Min_Temp] 
        Feasible_Routes_8    = Feasible_Routes_7[Feasible_Routes_7['Destination ' + month] < Max_Temp] 
        Infeasible_Routes_8  = Feasible_Routes_7[Feasible_Routes_7['Destination ' + month] > Max_Temp] 
        Feasible_Routes      = Feasible_Routes_8.head(int(len(Feasible_Routes_8)*percent_adoption/100 )) 
        Infeasible_Routes_9  = Feasible_Routes_8.tail(int(len(Feasible_Routes_8)*(100 - percent_adoption)/100 ))
        Infeasible_Routes    = pd.concat([Infeasible_Routes_1,Infeasible_Routes_2,Infeasible_Routes_3,
                                          Infeasible_Routes_4,Infeasible_Routes_5,Infeasible_Routes_6,
                                          Infeasible_Routes_7,Infeasible_Routes_8,Infeasible_Routes_9])    
   
        # EMISSIONS 
        # Compute emissions when electric aircraft are integrated into fleet  
        Infeasible_Route_fuel_volume          = np.sum(np.array(Infeasible_Routes['Total Fuel Per Route (Gal)'])) * gallons_to_Liters * liters_to_cubic_meters
        Emissions_w_Electric_Aircraft[m_i]    = kg_to_Megaton * JetA_GHG * Infeasible_Route_fuel_volume * density_JetA
        
        # Compute emissions without electric aircraft   
        Conventional_Air_Travel_fuel_volume   = np.sum(np.array(Routes_and_Temp_Mo['Total Fuel Per Route (Gal)'])) * gallons_to_Liters * liters_to_cubic_meters
        Emissions_wo_Electric_Aircraft[m_i]   = kg_to_Megaton * JetA_GHG * Conventional_Air_Travel_fuel_volume * density_JetA 
    
        # COST PER SEAT MILE
        # CASM for normal operations without electric aircraft 
        ASM_jet_A                = np.sum(Routes_and_Temp_Mo['Distance (miles)'] * Routes_and_Temp_Mo['Passengers'])
        Total_Fuel_Cost_jet_A    = np.sum(Routes_and_Temp_Mo['Fuel Cost']) 
        CASM_wo_E_Aircraft[m_i]  = 100*Total_Fuel_Cost_jet_A/ASM_jet_A    
        
        # CASM for when electric aircraft are integrated into fleet 
        if len(Feasible_Routes['E_Passengers']) == 0: 
            CASM_w_E_Aircraft[m_i]         = 0
        else:  
            Electricity_Cost   = sum((np.array(Feasible_Routes['Aircraft_Battery_Energy'])*joule_to_kWh*cost_of_electricity) * np.array(Feasible_Routes['No of Flights Per Month']) )
            ASM_electric       = sum(Feasible_Routes['Distance (miles)'] *Feasible_Routes['E_Passengers']  )
            CASM_w_E_Aircraft[m_i] = 100*Electricity_Cost/ASM_electric
            
    #================================================================================================================================================  
    # Plot Routes 
    #================================================================================================================================================          
    # Routes     
    fig_4                = go.Figure()
    color_1              = px.colors.qualitative.T10[7]
    color_2              = px.colors.qualitative.Pastel[10]
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
  
    fig_4.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            line = dict(width = 0.1,color = color_2),))
    
    lons       = np.empty(3 * len(Feasible_Routes))
    lons[::3]  = Feasible_Routes['Origin Longitude (Deg.)']
    lons[1::3] = Feasible_Routes['Destination Longitude (Deg.)']
    lons[2::3] = None
    lats       = np.empty(3 * len(Feasible_Routes))
    lats[::3]  = Feasible_Routes['Origin Latitude (Deg.)']
    lats[1::3] = Feasible_Routes['Destination Latitude (Deg.)']
    lats[2::3] = None    
  
    fig_4.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            line = dict(width = 2,color = "aquamarine"),)) 
 
    # Airports  
    fig_4.add_trace(go.Scattergeo( 
        lon = Routes_and_Temp['Destination Longitude (Deg.)'],
        lat = Routes_and_Temp['Destination Latitude (Deg.)'], 
        text = Routes_and_Temp['Destination City Name'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color,))) 

    fig_4.add_trace(go.Scattergeo( 
        lon = Routes_and_Temp['Origin Longitude (Deg.)'],
        lat = Routes_and_Temp['Origin Latitude (Deg.)'], 
        text = Routes_and_Temp['Origin City Name'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color, ))) 
     
    # Flight Paths 
    fig_4.update_layout(mapbox_style  = "open-street-map",      
                      showlegend    = False, 
                      height        = 400, 
                      geo_scope     ='usa',
                      margin        = {'t':0,'l':0,'b':0,'r':0},  
                      mapbox        = dict( accesstoken=mapbox_access_token,style=map_style,
                                            center=go.layout.mapbox.Center( lat=30, lon= 230 ))  )     
    #================================================================================================================================================      
    # Passenger vs Distance Traveled 
    #================================================================================================================================================     
    fig_5               = go.Figure() 
    fig_5.add_trace(go.Histogram(histfunc="sum",
                               x= Feasible_Routes['Distance (miles)'],
                               y = Feasible_Routes['E_Passengers'],
                               name='Electric', 
                               xbins=dict(start=0, end=4000, size=500),
                               marker_color=color_1,))
    fig_5.add_trace(go.Histogram(histfunc="sum",
                               x= Infeasible_Routes['Distance (miles)'],
                               y = Infeasible_Routes['Passengers'],
                               name='Fossil Fuel',
                               xbins=dict(start=0, end=4000, size=500),
                               marker_color=color_2,)) 
    
    # The two histograms are drawn on top of another
    fig_5.update_layout(barmode='stack', 
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
    fig_6 = go.Figure()
    Airport_Routes     = Feasible_Routes[['E_Passengers','Origin Airport','Destination City Name']]
    Cumulative_Flights = Airport_Routes.groupby(['Origin Airport']).sum()
    Busiest_Airports   = Cumulative_Flights.sort_values(by=['E_Passengers'], ascending = False).head(10) 
    Alphabetical_List  = Busiest_Airports.sort_values(by=['Origin Airport'])  
    fig_6.add_trace(go.Bar( x=list(Alphabetical_List['E_Passengers'].index),
                       y=np.array(Alphabetical_List['E_Passengers']),
                       marker_color=color_1)) 
    fig_6.update_layout(xaxis_title_text='Airport', 
                      yaxis_title_text='Passengers', 
                      height        = 300, 
                      width         = 600, 
                      margin        = {'t':0,'l':0,'b':0,'r':0},  
                      bargap        = 0.1,
                      font=dict(  size=font_size ))  

    #================================================================================================================================================      
    # Determine Ratio of Electrified to Jet-A Routes
    #================================================================================================================================================    
    fig_7                       = go.Figure()
    sector_colors               = [color_1,color_2]
    Feasible_Passenger_Miles    = np.sum(np.array(Feasible_Routes['E_Passengers'])* np.array(Feasible_Routes['Distance (miles)']))
    Infeasible_Passenger_Miles  = np.sum(np.array(Infeasible_Routes[['Passengers']])* np.array(Infeasible_Routes[['Distance (miles)']])) 
    labels                      = ["Electric", "Fossil Fuel"] 
    fig_7.add_trace(go.Pie(labels=labels,
                         values=[Feasible_Passenger_Miles, Infeasible_Passenger_Miles],
                         marker_colors=sector_colors)) 
    fig_7.update_traces(hole=.4, hoverinfo="label+percent+name") 
    fig_7.update_layout( height     = 400, 
                      width         = 600, 
                      margin        = {'t':50,'l':0,'b':0,'r':0},  
                      font=dict(  size=font_size ))  
     
       
    
    #================================================================================================================================================      
    # Monthly Emissions
    #================================================================================================================================================                  
    month_names         = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']      
    fig_8               = go.Figure() 
    fig_8.add_trace(go.Scatter(x=month_names, y=Emissions_w_Electric_Aircraft, name = 'Aircraft Fleet with Electric Aircraft',
                             line=dict(color=color_1, width=4)))  
    fig_8.add_trace(go.Scatter(x=month_names, y=Emissions_wo_Electric_Aircraft, name='Aircraft Fleet without Electric Aircraft',
                             line=dict(color=color_2, width=4)))   
    fig_8.update_layout( 
                      height           = 400, 
                      width            = 600, 
                      margin           = {'t':50,'l':0,'b':0,'r':0},
                      yaxis_title_text ='CO2 Emissions (MtCO2)', 
                      yaxis_range      = [0,20],
                      font=dict(  size=font_size ),
                      legend=dict(
                          yanchor="top",
                          y=0.99,
                          xanchor="center",
                          x=0.4 )) 
    
    #================================================================================================================================================      
    # Cost Per Seat Mile
    #================================================================================================================================================   
    fig_9 = go.Figure()       
    fig_9.add_trace(go.Scatter(x=month_names, y=CASM_w_E_Aircraft, name = 'Electric',
                             line=dict(color=color_1, width=4)))  
    fig_9.add_trace(go.Scatter(x=month_names, y=CASM_wo_E_Aircraft, name='Fossil Fuel',
                             line=dict(color=color_2, width=4)))  
    fig_9.update_layout( 
                      height           = 400, 
                      width            = 600, 
                      margin           = {'t':50,'l':0,'b':0,'r':0},
                      yaxis_title_text ='Energy Source CASM (cents)', 
                      yaxis_range      = [0,1.2 * np.maximum(max(CASM_wo_E_Aircraft), max(CASM_w_E_Aircraft))],
                      font=dict(  size=font_size ),
                      legend=dict(
                          yanchor="top",
                          y=0.99,
                          xanchor="center",
                          x=0.4 )) 
    
    #fig_4["layout"]["template"] = template   
    #fig_5["layout"]["template"] = template         
    #fig_6["layout"]["template"] = template    
    #fig_7["layout"]["template"] = template 
    #fig_8["layout"]["template"] = template  
    #fig_9["layout"]["template"] = template  
    
    return fig_4,fig_5,fig_6,fig_7,fig_8,fig_9
 

     
if __name__ == '__main__':
    main()
