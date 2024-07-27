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
    
    SAT_data                       = pd.read_excel(H2_filename,sheet_name=['Hydrogen']) 
    Hydrogen                       = SAT_data['Hydrogen']   
    a                              = Hydrogen['Feedstock']
    b                              = Hydrogen['Production Technology'] 
    Hydrogen["H2 Fuel Name"]  =  a  + ' via ' + b 
    Flight_Ops                     = pd.read_excel(routes_filename,sheet_name=['Sheet1']) 
    Flight_Ops                     = Flight_Ops['Sheet1']
    

    H2_1                  = Hydrogen['H2 Fuel Name'][2] 
    H2_2                  = Hydrogen['H2 Fuel Name'][1] 
    H2_3                  = Hydrogen['H2 Fuel Name'][5] 
    H2_4                  = Hydrogen['H2 Fuel Name'][10]     
    selected_fuels         = [H2_1,H2_2,H2_3,H2_4]     
    volume_fraction       = 10  
    selected_airpots      = "All Airports" # ["Top 5 Airports","Top 10 Airports","Top 20 Airports", "Top 50 Airports",  "All Airports"]
    percent_adoption      = 100 
    month_no              = 1
    altitude              = 36000
    mean_SFC_Imperial     = 0.6 # lbm/lbf-hr
    switch_off            = False  
    H2_dollars_per_gal    = 8.50
    
    
    slider_bar_options      = get_H2_distribution(Hydrogen,selected_fuels)
    fig  =  generate_H2_slider_bar(Hydrogen,selected_fuels,slider_bar_options)
    fig.show()
    
    fig_1, fig_2, fig_3, fig_4,fig_6 ,fig_8 = generate_electric_flight_operations_plots(Flight_Ops,Hydrogen,mean_SFC_Imperial, selected_fuels,altitude,selected_airpots,slider_bar_options, volume_fraction,percent_adoption,month_no,H2_dollars_per_gal ,switch_off)

    fig_1.show()
    fig_2.show()
    fig_3.show()
    fig_4.show()
    fig_6.show()
    fig_8.show()                
      
    return 

def get_H2_distribution(Hydrogen,selected_fuels): 
    mask                  = Hydrogen['H2 Fuel Name'].isin(selected_fuels)
    fuels_used            = Hydrogen[mask]
    
    # get unique colors 
    n =  len(fuels_used['Color'])
    
    # create list of slider bar options 
    percent_H2_color   =  np.linspace(0, 100, n+1)
    slider_bar_options =  list(percent_H2_color[1:-1]) 
    
    return  slider_bar_options

def Density_given_height(h_G): 
    r_e=3959*5280 #miles to feet
    R=1716.5 #ft2/R-sec
    g0=32.174 #ft/s**2
    T0=518.69 #deg R     
    h=(r_e/(r_e+h_G))*h_G
     
    if h  <= 36000:
        h0=0
        T0=518.69 #deg R 
        rho0=2.3769e-3 #slugs/ft3 
        a1=-3.57/1000 #deg R/ft
        T = T0 + a1*(h -h0) 
        rho  = rho0*(T/T0)**(-((g0/(R*a1))+1))
    else:
        h0=36000
        rho0=7.103559955456123e-04 #from running code at 36000
        T=389.99 
        rho = rho0*np.exp((-g0/(R*T))*(h-h0))  
     
    rho_SI =  rho *  515.378818 # conversion to kg / m3
    return rho_SI   


     
def generate_electric_flight_operations_plots(Flight_Ops,Hydrogen,mean_SFC_Imperial,selected_fuels,altitude,selected_airpots,percent_H2_color , volume_fraction_unnormalized,percent_adoption,month_no,H2_dollars_per_gal ,switch_off): 
    mapbox_access_token  = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ"     
    map_style            = None if switch_off else 'dark'  
    #template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"] 
    font_size            = 16      
 
    rho                = Density_given_height(altitude)
     
    density_JetA       = 820.0  # kg /m3  
    density_H2         = 70     # kg /m3  
    gallons_to_Liters  = 3.78541
    liters_to_cubic_meters =  0.001
    Jet_A_GHG          = 4.36466 # CO2e/kg fuel 
    mean_SFC_SI        = mean_SFC_Imperial * 2.832545029426566e-05
    SFC_H2             = mean_SFC_SI * 0.35
      
    #================================================================================================================================================  
    # Compute Feasible Routes 
    #================================================================================================================================================
    # Step 1: extract variables from data sheet  
    original_pax_capacity = np.array(Flight_Ops['Estimated Aircraft Capacity'])
    
    # Step 2: Use regression of passengers, available fuselage volume and weight to  compute h2 volume  
    aircraft_volume  = 0.0003 * (original_pax_capacity**3) - 0.0893*(original_pax_capacity**2) + 10.1*(original_pax_capacity **1) - 314.21
    volume_fraction  = volume_fraction_unnormalized / 100
    H2_volume        = volume_fraction *aircraft_volume
    Flight_Ops['H2_volume'] = H2_volume 
     
    # Step 3: Compute max range of each aircraft using range equation 
    g                    = 9.81
    Weight_empty         = (1.2636* (original_pax_capacity**2) + 145.58* (original_pax_capacity)) *g
    wing_area            = 0.0059* (original_pax_capacity**2)- 0.9942* (original_pax_capacity) + 132.03 
    CL                   = -6E-06* (original_pax_capacity**2) + 0.0028* (original_pax_capacity)+ 0.2695  
    CD                   = 8E-07* (original_pax_capacity**2) - 0.0003* (original_pax_capacity) + 0.0615 
    mass_pass            = 250
    Weight_Pass_total    = (1 -  volume_fraction) * mass_pass *g *  original_pax_capacity 
    Weight_H2            = H2_volume * density_H2 * g
    Weight_H2_0          = Weight_empty + Weight_Pass_total  + Weight_H2 
    Weight_H2_f          = Weight_H2_0 - (Weight_H2*0.95) # only 90% of fuel is used up 
    Range                = (1 /SFC_H2) *((CL**2)/CD) *np.sqrt(2/(rho*wing_area))*( np.sqrt(Weight_H2_0) -  np.sqrt(Weight_H2_f))
    Range_mi             = Range * 0.000621371 # conversion to mile 
    
    # Step 4: Compute modified passenger capacity   
    Flight_Ops['H2_Passengers']  = (1 - volume_fraction) *  np.array(Flight_Ops['Passengers']) 
       
    # Step 5: compute the percentages of different types of fuels used 
    H2_process_percentages_list  = [0]
    H2_process_percentages_list  += percent_H2_color
    H2_process_percentages_list  += [100]
    H2_process_percentages       = np.diff(np.array(H2_process_percentages_list))/100   
     
    # Step 6: determine the percentage of different H2 production processes 
    mask                  = Hydrogen['H2 Fuel Name'].isin(selected_fuels)
    fuels_used            = Hydrogen[mask] 
      
    # Step 7: Filter flight data based on option selected: i.e. top 10, top 20, top 50, all airpots 
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
    
    # Step 8: Filter airports that will support H2 and those that wont support H2 
    mask_1               = Flight_Ops['Origin Airport'].isin(Airport_List)
    H2_Airports          = Flight_Ops[mask_1]
    Non_H2_Airports      = Flight_Ops[~mask_1] 
    Feasible_Routes_0    = H2_Airports[H2_Airports['H2_Passengers'] > 0 ] 
    Infeasible_Routes_0  = H2_Airports[H2_Airports['H2_Passengers'] < 0 ]   
    Feasible_Routes_1    = Feasible_Routes_0[Feasible_Routes_0['Distance (miles)'] < Range_mi ] 
    Infeasible_Routes_1  = Feasible_Routes_0[Feasible_Routes_0['Distance (miles)'] > Range_mi ] 
    
    # Step 9: Out of H2 supporting airports, use the percent adoption to determine how many flights at that airport will use H2     
    Flight_at_H2_Airports_Using_H2  = Feasible_Routes_1.sample(frac=(percent_adoption/100))
    Infeasible_Routes_2             = Feasible_Routes_1[~Feasible_Routes_1.index.isin(Flight_at_H2_Airports_Using_H2.index)]
    Non_H2_Flights                  = pd.concat([Non_H2_Airports,Infeasible_Routes_0,Infeasible_Routes_1,Infeasible_Routes_2])# add list of flights from non supporting airports to non-H2 flights  
        
    # Steps 11 and 12 Determine Cost per Seat Mile and Emissions  
    CASM_jet_A      = np.zeros(12) 
    CASM_H2         = np.zeros(12) 
    CO2e_w_H2       = np.zeros(12) 
    CO2e_w_o_H2     = np.zeros(12) 
    for m_i in range(12): 
        Flight_at_H2_Airports_Using_H2_Mo  = Flight_at_H2_Airports_Using_H2.loc[Flight_at_H2_Airports_Using_H2['Month'] == m_i+1]
        Flight_Mo                          = Flight_Ops.loc[Flight_Ops['Month'] == m_i+1]
        Total_jet_A_no_H2                  = np.sum(np.array(Flight_Mo['Total Fuel Per Route (Gal)']))    
        
        H2_volumes_mo                      = np.sum(np.array(Flight_at_H2_Airports_Using_H2_Mo['H2_volume']))
        
        Non_H2_Flights_Mo                  = Non_H2_Flights.loc[Non_H2_Flights['Month'] == m_i+1 ]   
        Jet_A_fuel_volume_required_mo      = np.sum(np.array(Non_H2_Flights_Mo['Total Fuel Per Route (Gal)']))    
        
        
        if len(Non_H2_Flights_Mo ) == 0:
            CASM_jet_A[m_i] = 0
        else:
            ASM_jet_A               = np.sum(Non_H2_Flights_Mo['Distance (miles)'] * Non_H2_Flights_Mo['Passengers'])
            Total_Fuel_Cost_jet_A   = np.sum(Non_H2_Flights_Mo['Fuel Cost']) 
            CASM_jet_A[m_i]         = 100*Total_Fuel_Cost_jet_A/ASM_jet_A     
       
        if len(Flight_at_H2_Airports_Using_H2_Mo)  == 0:
            CASM_H2[m_i]  = 0 
        else:    
            ASM_H2               = np.sum(Flight_at_H2_Airports_Using_H2_Mo['Distance (miles)'] * Flight_at_H2_Airports_Using_H2_Mo['H2_Passengers']) 
            Total_Fuel_Cost_H2   = H2_volumes_mo /liters_to_cubic_meters /gallons_to_Liters  * H2_dollars_per_gal 
            CASM_H2[m_i]         = 100*Total_Fuel_Cost_H2/ASM_H2   
        
        CO2e_w_H2_only_H2    = np.sum(np.array(fuels_used['Direct GHG emissions [kg CO2e/kg H2]'])*H2_process_percentages)*H2_volumes_mo
        CO2e_w_H2_fuel_only  = Jet_A_GHG * Jet_A_fuel_volume_required_mo * gallons_to_Liters * liters_to_cubic_meters *  density_JetA 
        CO2e_w_H2[m_i]       = CO2e_w_H2_only_H2 +  CO2e_w_H2_fuel_only 
        CO2e_w_o_H2[m_i]     = Jet_A_GHG * Total_jet_A_no_H2 
           
    
    #================================================================================================================================================  
    # Plot Flight_Ops 
    #================================================================================================================================================     
    # Flight_Ops   
    colors               = px.colors.qualitative.Pastel 
    colors2              = px.colors.qualitative.Safe
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
            line = dict(width = 1,color = colors2[0]), ))
  
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
                               marker_color=colors2[0],))
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
                       marker_color=colors2[0])) 
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
    sector_colors               = [colors2[0],colors[2]]
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
                             line=dict(color= colors2[0], width=4)))  
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
    fig_8.add_trace(go.Scatter(x=month_names, y=CO2e_w_H2*1e-6, name = 'H2 Aircraft in Fleet',
                             line=dict(color=colors2[0], width=4)))  
    fig_8.add_trace(go.Scatter(x=month_names, y=CO2e_w_o_H2*1e-6, name='No H2 Aircraft in Fleet',
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
     
    #fig_1["layout"]["template"] = template  
    #fig_2["layout"]["template"] = template
 
    
    return fig_1, fig_2, fig_3, fig_4,fig_6 ,fig_8


# ---------------------------------------------------------------------------------------------------------------------------------------------------
# Motor Plots
# ---------------------------------------------------------------------------------------------------------------------------------------------------
def generate_H2_slider_bar(Hydrogen,selected_fuels,H2_values):


    mask                  = Hydrogen['H2 Fuel Name'].isin(selected_fuels)
    fuels_used            = Hydrogen[mask]
    
    # get unique colors 
    n =  len(fuels_used['Color'])
    
    v    = [0] 
    v    = v + H2_values
    v    = v + [100]
    vals = np.array(v) 
    
    fig = go.Figure(go.Histogram(
        y=np.random.randint(1, size= int(100)),
        marker_color=fuels_used['Color'].iloc[0],
        name        =fuels_used['H2 Fuel Name'].iloc[0], 
        bingroup=1))
    
    for i in  range(1, n): 
        fig.add_trace(go.Histogram(
            y=np.random.randint(1, size=int(100-vals[i])),
            marker_color=fuels_used['Color'].iloc[i],
            name        =fuels_used['H2 Fuel Name'].iloc[i], 
            bingroup=1))    
    
    fig.update_layout(
        barmode="overlay",
        bargap=0.1)
    
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                      xaxis_title=None,
                      yaxis_title=None, 
                      height     = 100, 
                      width      = 500, 
                      showlegend=False,
                      xaxis     =  {'title': 'x-label','visible': False,'showticklabels': True},
                      yaxis     =  {'title': 'y-label','visible': False,'showticklabels': False},
                              ) 
    return fig 


     
if __name__ == '__main__':
    main() 
