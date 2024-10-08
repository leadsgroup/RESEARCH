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
    
   
    separator                  = os.path.sep
    airline                    =  'Top_10_Major_US_Airlines'  
    H2_filename                = 'Data' + separator + 'Technology' + separator +  'Technology_Data.xlsx'
    route_temp_filename        = 'Data' + separator +  'Air_Travel' + separator + airline +  '_Flight_Ops_and_Climate.csv'   
    SAT_data                   = pd.read_excel(H2_filename,sheet_name=['Hydrogen']) 
    Hydrogen                   = SAT_data['Hydrogen']   
    a                          = Hydrogen['Feedstock']
    b                          = Hydrogen['Production Technology'] 
    c                          = Hydrogen['Production Process'] 
    Hydrogen["H2 Fuel Name"]   =  a  + ' via ' + b +  ' using ' + c
    Flight_Ops                 = pd.read_csv(route_temp_filename)  
    

    H2_1                  = Hydrogen['H2 Fuel Name'][2] 
    H2_2                  = Hydrogen['H2 Fuel Name'][1] 
    H2_3                  = Hydrogen['H2 Fuel Name'][5] 
    H2_4                  = Hydrogen['H2 Fuel Name'][10]     
    selected_h2           = [H2_1,H2_2,H2_3,H2_4]      
    h2_vol_percentage     = 15  
    h2_airports           = " Top 10 Airports" # DONE  # [" Top 5 Airports"," Top 10 Airports"," Top 20 Airports", " Top 50 Airports",  "All Airports"]
    h2_percent_adoption   = 100    
    h2_cruise_alt         = 35000  
    mean_SFC_Imperial     = 0.6   
    H2_dollars_per_gal    = 8  
    switch_off            = False  
    
    
    percent_H2_color      = get_H2_distribution(Hydrogen,selected_h2)
    fig  =  generate_H2_slider_bar(Hydrogen,selected_h2,percent_H2_color)
    fig.show()
    
    
    fig_1, fig_2, fig_3, fig_4,fig_5 ,fig_6 = generate_electric_flight_operations_plots(Flight_Ops,Hydrogen,selected_h2,mean_SFC_Imperial,h2_cruise_alt,h2_airports,percent_H2_color,h2_vol_percentage,h2_percent_adoption,H2_dollars_per_gal,switch_off)
 

    fig_1.show()
    fig_2.show()
    fig_3.show()
    fig_4.show()
    fig_5.show()
    fig_6.show()                
      
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
    r_e=3959*5280 # miles to feet
    R=1716.5      # ft2/R-sec
    g0=32.174     # ft/s**2
    T0=518.69     # deg R     
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

#def generate_electric_flight_operations_plots(Flight_Ops,Hydrogen,mean_SFC_Imperial,selected_fuels,altitude,selected_airpots,percent_H2_color,volume_fraction_unnormalized,percent_adoption,H2_dollars_per_gal,switch_off): 


def generate_electric_flight_operations_plots(Flight_Ops,Hydrogen,selected_h2,mean_SFC_Imperial,h2_cruise_alt,h2_airports,percent_H2_color,h2_vol_percentage,h2_percent_adoption,H2_dollars_per_gal,switch_off): 
    mapbox_access_token  = "pk.eyJ1IjoibWFjbGFya2UiLCJhIjoiY2xyanpiNHN6MDhsYTJqb3h6YmJjY2w5MyJ9.pQed7pZ9CnJL-mtqm1X8DQ"     
    map_style            = None if switch_off else 'dark'  
    #template             = pio.templates["minty"] if switch_off else pio.templates["minty_dark"] 
    font_size            = 16    
     
    #================================================================================================================================================  
    # Unit Conversions 
    #================================================================================================================================================    
    rho                    = Density_given_height(h2_cruise_alt) 
    density_JetA           = 820.0  # kg /m3  
    density_H2             = 70     # kg /m3  
    gallons_to_Liters      = 3.78541
    liters_to_cubic_meters = 0.001
    JetA_GHG               = 4.36466 # CO2e/kg fuel 
    imperial_to_SI_SFC     = 2.832545029426566e-05
    g                      = 9.81
    mile_to_meters         = 0.000621371 
    passenger_mass         = 250  # average mass of passenger and their luggage in kg
    kg_to_Megaton          = 1E-9
    SFC_H2                 = mean_SFC_Imperial * imperial_to_SI_SFC 
      
    #================================================================================================================================================  
    # Compute Feasible Routes 
    #================================================================================================================================================
    # Extract variables from data sheet  
    original_pax_capacity = np.array(Flight_Ops['Estimated Aircraft Capacity'])
    
    # Use regression of passengers, available fuselage volume and weight to  compute h2 volume  
    aircraft_volume         = 0.0003 * (original_pax_capacity**3) - 0.0893*(original_pax_capacity**2) + 10.1*(original_pax_capacity **1) - 314.21
    Weight_empty            = (1.2636* (original_pax_capacity**2) + 145.58* (original_pax_capacity)) *g
    wing_area               = 0.0059* (original_pax_capacity**2)- 0.9942* (original_pax_capacity) + 132.03 
    CL                      = -6E-06* (original_pax_capacity**2) + 0.0028* (original_pax_capacity)+ 0.2695  
    CD                      = 8E-07* (original_pax_capacity**2) - 0.0003* (original_pax_capacity) + 0.0615
    
    # Normalize volume fraction 
    volume_fraction         = h2_vol_percentage / 100
    
    # Apply volume fraction 
    H2_volume               = volume_fraction *aircraft_volume
    Flight_Ops['H2_volume'] = H2_volume 
     
    # Compute max range of each aircraft using range equation 
    Weight_Pass_total    = (1 -  volume_fraction) * passenger_mass *g *  original_pax_capacity 
    Weight_H2            = H2_volume * density_H2 * g
    Weight_H2_0          = Weight_empty + Weight_Pass_total  + Weight_H2 
    Weight_H2_f          = Weight_H2_0 - (Weight_H2*0.95) # only 90% of fuel is used up 
    Range                = (1 /SFC_H2) *((CL**2)/CD) *np.sqrt(2/(rho*wing_area))*( np.sqrt(Weight_H2_0) -  np.sqrt(Weight_H2_f))
    Range_mi             = Range * mile_to_meters # conversion to mile 
    
    # Compute the percentage of the aircraft with seated passengers
    percentage_capacity_of_aircraft =  np.array(Flight_Ops['Passengers']) / ( np.array(Flight_Ops['No of Flights Per Month']) * np.array(Flight_Ops['Estimated Aircraft Capacity']))
    
    # Determine unused passenger volume 
    unused_volume_of_aircraft =  np.ones_like(percentage_capacity_of_aircraft) - percentage_capacity_of_aircraft
    
    # Subtract unused volume of aircraft from the volume fraction specified
    additiona_volume_fractionl_required =  unused_volume_of_aircraft -  np.ones_like(percentage_capacity_of_aircraft)*volume_fraction
    
    # If required volume is negative, we do not need to remove passengers 
    additiona_volume_fractionl_required[additiona_volume_fractionl_required>0] =  0
    
    # If additiona_volume_fractionl_required is negative, we need to take off more passenger to accomodate H2 volume:
    # Determine how many additional passengers must be removed from aircraft to meet specified volume fraction
    Flight_Ops['H2_Passengers']  =   np.array(Flight_Ops['Passengers']) -  (1 + additiona_volume_fractionl_required) * np.array(Flight_Ops['Estimated Aircraft Capacity']) 
       
    # Compute the percentages of different types of fuels used 
    percentage_H2_process_vals  = [0]
    percentage_H2_process_vals  += percent_H2_color
    percentage_H2_process_vals  += [100]
    percentage_H2_process       = np.diff(np.array(percentage_H2_process_vals))/100   
     
    # Determine the percentage of different H2 production processes 
    selected_H2_mask      = Hydrogen['H2 Fuel Name'].isin(selected_h2)
    H2_used               = Hydrogen[selected_H2_mask] 
      
    # Filter flight data based on option selected: i.e. top 10, top 20, top 50, all airpots 
    Airport_Routes     = Flight_Ops[['Passengers','Origin Airport','Destination City Name']]
    Cumulative_Flights = Airport_Routes.groupby('Origin Airport', as_index=False).sum()[Airport_Routes.columns]
    if  h2_airports   == " Top 5 Airports":
        Busiest_H2_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(5) 
    if  h2_airports   == " Top 10 Airports":
        Busiest_H2_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(10) 
    elif  h2_airports ==" Top 20 Airports":
        Busiest_H2_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(20) 
    elif  h2_airports == " Top 50 Airports":
        Busiest_H2_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(50) 
    elif  h2_airports == " All Airports":
        Busiest_H2_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False) 
    Airport_List = list(Busiest_H2_Airports['Origin Airport'])
    
    # Filter airports that will support H2 and those that wont support H2  
    Feasible_Routes_1      = Flight_Ops[Flight_Ops['Origin Airport'].isin(Airport_List)]
    Infeasible_Routes_1    = Flight_Ops[~Flight_Ops['Origin Airport'].isin(Airport_List)]
    
    # Filter routes that do not meet the range of the aircraft 
    H2_Airports_Range    = Range_mi[Flight_Ops['Origin Airport'].isin(Airport_List)]   
    Feasible_Routes_2    = Feasible_Routes_1[Feasible_Routes_1['Distance (miles)'] < H2_Airports_Range ] 
    Infeasible_Routes_2  = Feasible_Routes_1[Feasible_Routes_1['Distance (miles)'] > H2_Airports_Range ]     

    # Filter out routes where there are no passengers 
    Feasible_Routes_3    = Feasible_Routes_2[Feasible_Routes_2['H2_Passengers'] > 0 ] 
    Infeasible_Routes_3  = Feasible_Routes_2[Feasible_Routes_2['H2_Passengers'] < 0 ]   
    
    # Out of H2 supporting airports, use the percent adoption to determine how many flights at that airport will use H2     
    H2_Flights            = Feasible_Routes_3.sample(frac=(h2_percent_adoption/100))
    Infeasible_Routes_4   = Feasible_Routes_3[~Feasible_Routes_3.index.isin(H2_Flights.index)] 
    Non_H2_Flights        = pd.concat([Infeasible_Routes_1,Infeasible_Routes_2,Infeasible_Routes_3,Infeasible_Routes_4])# add list of flights from non supporting airports to non-H2 flights  
        
    # Determine Cost per Seat Mile and Emissions  
    CASM_wo_H2_Aircraft     = np.zeros(12) 
    CASM_w_H2_Aircraft      = np.zeros(12) 
    CO2e_w_H2_Aircraft      = np.zeros(12) 
    CO2e_wo_H2_Aircraft     = np.zeros(12) 
    for m_i in range(12): 
        H2_Flights_Mo                = H2_Flights.loc[H2_Flights['Month'] == m_i+1] 
        Non_H2_Flights_Mo            = Non_H2_Flights.loc[Non_H2_Flights['Month'] == m_i+1 ]   
        Non_H2_Flights_JetA_vol_mo   = np.sum(np.array(Non_H2_Flights_Mo['Total Fuel Per Route (Gal)']))   * gallons_to_Liters * liters_to_cubic_meters   
        H2_volumes_mo                = np.sum(np.array(H2_Flights_Mo['H2_volume']))
        
        if len(Non_H2_Flights_Mo ) == 0:
            pass
        else:
            ASM_JetA               = np.sum(Non_H2_Flights_Mo['Distance (miles)'] * Non_H2_Flights_Mo['Passengers'])
            Total_Fuel_Cost_JetA   = np.sum(Non_H2_Flights_Mo['Fuel Cost']) 
            CASM_wo_H2_Aircraft[m_i]         = 100*Total_Fuel_Cost_JetA/ASM_JetA # convert to cents/pass-mile
       
        if len(H2_Flights_Mo)  == 0:
            pass
        else:    
            ASM_H2               = np.sum(H2_Flights_Mo['Distance (miles)'] * H2_Flights_Mo['H2_Passengers']) 
            Total_Fuel_Cost_H2   = H2_volumes_mo /liters_to_cubic_meters /gallons_to_Liters  * H2_dollars_per_gal 
            CASM_w_H2_Aircraft[m_i]         = 100*Total_Fuel_Cost_H2/ASM_H2   # convert to cents/pass-mile
        
        # CO2e of the Hydrogen-powered flights in the scenario with flight operations with H2 tech
        CO2e_H2_aircraft_in_JetAH2_Fleet    = kg_to_Megaton * np.sum(np.array(H2_used['Direct GHG emissions [kg CO2e/kg H2]'])*percentage_H2_process)*H2_volumes_mo
 
        # CO2e of the JetA-powered flights in the scenario with flight operations with H2 tech        
        CO2e_JetA_aircraft_in_JetAH2_Fleet  = kg_to_Megaton * JetA_GHG * Non_H2_Flights_JetA_vol_mo *  density_JetA 
        
        # Total CO2e flights in the scenario with flight operations with H2 tech
        CO2e_w_H2_Aircraft[m_i]       =  CO2e_H2_aircraft_in_JetAH2_Fleet +  CO2e_JetA_aircraft_in_JetAH2_Fleet
 
        # Total CO2e flights in the scenario with flight operations without H2 tech    
        Conventional_Air_Travel             = Flight_Ops.loc[Flight_Ops['Month'] == m_i+1]
        Conventional_Air_Travel_fuel_volume = np.sum(np.array(Conventional_Air_Travel['Total Fuel Per Route (Gal)'])) * gallons_to_Liters * liters_to_cubic_meters
        CO2e_wo_H2_Aircraft[m_i]            = kg_to_Megaton * JetA_GHG * Conventional_Air_Travel_fuel_volume * density_JetA        
         
    #================================================================================================================================================  
    # Plot Flight_Ops 
    #================================================================================================================================================     
    # Flight_Ops   
    color_1               = px.colors.qualitative.Pastel[2]
    color_2              = px.colors.qualitative.Alphabet[19]
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
            line = dict(width = 1,color = color_1), ))
    
    lons       = np.empty(3 * len(H2_Flights))
    lons[::3]  = H2_Flights['Origin Longitude (Deg.)']
    lons[1::3] = H2_Flights['Destination Longitude (Deg.)']
    lons[2::3] = None
    lats       = np.empty(3 * len(H2_Flights))
    lats[::3]  = H2_Flights['Origin Latitude (Deg.)']
    lats[1::3] = H2_Flights['Destination Latitude (Deg.)']
    lats[2::3] = None    
  
    fig_1.add_trace(
        go.Scattergeo( 
            lon = lons,
            lat = lats,
            mode = 'lines',
            opacity= 0.5,
            line = dict(width = 1,color = color_2), ))
  
    # Airports  
    fig_1.add_trace(go.Scattergeo( 
        lon = Flight_Ops['Destination Longitude (Deg.)'],
        lat = Flight_Ops['Destination Latitude (Deg.)'], 
        text = Flight_Ops['Destination City Name'],
        mode = 'markers',
        marker = dict(
            size = airport_marker_size,
            color = airport_marker_color, ))) 
 
    fig_1.add_trace(go.Scattergeo( 
        lon = Flight_Ops['Origin Longitude (Deg.)'],
        lat = Flight_Ops['Origin Latitude (Deg.)'], 
        text = Flight_Ops['Origin City Name'],
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
                               x= H2_Flights['Distance (miles)'],
                               y = H2_Flights['Passengers'],
                               name='Hydrogen', 
                               xbins=dict(start=0, end=4000, size=500),
                               marker_color=color_2,))
    fig_2.add_trace(go.Histogram(histfunc="sum",
                               x= Non_H2_Flights['Distance (miles)'],
                               y = Non_H2_Flights['Passengers'],
                               name='Jet-A',
                               xbins=dict(start=0, end=4000, size=500),
                               marker_color=color_1,)) 
    
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
    Airport_Routes     = H2_Flights[['Passengers','Origin Airport','Destination City Name']]
    Cumulative_Flights = Airport_Routes.groupby(['Origin Airport']).sum()
    Busiest_H2_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(10) 
    Alphabetical_List  = Busiest_H2_Airports.sort_values(by=['Origin Airport'])  
    fig_3.add_trace(go.Bar( x=list(Alphabetical_List['Passengers'].index),
                       y=np.array(Alphabetical_List['Passengers']),
                       marker_color=color_2)) 
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
    sector_colors               = [color_2,color_1]
    Feasible_Passenger_Miles    = np.sum(np.array(H2_Flights['Passengers'])* np.array(H2_Flights['Distance (miles)']))
    Infeasible_Passenger_Miles  = np.sum(np.array(Non_H2_Flights[['Passengers']])* np.array(Non_H2_Flights[['Distance (miles)']]))
    labels                      = ["Hydrogen", "Jet-A"] 
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
    fig_5 = go.Figure()       
    month_names         = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']    
    fig_5.add_trace(go.Scatter(x=month_names, y=CASM_w_H2_Aircraft, name = 'Hydrogen',
                             line=dict(color= color_2, width=4)))  
    fig_5.add_trace(go.Scatter(x=month_names, y=CASM_wo_H2_Aircraft, name='Jet-A',
                             line=dict(color= color_1, width=4)))  
    fig_5.update_layout( 
                      height           = 400, 
                      width            = 600, 
                      margin           = {'t':50,'l':0,'b':0,'r':0},
                      yaxis_title_text ='Energy Source CASM (cents)', 
                      font=dict(  size=font_size ),
                      yaxis_range      = [0,10],
                      legend=dict(
                          yanchor="top",
                          y=0.99,
                          xanchor="center",
                          x=0.4 )) 
 
    #================================================================================================================================================      
    # Emissions Comparison 
    #================================================================================================================================================   
    month_names         = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']      
    fig_6               = go.Figure() 
    fig_6.add_trace(go.Scatter(x=month_names, y=CO2e_w_H2_Aircraft, name = 'Domestic Air Travel with Hydrogen Aircraft',
                             line=dict(color=color_2, width=4)))  
    fig_6.add_trace(go.Scatter(x=month_names, y=CO2e_wo_H2_Aircraft, name='Domestic Air Travel without Hydrogen Aircraft',
                             line=dict(color=color_1, width=4)))   
    fig_6.update_layout( 
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
    
    return fig_1, fig_2, fig_3, fig_4,fig_5 ,fig_6 
     
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
