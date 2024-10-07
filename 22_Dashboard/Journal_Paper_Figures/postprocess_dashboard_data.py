
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTS 
# --------------------------------------------------------------------------------------------------------------------------------------------------- 

from RCAIDE.Framework.Core import  Units 
from RCAIDE.load import load as load_results
from RCAIDE.save import save as save_results 

import os 
import numpy as np           
import pandas as pd 
import json
from pyatmos import coesa76
import pickle
 
def main():
    # ---------------------------------------------------------------------------------------------------------------------------------------------------
    # RETRIEVE DATA 
    # --------------------------------------------------------------------------------------------------------------------------------------------------- 
    separator                  = os.path.sep 
    technology_filename        = '..' + separator +'Data'  + separator + 'Technology' + separator +  'Technology_Data.xlsx'
    crops_filename             = '..' + separator +'Data'  + separator +  'Crops'     + separator + 'All_Crops_2017.xlsx'   
    routes_filename            = '..' + separator +'Data'  + separator + 'Air_Travel' + separator + 'Top_10_Major_US_Airlines_Flight_Ops_and_Climate.csv' 
    temperature_filename       = '..' + separator +'Data' + separator  + 'US_Climate' + separator + 'Monthly_US_County_Temperature_2019.csv'
    SAT_data                   = pd.read_excel(technology_filename,sheet_name=['Commercial_Batteries','Battery_Development','Electric_Motor_Development','Commercial_SAF', 'Hydrogen']) 
    
    Commercial_Batteries       = SAT_data['Commercial_Batteries'] 
    Electric_Motor_Development = SAT_data['Electric_Motor_Development']
    a                          = Commercial_Batteries['Brand']  
    b                          = Commercial_Batteries['Abbreviation']
    c                          = Commercial_Batteries['Model']
    d                          = a + ': ' + b + '-' + c 
    Commercial_Batteries["Battery Name"] = d     
    Battery_Development        = SAT_data['Battery_Development']
    
    Commercial_SAF             = SAT_data['Commercial_SAF']  
    a                          = Commercial_SAF['Brand']  
    b                          = Commercial_SAF['Fuel Type']
    c                          = Commercial_SAF['Process']
    d                          = Commercial_SAF['Source']
    e                          = Commercial_SAF['Feedstock']
    Commercial_SAF["Fuel Name"]= ' ' + a + ' ' + c + '-' + b + ' from ' + e  + ' (' + d  + ')'
     
    Hydrogen                   = SAT_data['Hydrogen']   
    a                          = Hydrogen['Feedstock']
    b                          = Hydrogen['Production Technology'] 
    c                          = Hydrogen['Production Process'] 
    Hydrogen["H2 Fuel Name"]   =  a  + ' via ' + b +  ' using ' + c 
     
    Flight_Ops                 = pd.read_csv(routes_filename)    
    US_Temperature_F           = pd.read_csv(temperature_filename)   
    feedstocks                 = pd.read_excel(crops_filename,sheet_name=['Corn','Soybean','Canola','Sunflower','Sorghum','Wheat'])  

 
    
    generate_saf_results_data(Commercial_SAF,Flight_Ops,feedstocks)
    generate_eletrification_results_data(Flight_Ops)
    generate_hydrogen_results_data(Hydrogen,Flight_Ops) 
     
    return


def generate_saf_results_data(Commercial_SAF,Flight_Ops,feedstocks): 

    saf_1                           = Commercial_SAF['Fuel Name'][0] 
    saf_2                           = Commercial_SAF['Fuel Name'][1] 
    saf_3                           = Commercial_SAF['Fuel Name'][2] 
    saf_4                           = Commercial_SAF['Fuel Name'][3]     
    selected_fuels                  = [saf_1,saf_2,saf_3,saf_4]
    blend_ratios                    = [100,40,50,25]
    percent_fuel_use                = [75,80,95]   
    State_List_1                    = ["California","Colorado"]  
    State_List_2                    = ["Georgia","Idaho","Illinois", "Indiana","Iowa","Kansas","Kentucky","Louisiana"] 
    State_List_3                    = ["Missouri", "Nebraska"] 
    State_List_4                    = ["Ohio","Oklahoma"] 
    State_List_5                    = [""]         
    aircraft                        = "Boeing 787-8"
    selected_airpots                = "All Airports" # ["Top 5 Airports","Top 10 Airports","Top 20 Airports", "Top 50 Airports",  "All Airports"]
    percent_adoption                = 100 
    month_no                        = 1
    selected_feedstock              = 'Canola' 
    feedstock_producing_states      = State_List_1+State_List_2+State_List_3+State_List_4+State_List_5
     
    #================================================================================================================================================  
    # Unit Conversions 
    #================================================================================================================================================     
    JetA_GHG               = 4.36466 # CO2e/kg fuel 
    gallons_to_Liters      = 3.78541
    liters_to_cubic_meters = 0.001
    Jet_A_density          = 800.0  
    density_JetA           = 820.0  # kg /m3  
    kg_to_Megaton          = 1E-9 
    g_to_kg                = 0.001
    
    # Compute the percentages of different types of fuels used 
    fuel_percentages_list = [0]
    fuel_percentages_list += percent_fuel_use
    fuel_percentages_list += [100]
    fuels_percentages     = np.diff(np.array(fuel_percentages_list))/100  
    
    # Determine the percentage of neat (pure) saf and Jet-A1 using blending ratios   
    if Commercial_SAF['Fuel Name'][0] not in selected_fuels:  
        selected_fuels    = [Commercial_SAF['Fuel Name'][0]] + selected_fuels
        fuels_percentages = np.hstack((np.array([0]),fuels_percentages)) 
    mask                = Commercial_SAF['Fuel Name'].isin(selected_fuels)
    fuels_used          = Commercial_SAF[mask] 
    
    num_fuels = len(selected_fuels)
    cumulative_fuel_use   = np.zeros(num_fuels) 
    SAF_LCA_val           = np.zeros(num_fuels) 
    Jet_A_LCA_val         = np.zeros(num_fuels)  

    # Loop through fuels and get percentage of fuel used by each type    
    for i in range(1,num_fuels):
        blend_ratio             = np.array(Commercial_SAF.loc[Commercial_SAF['Fuel Name'] == selected_fuels[i]]['Maximum Blend Ratio']) 
        cumulative_fuel_use[i]  = fuels_percentages[i]* blend_ratio/100
        SAF_LCA_val[i]          = Commercial_SAF[Commercial_SAF['Fuel Name'] == selected_fuels[i]]['LCA Value']
    cumulative_fuel_use[0]      = 1 - np.sum(cumulative_fuel_use[1:])
    SAF_LCA_val[0]   = 89
    Jet_A_LCA_val[0] = 89
        
    # Filter flight data based on option selected: i.e. top 10, top 20, top 50, all airpots 
    Airport_Routes     = Flight_Ops[['Passengers','Origin Airport','Destination City Name']]
    Cumulative_Flights = Airport_Routes.groupby('Origin Airport', as_index=False).sum()[Airport_Routes.columns]
    if  selected_airpots == " Top 5 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(5) 
    if  selected_airpots == " Top 10 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(10) 
    elif  selected_airpots ==" Top 20 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(20) 
    elif  selected_airpots == " Top 50 Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False).head(50) 
    elif  selected_airpots == " All Airports":
        Busiest_Airports   = Cumulative_Flights.sort_values(by=['Passengers'], ascending = False) 
    Airport_List = list(Busiest_Airports['Origin Airport'])
    
    # Filter airports that will support SAF and those that wont support SAF 
    mask_1                = Flight_Ops['Origin Airport'].isin(Airport_List)
    SAF_Airports          = Flight_Ops[mask_1]
    Non_SAF_Airports      = Flight_Ops[~mask_1]  
    
    # Out of SAF supporting airports, use the percent adoption to determine how many flights at that airport will use SAF 
    Flight_at_SAF_Airports_Using_SAF      = SAF_Airports.sample(frac=(percent_adoption/100))
    Flight_at_SAF_Airports_Using_Jet_A    = SAF_Airports[~SAF_Airports.index.isin(Flight_at_SAF_Airports_Using_SAF.index)]
    Non_SAF_Flights                       = pd.concat([Non_SAF_Airports, Flight_at_SAF_Airports_Using_Jet_A] )   # add list of flights from non supporting airports to non-SAF flights  
     
    # Get total volume of each SAF required at the airports 
    total_fuel_volume_required = np.sum(np.array(Flight_at_SAF_Airports_Using_SAF['Total Fuel Per Route (Gal)']))
    fuel_volumes               = cumulative_fuel_use*total_fuel_volume_required  
    
    # Sort SAF's by feedstock and sum all fuel volumes based on fuel   
    fuels_used['Total Fuel Volume']   = fuel_volumes 
    relative_crop_area                = fuels_used['Total Fuel Volume']/fuels_used['SAF Gallons per Acre']
    fuels_used['Requires Acres']      = relative_crop_area
    select_fuels                      = fuels_used[['Source','Total Fuel Volume','Requires Acres']].groupby('Source', as_index=False).sum()[fuels_used[['Source','Total Fuel Volume','Requires Acres']].columns]
    required_crop_area                = np.array(select_fuels.loc[select_fuels['Source'] == selected_feedstock]['Requires Acres'])
    
    # Determine how many states will source the feedstock under consideration  
    crop_data  = feedstocks[selected_feedstock]   
    crop_data['FIPS'] = crop_data['FIPS'].apply('{:0>5}'.format)
    
    # Filter out states have not been selected and get number of states  
    mask                  = crop_data['State'].isin(feedstock_producing_states)
    feedstock_states      = crop_data[mask]    
    non_feedstock_states  = crop_data[~mask]   
    non_feedstock_states["Feedstock Usage"] = list(np.ones(len(non_feedstock_states)))

    # Randomize tracts in terms of crop area/usage then Recursively add rows until requied volume is met 
    Used_Feedstock        = feedstock_states.sample(frac = 1) 
    Used_Feedstock["Feedstock Usage"] =  np.ones(len(feedstock_states))
    
    idx        = 0
    total_vol  = 0    
    if len(required_crop_area) == 0:
        RCA  = 0
    else:
        RCA  = required_crop_area[0]  
    available_tracts =len(Used_Feedstock)
    while total_vol<RCA: 
        total_vol += Used_Feedstock.loc[Used_Feedstock.index[idx]]['Acres Harvested']
        Used_Feedstock["Feedstock Usage"][Used_Feedstock.index[idx]] = 0.1 
        idx += 1      
        if available_tracts == idx:
            total_vol = 1E9  
     
    # Determine Cost per Seat Mile and Emissions  
    CASM_wo_SAF_Aircraft  = np.zeros(12) 
    CASM_w_SAF_Aircraft    = np.zeros(12) 
    Emissions_w_SAF_Aircraft       = np.zeros(12) 
    Emissions_wo_SAF_Aircraft     = np.zeros(12) 
    gallons_to_Liters  = 3.78541
    for m_i in range(12): 
        Routes_and_Temp_Mo                  = Flight_Ops[Flight_Ops['Month'] == m_i+1 ]   
            
        # EMISSIONS 
        # Emissions of SAF-Scenario from non-SAF Flights 
        Non_SAF_Flights_Mo                  = Non_SAF_Flights.loc[Non_SAF_Flights['Month'] == m_i+1 ]   
        Infeasible_Routes_fuel_volume          = np.sum(np.array(Non_SAF_Flights_Mo['Total Fuel Per Route (Gal)'])) * gallons_to_Liters * liters_to_cubic_meters
        Infeasible_Routes_Emissions            = kg_to_Megaton * JetA_GHG * Infeasible_Routes_fuel_volume * density_JetA
        
        # Emissions of SAF-Scenario from SAF Flights 
        Flight_at_SAF_Airports_Using_SAF_Mo = Flight_at_SAF_Airports_Using_SAF.loc[Flight_at_SAF_Airports_Using_SAF['Month'] == m_i+1 ]    
        total_SAF_fuel_volume_required_mo   = np.sum(np.array(Flight_at_SAF_Airports_Using_SAF_Mo['Total Fuel Per Route (Gal)']))  
        SAF_volumes_mo                      = cumulative_fuel_use*total_SAF_fuel_volume_required_mo         
        SAF_Emissions                   = g_to_kg * kg_to_Megaton * np.sum(fuels_used['LCEF (gCO2e/MJ)']*fuels_used['Volumetric Energy Density (MJ/L)']*SAF_volumes_mo*gallons_to_Liters)
        Emissions_w_SAF_Aircraft[m_i]   = SAF_Emissions +  Infeasible_Routes_Emissions

        # Compute emissions without SAF integration        
        Conventional_Air_Travel_fuel_volume   = np.sum(np.array(Routes_and_Temp_Mo['Total Fuel Per Route (Gal)'])) * gallons_to_Liters * liters_to_cubic_meters
        Emissions_wo_SAF_Aircraft[m_i]        = kg_to_Megaton * JetA_GHG * Conventional_Air_Travel_fuel_volume * density_JetA

        # COST PER SEAT MILE 
        # CASM for normal operations without SAF aircraft  
        if len(Non_SAF_Flights_Mo ) == 0:
            pass
        else: 
            ASM_jet_A                = np.sum(Routes_and_Temp_Mo['Distance (miles)'] * Routes_and_Temp_Mo['Passengers'])
            Total_Fuel_Cost_jet_A    = np.sum(Routes_and_Temp_Mo['Fuel Cost']) 
            CASM_wo_SAF_Aircraft[m_i]  = 100*Total_Fuel_Cost_jet_A/ASM_jet_A    
       
        if len(Flight_at_SAF_Airports_Using_SAF_Mo)  == 0:
            pass
        else:    
            ASM_SAF                    = np.sum(Flight_at_SAF_Airports_Using_SAF_Mo['Distance (miles)'] * Flight_at_SAF_Airports_Using_SAF_Mo['Passengers']) 
            Total_Fuel_Cost_SAF        = np.sum(Flight_at_SAF_Airports_Using_SAF_Mo['Total Fuel Per Route (Gal)'] ) * SAF_dollars_per_gal 
            CASM_w_SAF_Aircraft[m_i]   = 100*Total_Fuel_Cost_SAF/ASM_SAF   
        
    return


def generate_eletrification_results_data(Flight_Ops):
     
     
    aircraft              = "Short-Haul (120 Pax)"
    airline               = 'All'
    battery_choice        = Commercial_Batteries['Battery Name'][13] 
    weight_fraction       = 35
    system_voltage        = 800 
    propulsive_efficiency = 90
    percent_adoption      = 100 
    month_no              = 1  
    cost_of_electricity   = 0.3 
        
    cell_V,
    capacity
    cell_C_max
    cell_e0
    cell_Temp    
 
        
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
    Routes_and_Temp_Mo  = Flight_Ops[Flight_Ops['Month'] == month_no+1 ]    
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
        Routes_and_Temp_Mo  = Flight_Ops.loc[Flight_Ops['Month'] == m_i+1 ] 
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
            
    
    return 
     
def generate_hydrogen_results_data(Hydrogen,Flight_Ops):
     

    H2_1                  = Hydrogen['H2 Fuel Name'][2] 
    H2_2                  = Hydrogen['H2 Fuel Name'][1] 
    H2_3                  = Hydrogen['H2 Fuel Name'][5] 
    H2_4                  = Hydrogen['H2 Fuel Name'][10]     
    selected_h2           = [H2_1,H2_2,H2_3,H2_4]      
    h2_vol_percentage     = 15  
    h2_airports           = " Top 10 Airports" # DONE  # [" Top 5 Airports"," Top 10 Airports"," Top 20 Airports", " Top 50 Airports",  "All Airports"]
    h2_percent_adoption   = 100    
    h2_cruise_alt         = 35000  # ft 
    mean_SFC_Imperial     = 0.6   
    H2_dollars_per_gal    = 8   
       
    #================================================================================================================================================  
    # Unit Conversions 
    #================================================================================================================================================
    coesa76_geom           = coesa76([ (h2_cruise_alt / Units.feet) /1000 ]) # conversion to km
    rho                    = Density_given_height(h2_cruise_alt)  # coesa76_geom.rho
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
          
    
          
        return  

 
    
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
    axes.grid(which='major', linestyle='-', linewidth=0.5, color='grey')
    axes.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')
    axes.grid(True)
    axes.get_yaxis().get_major_formatter().set_scientific(False)
    axes.get_yaxis().get_major_formatter().set_useOffset(False)

    return
 
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
 
     
if __name__ == '__main__':
    main()
    