
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTS 
# --------------------------------------------------------------------------------------------------------------------------------------------------- 

from RCAIDE.Framework.Core import  Data 
from RCAIDE.load import load as load_results
from RCAIDE.save import save as save_results 

import os 
import numpy as np           
import pandas as pd 
import json
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

  
    saf_data       =  generate_saf_results_data(Commercial_SAF,Flight_Ops,feedstocks) 
    save_results(saf_data, 'saf_data.res')
    
    electric_data  = generate_eletrification_results_data(Flight_Ops)
    save_results(electric_data, 'electric_data.res')
    
    hydrogen_data  = generate_hydrogen_results_data(Hydrogen,Flight_Ops) 
    save_results(hydrogen_data, 'hydrogen_data.res')    
     
    return  


def generate_saf_results_data(Commercial_SAF,Flight_Ops,feedstocks): 


    states_1  = ["Alabama","Arizona","Arkansas","California","Colorado", "Connecticut","Delaware","Florida"]  
    states_2  = ["Georgia","Idaho","Illinois", "Indiana","Iowa","Kansas","Kentucky","Louisiana"] 
    states_3  = ["Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana"] 
    states_4  = ["Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York", "North Carolina","Ohio"] 
    states_5  = ["Oklahoma","Oregon","Pennsylvania", "Rhode Island","South Carolina","South Dakota","North Dakota","Tennessee"]
    states_6  = ["Texas","Utah",  "Vermont","Virginia","Washington","West Virginia","Wisconsin","Wyoming"]
     
    feedstock_producing_states      = states_1 + states_2 + states_3 + states_4 +  states_5 +  states_6 
    
    selected_fuels_list                  = [[Commercial_SAF['Fuel Name'][4]],
                                            [Commercial_SAF['Fuel Name'][5]],
                                            [Commercial_SAF['Fuel Name'][10]],
                                            [Commercial_SAF['Fuel Name'][11]],
                                            [Commercial_SAF['Fuel Name'][19]],
                                            [Commercial_SAF['Fuel Name'][20]] ]   #[ HEFA-Soyabean-UCO , HEFA-Canola-UCO  HEFA-Soyabean-Vegetable oil , HEFA-Canola-Vegetable oil ,  ATJ,Corn-Corn ,  ATJ-Wheat-Wheat ]
    selected_feedstock_list         = ['Soybean', 'Canola', 'Soybean', 'Canola', 'Corn', 'Wheat']    
    selected_airpots_list           = [" Top 5 Airports"," Top 10 Airports", " Top 20 Airports"," Top 50 Airports"]
    percent_adoption_list           = np.linspace(0, 100, 11)   
    SAF_dollars_per_gal  =  8
    land_area                     = np.zeros((len(selected_fuels_list), len(selected_airpots_list), len(percent_adoption_list)))
    CASM_wo_SAF_Aircraft          = np.zeros((len(selected_fuels_list), len(selected_airpots_list), len(percent_adoption_list),12))
    CASM_w_SAF_Aircraft           = np.zeros((len(selected_fuels_list), len(selected_airpots_list), len(percent_adoption_list),12))
    
    for saf_i1 in  range(len(selected_fuels_list)):
        for saf_i2 in  range(len(selected_airpots_list)):
            for saf_i3 in  range(len(percent_adoption_list)): 
                selected_fuels     =  selected_fuels_list[saf_i1]
                selected_feedstock =  selected_feedstock_list[saf_i1]
                selected_airpots   =  selected_airpots_list[saf_i2]
                percent_adoption   =  percent_adoption_list[saf_i3]
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
                #fuel_percentages_list += percent_fuel_use # only saf 
                fuel_percentages_list += [100]
                fuels_percentages     = np.diff(np.array(fuel_percentages_list))/100  
                
                # Determine the percentage of neat (pure) saf and Jet-A1 using blending ratios   
                if Commercial_SAF['Fuel Name'][0] not in selected_fuels:  
                    selected_fuels   = [Commercial_SAF['Fuel Name'][0]] + selected_fuels
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
                
                land_area[saf_i1,saf_i2,saf_i3] = RCA
                
                while total_vol<RCA: 
                    total_vol += Used_Feedstock.loc[Used_Feedstock.index[idx]]['Acres Harvested']
                    Used_Feedstock["Feedstock Usage"][Used_Feedstock.index[idx]] = 0.1 
                    idx += 1      
                    if available_tracts == idx:
                        total_vol = 1E9
                        
                # Determine Cost per Seat Mile and Emissions  
                gallons_to_Liters             = 3.78541
            
                Emissions_w_SAF_Aircraft      = np.zeros(12)
                Emissions_wo_SAF_Aircraft     = np.zeros(12)
                
                for m_i in range(12): 
                    Routes_and_Temp_Mo                  = Flight_Ops[Flight_Ops['Month'] == m_i+1 ]   
                        
                    # EMISSIONS 
                    # Emissions of SAF-Scenario from non-SAF Flights 
                    Non_SAF_Flights_Mo                  = Non_SAF_Flights.loc[Non_SAF_Flights['Month'] == m_i+1 ]   
                    Infeasible_Routes_fuel_volume          = np.sum(np.array(Non_SAF_Flights_Mo['Total Fuel Per Route (Gal)'])) * gallons_to_Liters * liters_to_cubic_meters
                    Infeasible_Routes_Emissions            = kg_to_Megaton * JetA_GHG * Infeasible_Routes_fuel_volume * density_JetA
                    
                    # Emissions of SAF-Scenario from SAF Flights 
                    Flight_at_SAF_Airports_Using_SAF_Mo = Flight_at_SAF_Airports_Using_SAF.loc[Flight_at_SAF_Airports_Using_SAF['Month'] == m_i+1 ]
                    
                    # COST PER SEAT MILE 
                    # CASM for normal operations without SAF aircraft  
                    if len(Non_SAF_Flights_Mo ) == 0:
                        pass
                    else: 
                        ASM_jet_A                = np.sum(Routes_and_Temp_Mo['Distance (miles)'] * Routes_and_Temp_Mo['Passengers'])
                        Total_Fuel_Cost_jet_A    = np.sum(Routes_and_Temp_Mo['Fuel Cost']) 
                        CASM_wo_SAF_Aircraft[saf_i1,saf_i2,saf_i3,m_i]  = 100*Total_Fuel_Cost_jet_A/ASM_jet_A    
                   
                    if len(Flight_at_SAF_Airports_Using_SAF_Mo)  == 0:
                        pass
                    else:    
                        ASM_SAF                    = np.sum(Flight_at_SAF_Airports_Using_SAF_Mo['Distance (miles)'] * Flight_at_SAF_Airports_Using_SAF_Mo['Passengers']) 
                        Total_Fuel_Cost_SAF        = np.sum(Flight_at_SAF_Airports_Using_SAF_Mo['Total Fuel Per Route (Gal)'] ) * SAF_dollars_per_gal 
                        CASM_w_SAF_Aircraft[saf_i1,saf_i2,saf_i3,m_i]   = 100*Total_Fuel_Cost_SAF/ASM_SAF   
                                                    
    
    saf_data = Data()
    saf_data.saf_plot_1_selected_fuels_list       = selected_fuels_list        
    saf_data.saf_plot_1_selected_feedstock_list   = selected_feedstock_list     
    saf_data.saf_plot_1_selected_airpots_list     = selected_airpots_list        
    saf_data.saf_plot_1_percent_adoption_list     = percent_adoption_list     
    saf_data.saf_plot_1_land_area                 = land_area               
    saf_data.saf_plot_1_CASM_wo_SAF_Aircraft      = CASM_wo_SAF_Aircraft    
    saf_data.saf_plot_1_CASM_w_SAF_Aircraft       = CASM_w_SAF_Aircraft 
     
    states_1                       = ["Alabama","Arizona","Arkansas","California","Colorado", "Connecticut","Delaware","Florida"]  
    states_2                       = ["Georgia","Idaho","Illinois", "Indiana","Iowa","Kansas","Kentucky","Louisiana"] 
    states_3                       = ["Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana"] 
    states_4                       = ["Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York", "North Carolina","Ohio"] 
    states_5                       = ["Oklahoma","Oregon","Pennsylvania", "Rhode Island","South Carolina","South Dakota","North Dakota","Tennessee"]
    states_6                       = ["Texas","Utah",  "Vermont","Virginia","Washington","West Virginia","Wisconsin","Wyoming"]
     
    feedstock_producing_states     = states_1 + states_2 + states_3 + states_4 +  states_5 +  states_6 
    
    selected_fuels_list            = [[Commercial_SAF['Fuel Name'][0] , Commercial_SAF['Fuel Name'][10], Commercial_SAF['Fuel Name'][19],Commercial_SAF['Fuel Name'][21]]]
                      
                                   #[ Conventional Jet-A,  HEFA-Soyabean-Vegetable oil , ATJ,Corn-Corn , Gasification/FT Landfill Waste  ]
    selected_feedstock_list        = ['Soybean', 'Canola', 'Soybean', 'Canola', 'Corn', 'Wheat']    
    selected_airpots_list          = [" All Airports"]
    percent_adoption_list          = np.linspace(0, 100, 5)  
    percent_fuel_use_list          = [[25,50,75],[40,60,80],[20,40,60] ,[20,60,80]]
       
    Emissions_w_SAF_Aircraft      = np.zeros((len(selected_airpots_list), len(percent_adoption_list), len(percent_fuel_use_list),12))
    Emissions_wo_SAF_Aircraft     = np.zeros((len(selected_airpots_list), len(percent_adoption_list), len(percent_fuel_use_list),12))
    
    for saf_j1 in  range(len(selected_airpots_list)): 
        for saf_j2 in  range(len(percent_adoption_list)):
            for saf_j3 in  range(len(percent_fuel_use_list)):
                selected_fuels     =  selected_fuels_list[0]
                selected_feedstock =  selected_feedstock_list[0]
                selected_airpots   =  selected_airpots_list[saf_j1]
                percent_adoption   =  percent_adoption_list[saf_j2]
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
                fuel_percentages_list += percent_fuel_use_list[saf_j3] # only saf 
                fuel_percentages_list += [100]
                fuels_percentages     = np.diff(np.array(fuel_percentages_list))/100  
                
                # Determine the percentage of neat (pure) saf and Jet-A1 using blending ratios   
                if Commercial_SAF['Fuel Name'][0] not in selected_fuels:  
                    selected_fuels   = [Commercial_SAF['Fuel Name'][0]] + selected_fuels
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
                gallons_to_Liters             = 3.78541
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
                    Emissions_w_SAF_Aircraft[saf_j1,saf_j2,saf_j3,m_i]   = SAF_Emissions +  Infeasible_Routes_Emissions
            
                    # Compute emissions without SAF integration        
                    Conventional_Air_Travel_fuel_volume   = np.sum(np.array(Routes_and_Temp_Mo['Total Fuel Per Route (Gal)'])) * gallons_to_Liters * liters_to_cubic_meters
                    Emissions_wo_SAF_Aircraft[saf_j1,saf_j2,saf_j3,m_i]        = kg_to_Megaton * JetA_GHG * Conventional_Air_Travel_fuel_volume * density_JetA
                      
    saf_data.saf_plot_2_selected_feedstock_list    = selected_feedstock_list      
    saf_data.saf_plot_2_selected_airpots_list      = selected_airpots_list        
    saf_data.saf_plot_2_percent_adoption_list      = percent_adoption_list        
    saf_data.saf_plot_2_percent_fuel_use_list      = percent_fuel_use_list         
    saf_data.saf_plot_2_percent_fuel_use_list      = percent_fuel_use_list 
    saf_data.saf_plot_2_Emissions_w_SAF_Aircraft   = Emissions_w_SAF_Aircraft  
    saf_data.saf_plot_2_Emissions_wo_SAF_Aircraft  = Emissions_wo_SAF_Aircraft
    
    return saf_data


def generate_eletrification_results_data(Flight_Ops):
     
      
    weight_fraction   = np.linspace(5,60,12)
    cell_e0           = np.linspace(200,1000,9) 
    aircraft_capacity = np.array([19, 88, 120, 189, 368])
    aircraft_range    = np.zeros((len(aircraft_capacity),12,len(weight_fraction),len(cell_e0 ) )) 
    passenger_volume  = np.zeros((len(aircraft_capacity),12,len(weight_fraction),len(cell_e0 ) )) 
    CASM_electric     = np.zeros((len(aircraft_capacity),12,len(weight_fraction),len(cell_e0 ) )) 
    CASM_Jet_A        = np.zeros((len(aircraft_capacity),12,len(weight_fraction),len(cell_e0 ) )) 
    
    for i in  range(len(weight_fraction)):
        for j in range(len(cell_e0)):
            for ac in  range(len(aircraft_capacity)): 
                for m_i in  range(12):
                    system_voltage        = 800 
                    propulsive_efficiency = 90
                    percent_adoption      = 100 
                    cost_of_electricity   = 0.3 
                        
                    cell_V       = 4.2         
                    capacity_Wh  = 4000
                    capacity     = capacity_Wh /cell_V 
                    cell_C_max   = 6
                    cell_Temp    = [-20,50] 
                        
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
                    month               = months[m_i] 
                    Routes_and_Temp_Mo  = Flight_Ops[Flight_Ops['Month'] == m_i+1 ]    
                    V_cell              = cell_V
                    e_cell              = cell_e0[j] *Wh_per_kg_to_J
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
                    W_bat           = (weight_fraction[i]/100) * Takeoff_Weight
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
                    Routes_and_Temp_Mo['Aircraft_Battery_Energy']  = E_bat
                    
                    # compute range 
                    Range_mi = meters_to_miles * (e_cell/9.81) * Lift_to_Drag_ratio * (weight_fraction[i]/100)* eta_0
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
                    
                    Routes_and_Temp_Mo_1   = Routes_and_Temp_Mo[Routes_and_Temp_Mo['Estimated Aircraft Capacity'] == aircraft_capacity[ac] ] 
                     
                    # COST PER SEAT MILE
                    # CASM for normal operations without electric aircraft 
                    ASM_jet_A_1                = np.sum(Routes_and_Temp_Mo_1['Distance (miles)'] * Routes_and_Temp_Mo_1['Passengers']) 
                    Total_Fuel_Cost_jet_A_1    = np.sum(Routes_and_Temp_Mo_1['Fuel Cost'])   
                    CASM_wo_E_Aircraft_1       = 100*Total_Fuel_Cost_jet_A_1/ASM_jet_A_1     
                    
                    # CASM for when electric aircraft are integrated into fleet
    
                    Feasible_Routes_e1   = Feasible_Routes[Feasible_Routes['Estimated Aircraft Capacity'] ==  aircraft_capacity[ac]  ] 
                    
                    if len(Feasible_Routes_e1['E_Passengers']) == 0:  
                        CASM_w_E_Aircraft_1  = 0 
                    else:  
                        Electricity_Cost_1   = sum((np.array(Feasible_Routes_e1['Aircraft_Battery_Energy'])*joule_to_kWh*cost_of_electricity) * np.array(Feasible_Routes_e1['No of Flights Per Month'])) 
                        ASM_electric_1       = sum(Feasible_Routes_e1['Distance (miles)'] *Feasible_Routes_e1['E_Passengers']  )  
                        CASM_w_E_Aircraft_1 =  100*Electricity_Cost_1/ASM_electric_1  
     
                    # Figyre 1 Data 
                    aircraft_range[ac,m_i,i, j]  = Range_mi[est_pax_capacity==aircraft_capacity[ac] ][0] 
                    
                    # Figure 2 Data
                    passenger_volume[ac,m_i,i, j]  = sum(Feasible_Routes_e1['E_Passengers']  )   
                    # Figure 2 Data
                    CASM_electric[ac,m_i,i, j]  = CASM_w_E_Aircraft_1  
    
                    CASM_Jet_A[ac,m_i,i, j]  = CASM_wo_E_Aircraft_1
                    
                    
    electric_data = Data()            
    electric_data.weight_fraction   =  weight_fraction   
    electric_data.cell_e0           =  cell_e0          
    electric_data.aircraft_range    =  aircraft_range   
    electric_data.passenger_volume  =  passenger_volume 
    electric_data.CASM_electric     =  CASM_electric    
    electric_data.CASM_Jet_A        =  CASM_Jet_A       
    electric_data.aircraft_capacity =  aircraft_capacity
    return electric_data
     
def generate_hydrogen_results_data(Hydrogen,Flight_Ops):
    
    
    volume_fraction_list    = np.linspace(5, 60,12)
    aircraft_capacity       = np.array([19, 88, 120, 189, 368])
    aircraft_range          = np.zeros((len(aircraft_capacity),12,len(volume_fraction_list))) 
    passenger_volume        = np.zeros((len(aircraft_capacity),12,len(volume_fraction_list)))  
    CASM_wo_H2_Aircraft     = np.zeros((len(aircraft_capacity),12,len(volume_fraction_list)))
    CASM_w_H2_Aircraft      = np.zeros((len(aircraft_capacity),12,len(volume_fraction_list)))
    CO2e_w_H2_Aircraft      = np.zeros((len(aircraft_capacity),12,len(volume_fraction_list)))
    CO2e_wo_H2_Aircraft     = np.zeros((len(aircraft_capacity),12,len(volume_fraction_list)))
       
    for i in range(len(volume_fraction_list)):
        for ac in  range(len(aircraft_capacity)): 
            H2_1                  = Hydrogen['H2 Fuel Name'][2]    
            selected_h2           = [H2_1]      
            h2_vol_percentage     = volume_fraction_list[i]
            h2_airports           = " All Airports"
            h2_percent_adoption   = 100    
            h2_cruise_alt         = 35000  # ft 
            mean_SFC_Imperial     = 0.6   
            H2_dollars_per_gal    = 8   
               
            #================================================================================================================================================  
            # Unit Conversions 
            #================================================================================================================================================ 
            r_e=3959*5280 #miles to feet
            R=1716.5 #ft2/R-sec
            g0=32.174 #ft/s**2
            T0=518.69 #deg R     
            h=(r_e/(r_e+h2_cruise_alt))*h2_cruise_alt
             
            if h  <= 36000:
                h0=0
                T0=518.69 #deg R 
                rho0=2.3769e-3 #slugs/ft3 
                a1=-3.57/1000 #deg R/ft
                T = T0 + a1*(h -h0) 
                rho_Imperial  = rho0*(T/T0)**(-((g0/(R*a1))+1))
            else:
                h0=36000
                rho0=7.103559955456123e-04 #from running code at 36000
                T=389.99 
                rho_Imperial = rho0*np.exp((-g0/(R*T))*(h-h0))  
             
            rho =  rho_Imperial *  515.378818 # conversion to kg / m3 
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
            additional_volume_fractionl_required =  unused_volume_of_aircraft -  np.ones_like(percentage_capacity_of_aircraft)*volume_fraction
            
            # If required volume is negative, we do not need to remove passengers 
            additional_volume_fractionl_required[additional_volume_fractionl_required>0] =  0
            
            # If additional_volume_fractionl_required is negative, we need to take off more passenger to accomodate H2 volume:
            # Determine how many additional passengers must be removed from aircraft to meet specified volume fraction
            Flight_Ops['H2_Passengers']  =   np.array(Flight_Ops['Passengers']) -  (1 + additional_volume_fractionl_required) * np.array(Flight_Ops['Estimated Aircraft Capacity']) 
               
            # Compute the percentages of different types of fuels used 
            percentage_H2_process_vals  = [0] 
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
            for m_i in range(12):
    
                 
                # Filter out aircraft 
                H2_Flights_AC                = H2_Flights[H2_Flights['Estimated Aircraft Capacity'] == aircraft_capacity[ac]]                
                H2_Flights_Mo                = H2_Flights_AC.loc[H2_Flights_AC['Month'] == m_i+1] 
                Non_H2_Flights_Mo            = Non_H2_Flights.loc[Non_H2_Flights['Month'] == m_i+1 ]   
                Non_H2_Flights_JetA_vol_mo   = np.sum(np.array(Non_H2_Flights_Mo['Total Fuel Per Route (Gal)']))   * gallons_to_Liters * liters_to_cubic_meters   
                H2_volumes_mo                = np.sum(np.array(H2_Flights_Mo['H2_volume']))
                
                if len(Non_H2_Flights_Mo ) == 0:
                    pass
                else:
                    ASM_JetA               = np.sum(Non_H2_Flights_Mo['Distance (miles)'] * Non_H2_Flights_Mo['Passengers'])
                    Total_Fuel_Cost_JetA   = np.sum(Non_H2_Flights_Mo['Fuel Cost']) 
                    CASM_wo_H2_Aircraft[ac,m_i,i]         = 100*Total_Fuel_Cost_JetA/ASM_JetA # convert to cents/pass-mile
               
                if len(H2_Flights_Mo)  == 0:
                    pass
                else:    
                    ASM_H2               = np.sum(H2_Flights_Mo['Distance (miles)'] * H2_Flights_Mo['H2_Passengers']) 
                    Total_Fuel_Cost_H2   = H2_volumes_mo /liters_to_cubic_meters /gallons_to_Liters  * H2_dollars_per_gal 
                    CASM_w_H2_Aircraft[ac,m_i,i]         = 100*Total_Fuel_Cost_H2/ASM_H2   # convert to cents/pass-mile
                
                # CO2e of the Hydrogen-powered flights in the scenario with flight operations with H2 tech
                CO2e_H2_aircraft_in_JetAH2_Fleet    = kg_to_Megaton * np.sum(np.array(H2_used['Direct GHG emissions [kg CO2e/kg H2]'])*percentage_H2_process)*H2_volumes_mo
         
                # CO2e of the JetA-powered flights in the scenario with flight operations with H2 tech        
                CO2e_JetA_aircraft_in_JetAH2_Fleet  = kg_to_Megaton * JetA_GHG * Non_H2_Flights_JetA_vol_mo *  density_JetA 
                
                # Total CO2e flights in the scenario with flight operations with H2 tech
                CO2e_w_H2_Aircraft[ac,m_i,i]       =  CO2e_H2_aircraft_in_JetAH2_Fleet +  CO2e_JetA_aircraft_in_JetAH2_Fleet
         
                # Total CO2e flights in the scenario with flight operations without H2 tech    
                Conventional_Air_Travel             = Flight_Ops.loc[Flight_Ops['Month'] == m_i+1]
                Conventional_Air_Travel_fuel_volume = np.sum(np.array(Conventional_Air_Travel['Total Fuel Per Route (Gal)'])) * gallons_to_Liters * liters_to_cubic_meters
                CO2e_wo_H2_Aircraft[ac,m_i,i]            = kg_to_Megaton * JetA_GHG * Conventional_Air_Travel_fuel_volume * density_JetA        
                      

                passenger_volume[ac,m_i,i]   = np.sum(H2_Flights_Mo['H2_Passengers'])

                aircraft_range[ac,m_i,i]     = np.maximum(0, Range_mi[original_pax_capacity==aircraft_capacity[ac] ][0])
                
            
    hydrogen_data = Data()     
    hydrogen_data.volume_fraction         = volume_fraction_list     
    hydrogen_data.aircraft_capacity       = aircraft_capacity   
    hydrogen_data.aircraft_range          = np.nan_to_num( aircraft_range )
    hydrogen_data.passenger_volume        = passenger_volume    
    hydrogen_data.CASM_wo_H2_Aircraft     = CASM_wo_H2_Aircraft 
    hydrogen_data.CASM_w_H2_Aircraft      = CASM_w_H2_Aircraft  
    hydrogen_data.CO2e_w_H2_Aircraft      = CO2e_w_H2_Aircraft  
    hydrogen_data.CO2e_wo_H2_Aircraft     = CO2e_wo_H2_Aircraft         
    return  hydrogen_data 

 
    

 
     
if __name__ == '__main__':
    main()
    