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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.colors as mcolors
 
def main():
    
    # ---------------------------------------------------------------------------------------------------------------------------------------------------
    # Data
    # ---------------------------------------------------------------------------------------------------------------------------------------------------
    
    separator            = os.path.sep
    technology_filename  = 'Data' + separator + 'Technology' + separator + 'Technology_Data.xlsx'
    SAT_data             = pd.read_excel(technology_filename,sheet_name=['Commercial_Batteries','Battery_Development','Electric_Motor_Development']) 
    Commercial_Batteries = SAT_data['Commercial_Batteries'] 
    Electric_Motor_Development = SAT_data['Electric_Motor_Development']
    a                    = Commercial_Batteries['Brand']  
    b                    = Commercial_Batteries['Chemistry']
    c                    = Commercial_Batteries['Model']
    Commercial_Batteries["Battery Name"] = a + ': ' + b + '-' + c  
    Battery_Development  = SAT_data['Battery_Development']   
    route_temp_filename  = 'Data' + separator +  'Air_Travel' + separator + 'American_Airlines_Flight_Ops_and_Climate.xlsx'
    Routes_and_Temp      = pd.read_excel(route_temp_filename,sheet_name=['Sheet1']) 
    Routes_and_Temp      = Routes_and_Temp['Sheet1']     
    filename             = 'Data' + separator  +  'US_Climate' + separator +  'Monthly_US_County_Temperature_2019.xlsx'
    Temeperature_data    = pd.read_excel(filename,sheet_name=['US_County_Temperature_F'])  
    US_Temperature_F     = Temeperature_data['US_County_Temperature_F'] 
    
    
    selected_brand     = 'All'
    selected_chemistry = 'All'
    selected_x_axis    = list(Commercial_Batteries.columns.values)[4:18][7]
    selected_y_axis    = list(Commercial_Batteries.columns.values)[4:18][9]
    switch_off         = False      
    
    bat_1           = Commercial_Batteries['Battery Name'][5] 
    bat_2           = Commercial_Batteries['Battery Name'][13] 
    bat_3           = Commercial_Batteries['Battery Name'][6] 
    switch_off         = False          

    selected_sector = 'All'
    selected_type   = 'All'
    switch_off         = False  

    switch_off            = False     
    month_no              = 1     
  
    aircraft              = 'Boeing 737 MAX-8'
    battery_choice        = Commercial_Batteries['Battery Name'][13] 
    weight_fraction       = 35
    system_voltage        = 400 
    propulsive_efficiency = 90
    percent_adoption      = 100 
    month_no              = 1 
    switch_off            = False 
    cost_of_electricity   = 0.3 # per Kw 
    fig = generate_electric_flight_operations_plots(Routes_and_Temp,Commercial_Batteries,aircraft,battery_choice,weight_fraction,system_voltage,propulsive_efficiency,percent_adoption,month_no,cost_of_electricity,switch_off)
                
    plt.show()
    
    return 
    
def generate_electric_flight_operations_plots(Routes_and_Temp,Commercial_Batteries,aircraft,battery_choice,weight_fraction,system_voltage,propulsive_efficiency,percent_adoption,month_no,cost_of_electricity ,switch_off):    
    
    #if aircraft == 'ATR 72-600': 
        #P_max           = 1568083 
        #W_0             = 23000    
        #fuel_volume     = 6410.2
        #fuel_economy    = 3.27 
        #passengers      = 72 
        #CL_cruise       =  
        #fuselage_volume =  
        #CD_cruise       =  
        #S_ref           = 61.0   
        #L_div_D         = CL_cruise/CD_cruise 
        
    if aircraft == 'Embraer 190': 
        P_max           = 15.27 * 1E6
        W_0             = 52290    
        fuel_volume     = 16629 
        fuel_economy    = 3.54 
        passengers      = 106
        fuselage_volume = 75.56964 
        SFC_lb_lbfhr    = 0.624 
        CL_cruise       = 0.501
        CD_cruise       = 0.0342 
        S_ref           = 92.53 
        L_div_D         = CL_cruise/CD_cruise
        
    elif aircraft == "Boeing 737 MAX-8":
        P_max           = 15000000 
        W_0             = 79015.8    
        fuel_volume     = 26024 
        fuel_economy    = 2.04 
        passengers      = 162  
        SFC_lb_lbfhr    = 0.6262 # (lb / lbf - hr) , in cruise 
        fuselage_volume = 129.742 # fuselage cabin length * fuselage radius^2 * pi * 0.5 (where passengers are)
        CD_cruise       = 0.027
        CL_cruise       = 0.569
        S_ref           = 124.862
        L_div_D         = CL_cruise/CD_cruise
        
    elif aircraft == 'Airbus A320 neo': 
        P_max           = 15000000 
        W_0             = 78000 
        fuel_volume     = 26024  
        fuel_economy    = 2.04  
        passengers      = 162
        W_0             = 79015.8  
        SFC_lb_lbfhr    = 0.6262 # (lb / lbf - hr), in cruise 
        fuselage_volume = 129.742 # fuselage cabin length * fuselage radius^2 * pi * 0.5 (where passengers are)
        CL_cruise       = 0.55 
        CD_cruise       = 0.035
        S_ref           = 127   
        L_div_D         = CL_cruise/CD_cruise
        
    #elif aircraft == "Boeing 787-8": 
        #P_max           =   
        #fuel_volume     = 126206  
        #fuel_economy    = 2.68    
        #passengers      = 238
        #SFC_lb_lbfhr    = 
        #W_0             = 227900   
        #CL_cruise       = 0.055 
        #fuselage_volume = 480.8650
        #CD_cruise       = 0.035
        #S_ref           = 377
        #L_div_D         = CL_cruise/CD_cruise
        

    elif aircraft == "Boeing 777-300": 
        P_max           = 68 * 1E6
        fuel_volume     = 126206  
        fuel_economy    = 2.68    
        passengers      = 238
        SFC_lb_lbfhr    = 0.51 
        W_0             = 227900   
        CL_cruise       = 0.55 
        fuselage_volume = 682.008 
        CD_cruise       = 0.035
        S_ref           = 377
        L_div_D         = CL_cruise/CD_cruise
        
    #elif aircraft == 'Airbus A350-1000':
        #P_max          =   
        #fuel_volume    = 166488        # liters 
        #fuel_economy   = 2.85          # L/100km per Pax
        #SFC_lb_lbfhr   =  
        #passengers     = 327  
        #W_0            = 322050 
        #fuselage_volume = 
        #CL_cruise      = 0.055 
        #CD_cruise      = 0.035
        #S_ref          = 464.3
        #L_div_D        = CL_cruise/CD_cruise       
             
    CO2e_per_mile      = 9.0736                 
    Wh_per_kg_to_J     = 3600.0
    Ah_to_C            = 3600.0  
    V_bat              = system_voltage
    eta_0              = propulsive_efficiency/100 
    
    #================================================================================================================================================  
    # Compute Feasible Routes 
    #================================================================================================================================================    
    # Compute Range  
    months          = ['January', 'February', 'March', 'April', 'May', 'June', 'July','August', 'September', 'October', 'November', 'December']      
    month           =  months[month_no] 
    data            = Commercial_Batteries[Commercial_Batteries['Battery Name'] == battery_choice] 
    V_cell          = np.array(data['Nominal Voltage (V)'])[0]
    #e_cell          = np.array(data['Gravimetric Energy Density (Wh/kg)'])[0] *Wh_per_kg_to_J
    e_cell          = np.arange(200, 401, 10) *Wh_per_kg_to_J
    q_cell          = np.array(data['Capacity (mAh)'])[0]/1000 * Ah_to_C  # conversion to coulombs
    i_max           = np.array(data['Maximum Discharge Current (A)'])[0] # amps   
    Min_Temp        = np.array(data['Minimum Disharge Temperature (°C)'])[0]*9/5 + 32
    Max_Temp        = np.array(data['Maximum Discharge Temperature (°C)'])[0]*9/5 + 32  
    I_bat           = P_max/ V_bat
    n_series        = V_bat/V_cell
    wf              = np.arange(0.1, 0.8, 0.1)
    Routes_and_Temp_Mo                    = Routes_and_Temp[Routes_and_Temp['Month'] == month_no+1 ]  
    W_bat           = np.zeros(len(wf))
    W_residual      = np.zeros((len(wf), len(Routes_and_Temp_Mo)))
    passenger_reductions = np.zeros((len(wf), len(Routes_and_Temp_Mo)))
    remaining_Pax   = np.zeros((len(wf), len(Routes_and_Temp_Mo)))
    E_bat           = np.zeros((len(wf), len(e_cell)))
    Q_bat           = np.zeros((len(wf), len(e_cell)))
    n_parallel      = np.zeros((len(wf), len(e_cell)))
    n_parallel_min  = np.zeros((len(wf), len(e_cell)))
    Range           = np.zeros((len(wf), len(e_cell)))
    Range_mi        = np.zeros((len(wf), len(e_cell)))
    v_remaining_Pax = np.zeros(len(wf))
    
    for i in range(len(wf)):
        W_bat[i]           = (wf[i]) * W_0
        for j in range(len(e_cell)):
            E_bat[i, j]    = W_bat[i] * e_cell[j]
            Q_bat[i, j]    = E_bat[i, j] / V_bat
            n_parallel[i, j] = Q_bat[i, j] / q_cell 
            n_parallel_min[i, j] = I_bat / i_max 
        
            if n_parallel_min[i, j] < n_parallel[i, j]: 
                Range[i, j] = (e_cell[j]/9.81) * L_div_D * (wf[i]) * eta_0
            else:  
                Range[i, j] = 0 
            Range_mi[i, j] = Range[i, j] * 0.000621371
    
            # Compute distances between departure and destimation points 
       
        
            Jet_A_density                         = 800.0 #  kg/m3
            fuel_volume_L                         = Routes_and_Temp_Mo['Fuel Consumed Per Flight (Liters)']  
            fuel_volume_m_3                       = fuel_volume_L*0.001  
            W_f                                   = Jet_A_density*fuel_volume_m_3
            #for k in range(len(W_f)):
            W_residual[i]                         = W_bat[i]-W_f  
            weight_per_pass                       = 158.757 # in kg  (250 lb for person, 100 lb for luggage) 
            passenger_reductions[i]               = np.ceil(np.array(W_residual[i])/weight_per_pass)   
            passenger_reductions[passenger_reductions < 0] = 0
            original_Pax_volume                   = np.array(Routes_and_Temp_Mo['Passengers'])
            remaining_Pax[i]                      = original_Pax_volume - passenger_reductions[i]*np.array(Routes_and_Temp_Mo['No of Flights Per Month'])
            remaining_Pax[remaining_Pax<0]        = 0
            v_remaining_Pax[i]                    = np.sum(remaining_Pax[i])
            
    #================================================================================================================================================  
    # Plots
    #================================================================================================================================================          

    # Range VS Spec En dens VS battery weight fraction

    fig_4 = plt.figure()

    plt.plot(e_cell/Wh_per_kg_to_J, Range_mi[0,:], '.-', label='Battery weight fraction: 0.1', color='C0')
    plt.plot(e_cell/Wh_per_kg_to_J, Range_mi[1,:], '.-', label='Battery weight fraction: 0.2', color='C1')
    plt.plot(e_cell/Wh_per_kg_to_J, Range_mi[2,:], '.-', label='Battery weight fraction: 0.3', color='C2')
    plt.plot(e_cell/Wh_per_kg_to_J, Range_mi[3,:], '.-', label='Battery weight fraction: 0.4', color='C3') 
    plt.plot(e_cell/Wh_per_kg_to_J, Range_mi[4,:], '.-', label='Battery weight fraction: 0.5', color='C4') 
    plt.plot(e_cell/Wh_per_kg_to_J, Range_mi[5,:], '.-', label='Battery weight fraction: 0.6', color='C5') 
    plt.plot(e_cell/Wh_per_kg_to_J, Range_mi[6,:], '.-', label='Battery weight fraction: 0.7', color='C6')     
    plt.xlabel('Specific energy density [Wh/kg]')
    plt.ylabel('Range [mi]')
    plt.title('Range vs. specific energy density')
    plt.legend()
    plt.grid(True) 

    # -------------------------------------------------------------------------------------------------------

    # Create figure and 3D axes
    fig_5 = plt.figure()
    #ax = fig_5.add_subplot(111, projection='3d')

    ## Create a color map and norm for the colorbar
    #cmap = cm.get_cmap('coolwarm_r')  # 'coolwarm' is a red to blue colormap
    #norm = mcolors.Normalize(vmin=200, vmax=400)  # Normalize between 200 and 400 Wh/kg    

    ## Create rings
    #for i in range(len(wf)):
        #theta = np.linspace(0, 2 * np.pi, 100)
        #r = np.linspace(min(Range_mi[i,:]), max(Range_mi[i,:]), 100)
        #r_mesh, theta_mesh = np.meshgrid(r, theta)

        #X = r_mesh * np.cos(theta_mesh)
        #Y = r_mesh * np.sin(theta_mesh)
        #Z = np.ones_like(X) * wf[i]  # Align along the Battery Weight Fraction axis

        ## Color gradient from blue to red using the inverted color map
        #colors = cmap(np.linspace(0, 1, 100))
        #colors = np.repeat(colors[np.newaxis, :, :], 100, axis=0)  # Repeat along the first axis
        #colors[:, :, 3] = 0.1  # Set alpha channel for transparency

        #ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1, linewidth=0, antialiased=False)

    ## Set limits and labels
    #ax.set_xlim([-1500, 1500])
    #ax.set_ylim([-1500, 1500])
    #ax.set_zlim([0, 1])
    #ax.set_xlabel('Range [mi]')
    #ax.set_ylabel('Range [mi]')
    #ax.set_zlabel('Battery Weight Fraction')    

    fig_6 = plt.figure()

    plt.plot(wf, v_remaining_Pax, '.-', color='C0')
    plt.xlabel('Battery weight fraction [-]')
    plt.ylabel('No of passengers [-]')
    plt.title('No of passengers per month vs. Battery weight fraction')
    plt.legend()
    plt.grid(True) 

    #================================================================================================================================================  
    # Compute Feasible Routes 
    #================================================================================================================================================    
    # Compute Range  
    months          = ['January', 'February', 'March', 'April', 'May', 'June', 'July','August', 'September', 'October', 'November', 'December']      
    month           =  months[month_no] 
    data            = Commercial_Batteries[Commercial_Batteries['Battery Name'] == battery_choice] 
    V_cell          = np.array(data['Nominal Voltage (V)'])[0]
    #e_cell          = np.array(data['Gravimetric Energy Density (Wh/kg)'])[0] *Wh_per_kg_to_J
    e_cell          = np.arange(300, 1401, 100) *Wh_per_kg_to_J
    q_cell          = np.array(data['Capacity (mAh)'])[0]/1000 * Ah_to_C  # conversion to coulombs
    i_max           = np.array(data['Maximum Discharge Current (A)'])[0] # amps   
    Min_Temp        = np.array(data['Minimum Disharge Temperature (°C)'])[0]*9/5 + 32
    Max_Temp        = np.array(data['Maximum Discharge Temperature (°C)'])[0]*9/5 + 32  
    I_bat           = P_max/ V_bat
    n_series        = V_bat/V_cell
    Routes_and_Temp_Mo                    = Routes_and_Temp[Routes_and_Temp['Month'] == month_no+1 ]  
    W_bat           = np.zeros(len(e_cell))
    wf           = np.zeros(len(e_cell))
    W_residual      = np.zeros((len(e_cell), len(Routes_and_Temp_Mo)))
    passenger_reductions = np.zeros((len(e_cell), len(Routes_and_Temp_Mo)))
    remaining_Pax   = np.zeros((len(e_cell), len(Routes_and_Temp_Mo)))
    Range           = np.zeros(len(e_cell))
    Range_mi        = np.zeros(len(e_cell))
    v_remaining_Pax = np.zeros(len(e_cell))
    Infeasible_Passenger_Miles = np.zeros(len(e_cell))
    Total_Passenger_Miles = np.zeros(len(e_cell)) 
    battery = np.zeros(len(e_cell))               
    no_battery = np.zeros(len(e_cell))         
    ASM_jet_A = np.zeros(len(e_cell))          
    Total_Fuel_Cost_jet_A = np.zeros(len(e_cell)) 
    CASM_jet_A = np.zeros(len(e_cell)) 
    Total_Fuel_Cost_electric = np.zeros(len(e_cell))  
    ASM_electric = np.zeros(len(e_cell))              
    CASM_electric = np.zeros(len(e_cell))          
    E_bat           = 23396578380.0
    
    for j in range(len(e_cell)):
        W_bat[j]    =  E_bat/e_cell[j]
        wf[j]       = W_bat[j]/ W_0
        Q_bat    = E_bat / V_bat
        n_parallel = Q_bat / q_cell 
        n_parallel_min = I_bat / i_max 
    
        if n_parallel_min < n_parallel: 
            Range[j] = (e_cell[j]/9.81) * L_div_D * (wf[j]) * eta_0
        else:  
            Range[j] = 0 
        Range_mi[j] = Range[j] * 0.000621371

        # Compute distances between departure and destimation points 
   
        Jet_A_density                         = 800.0 #  kg/m3
        fuel_volume_L                         = Routes_and_Temp_Mo['Fuel Consumed Per Flight (Liters)']  
        fuel_volume_m_3                       = fuel_volume_L*0.001  
        W_f                                   = Jet_A_density*fuel_volume_m_3
        #for k in range(len(W_f)):
        W_residual[j]                         = W_bat[j]-W_f  
        weight_per_pass                       = 158.757 # in kg  (250 lb for person, 100 lb for luggage) 
        passenger_reductions[j]               = np.ceil(np.array(W_residual[j])/weight_per_pass)   
        passenger_reductions[passenger_reductions < 0] = 0
        original_Pax_volume                   = np.array(Routes_and_Temp_Mo['Passengers'])
        remaining_Pax[j]                      = original_Pax_volume - passenger_reductions[j]*np.array(Routes_and_Temp_Mo['No of Flights Per Month'])
        remaining_Pax[remaining_Pax<0]        = 0
        v_remaining_Pax[j]                    = np.sum(remaining_Pax[j])
        Routes_and_Temp_Mo['E_Passengers']   = remaining_Pax[j] 
        
        Feasible_Routes_0    = Routes_and_Temp_Mo[Routes_and_Temp_Mo['E_Passengers'] > 0 ] 
        Infeasible_Routes_0  = Routes_and_Temp_Mo[Routes_and_Temp_Mo['E_Passengers'] < 0 ]   
        Feasible_Routes_1    = Feasible_Routes_0[Feasible_Routes_0['Distance (miles)'] < Range_mi[j] ] 
        Infeasible_Routes_1  = Feasible_Routes_0[Feasible_Routes_0['Distance (miles)'] > Range_mi[j] ]  
        Feasible_Routes_2    = Feasible_Routes_1[Feasible_Routes_1['Origin ' + month] > Min_Temp] 
        Infeasible_Routes_2  = Feasible_Routes_1[Feasible_Routes_1['Origin ' + month] < Min_Temp] 
        Feasible_Routes_3    = Feasible_Routes_2[Feasible_Routes_2['Origin ' + month] < Max_Temp] 
        Infeasible_Routes_3  = Feasible_Routes_2[Feasible_Routes_2['Origin ' + month] > Max_Temp]  
        Feasible_Routes_4    = Feasible_Routes_3[Feasible_Routes_3['Destination ' + month] > Min_Temp] 
        Infeasible_Routes_4  = Feasible_Routes_3[Feasible_Routes_3['Destination ' + month] < Min_Temp] 
        Feasible_Routes_5    = Feasible_Routes_4[Feasible_Routes_4['Destination ' + month] < Max_Temp] 
        Infeasible_Routes_5  = Feasible_Routes_4[Feasible_Routes_4['Destination ' + month] > Max_Temp] 
        Feasible_Routes      = Feasible_Routes_5.head(int(len(Feasible_Routes_5)*percent_adoption/100 )) 
        Infeasible_Routes_6  = Feasible_Routes_5.tail(int(len(Feasible_Routes_5)*(100 - percent_adoption)/100 ))
        Infeasible_Routes    = pd.concat([Infeasible_Routes_0,Infeasible_Routes_1,Infeasible_Routes_2,Infeasible_Routes_3,Infeasible_Routes_4,Infeasible_Routes_5,Infeasible_Routes_6])     
    
    
        Infeasible_Passenger_Miles[j]  = np.sum(np.array(Infeasible_Routes[['Distance (miles)']]))
        Total_Passenger_Miles[j]       = np.sum( np.array(Feasible_Routes[['Distance (miles)']]))  
        battery[j]                 = Total_Passenger_Miles[j] * CO2e_per_mile
        no_battery[j]              = Infeasible_Passenger_Miles[j] * CO2e_per_mile 
    
        # Infeasible Routes (Fuel) Energy Carrier Cost Per Seat Mile 
        ASM_jet_A[j]             = np.sum(Infeasible_Routes['Distance (miles)'] * Infeasible_Routes['Passengers'])
        Total_Fuel_Cost_jet_A[j] = np.sum(Infeasible_Routes['Fuel Cost'])
        
        # Compute electric CASM
        CASM_jet_A[j]        = 100*Total_Fuel_Cost_jet_A[j]/ASM_jet_A[j]    
        
        # Compute electric CASM
        if len(Feasible_Routes['E_Passengers']) == 0:
            electric_flight_passengers = 0
            CASM_electric          = 0
        else: 
            electric_flight_passengers = Feasible_Routes['E_Passengers']  
            joule_to_kWh               = 2.77778e-7
            Total_Fuel_Cost_electric[j]   = sum( (E_bat*joule_to_kWh*cost_of_electricity) * np.array(Feasible_Routes['No of Flights Per Month']) )
            ASM_electric[j]               = sum(Feasible_Routes['Distance (miles)'] * electric_flight_passengers)
            CASM_electric[j]         = 100*Total_Fuel_Cost_electric[j]/ASM_electric[j]               
        
        
    fig_7 = plt.figure()

    plt.plot(e_cell/Wh_per_kg_to_J, v_remaining_Pax, '.-', color='C0')
    plt.xlabel('Specific energy density [Wh/kg]')
    plt.ylabel('No of passengers [-]')
    plt.title('No of passengers per month vs. Battery specific energy density')
    plt.legend()
    plt.grid(True)    
    
    fig_8 = plt.figure()

    plt.plot(e_cell/Wh_per_kg_to_J, CASM_electric, '.-', label='CASM electric', color='C0')
    plt.plot(e_cell/Wh_per_kg_to_J, CASM_jet_A, '.-', label='CASM Jet A', color='C1')
    plt.xlabel('Specific energy density [Wh/kg]')
    plt.ylabel('CASM[-]')
    plt.title('CASM vs. Battery specific energy density')
    plt.legend()
    plt.grid(True)     
            

        
       
 
               
    
    ##================================================================================================================================================      
    ## Monthly Emissions
    ##================================================================================================================================================                  
    #colors              = px.colors.qualitative.Pastel 
    #sector_colors       = [colors[0],colors[2]] 
    #month_names         = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']      
    #fig_8               = go.Figure() 
    #fig_8.add_trace(go.Scatter(x=month_names, y=battery, name = 'Aircraft Fleet with Electric Aircraft',
                             #line=dict(color=sector_colors[0], width=4)))  
    #fig_8.add_trace(go.Scatter(x=month_names, y=no_battery, name='Aircraft Fleet without Electric Aircraft',
                             #line=dict(color=sector_colors[1], width=4)))   
    #fig_8.update_layout( 
                      #height           = 400, 
                      #width            = 600, 
                      #margin           = {'t':50,'l':0,'b':0,'r':0},
                      #yaxis_title_text ='CO2e (kg)', # yaxis label
                      #yaxis_range      = [0,10000000],
                      #font=dict(  size=font_size ),
                      #legend=dict(
                          #yanchor="top",
                          #y=0.99,
                          #xanchor="center",
                          #x=0.4 )) 
    
    ##================================================================================================================================================      
    ## Cost Per Seat Mile
    ##================================================================================================================================================   
    #fig_9 = go.Figure()       
    #fig_9.add_trace(go.Scatter(x=month_names, y=CASM_electric, name = 'Electric',
                             #line=dict(color=sector_colors[0], width=4)))  
    #fig_9.add_trace(go.Scatter(x=month_names, y=CASM_jet_A, name='Jet-A',
                             #line=dict(color=sector_colors[1], width=4)))  
    #fig_9.update_layout( 
                      #height           = 400, 
                      #width            = 600, 
                      #margin           = {'t':50,'l':0,'b':0,'r':0},
                      #yaxis_title_text ='Cost Per Seat Mile* (cents)', # yaxis label
                      #yaxis_range      = [0,10],
                      #font=dict(  size=font_size ),
                      #legend=dict(
                          #yanchor="top",
                          #y=0.99,
                          #xanchor="center",
                          #x=0.4 )) 
    
    #fig_4["layout"]["template"] = template   
    #fig_5["layout"]["template"] = template         
    #fig_6["layout"]["template"] = template    
    #fig_7["layout"]["template"] = template 
    #fig_8["layout"]["template"] = template  
    #fig_9["layout"]["template"] = template  
    
    return fig_4, fig_5, fig_6, fig_7, fig_8
     
if __name__ == '__main__':
    main()
    