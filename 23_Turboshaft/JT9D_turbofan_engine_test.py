# RCAIDE imports 
import RCAIDE
from RCAIDE.Core import Units    
from RCAIDE.Methods.Propulsion.Design           import design_turbofan ,design_propeller
from RCAIDE.Visualization                       import *     
from RCAIDE.Methods.Propulsion.Performance.turbofan_propulsor  import compute_propulsor_performance   , compute_unique_propulsor_groups 

# python imports 
import numpy as np  
import pickle
from copy import deepcopy
import matplotlib.pyplot as plt  
import os   
import matplotlib.cm as cm

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main(): 

    comparison_of_perfomrance()
    #first_order_analysis()

    #second_order_analysis()
    return 

def comparison_of_perfomrance():
    baseline_file    =  'baseline'
    baseline_res     = load_results(baseline_file)

    hybrid_file   = 'series_parallel_hybrid_aircraft'
    hybrid_res  = load_results(hybrid_file) 


    # get plotting style 
    ps      = plot_style()  
    ps.axis_font_size = 14

    parameters = {'axes.labelsize': ps.axis_font_size,
                  'xtick.labelsize': ps.axis_font_size,
                  'ytick.labelsize': ps.axis_font_size,
                  'axes.titlesize': ps.title_font_size}
    plt.rcParams.update(parameters)

    # get line colors for plots     

    fig   = plt.figure()
    fig.set_size_inches(5,5) 
    axes_1 = fig.add_subplot(3,1,1)
    axes_2 = axes_1.twinx() 
    axes_3 = fig.add_subplot(3,1,2)
    axes_4 = fig.add_subplot(3,1,3)
    for i in range(len(baseline_res.segments)): 
        time     = baseline_res.segments[i].conditions.frames.inertial.time[:,0] / Units.min
        airspeed = baseline_res.segments[i].conditions.freestream.velocity[:,0] /   Units['mph']  
        altitude = baseline_res.segments[i].conditions.freestream.altitude[:,0]/Units.feet   

        Weight_baseline   = baseline_res.segments[i].conditions.weights.total_mass[:, 0] * 9.81   
        mdot_baseline     = baseline_res.segments[i].conditions.weights.vehicle_mass_rate[:, 0]
        thrust_baseline   = baseline_res.segments[i].conditions.frames.body.thrust_force_vector[:, 0]
        sfc_baseline      = (mdot_baseline / Units.lb) / (thrust_baseline / Units.lbf) * Units.hr    

        Weight_hybrid   = hybrid_res.segments[i].conditions.weights.total_mass[:, 0] * 9.81   
        mdot_hybrid     = hybrid_res.segments[i].conditions.weights.vehicle_mass_rate[:, 0]
        thrust_hybrid   = hybrid_res.segments[i].conditions.frames.body.thrust_force_vector[:, 0]
        sfc_hybrid      = (mdot_hybrid / Units.lb) / (thrust_hybrid / Units.lbf) * Units.hr      

        axes_1.plot(time, altitude, color = 'firebrick', marker = ps.marker, linewidth = ps.line_width)
        axes_1.set_ylabel(r'Altitude (ft)') 
        axes_1.set_ylim([0, 40000]) 
        set_axes(axes_1)     
        axes_2.plot(time, airspeed, color = 'mediumblue', marker = ps.marker, linewidth = ps.line_width)   
        axes_2.tick_params(axis="y", which="both", right=False, left=False, labelright=True,
                           color='mediumblue', labelcolor='mediumblue')        
        axes_2.set_ylabel(r"Airspeed (mph)" ,color='mediumblue')  
        axes_2.set_ylim([100, 1000])      

        if i == 0:
            axes_3.plot(time, Weight_baseline/1000 , color = 'black', marker = ps.marker, linewidth = ps.line_width, label = 'Baseline') 
            axes_3.plot(time, Weight_hybrid/1000 , color = 'green', marker = ps.marker, linewidth = ps.line_width, label = 'Hybrid') 
        else:
            axes_3.plot(time, Weight_baseline/1000 , color = 'black', marker = ps.marker, linewidth = ps.line_width) 
            axes_3.plot(time, Weight_hybrid/1000 , color = 'green', marker = ps.marker, linewidth = ps.line_width) 
        axes_3.set_ylabel(r'Weight (kN)') 
        set_axes(axes_3) 

        axes_4.plot(time, sfc_baseline, color = 'black', marker = ps.marker, linewidth = ps.line_width)
        axes_4.plot(time, sfc_hybrid, color = 'green', marker = ps.marker, linewidth = ps.line_width)
        axes_4.set_xlabel('Time (mins)')
        axes_4.set_ylabel(r'SFC (lb/lbf-hr)')
        axes_4.set_ylim([0, 1]) 
        set_axes(axes_4) 
        axes_3.legend()   

    fig.tight_layout()
    return



def test_propeller():


    #------------------------------------------------------------------------------------------------------------------------------------  
    # Propeller    
    #------------------------------------------------------------------------------------------------------------------------------------           
    propeller                                        = RCAIDE.Energy.Converters.Propeller() 
    propeller.tag                                    = 'propeller_1'  
    propeller.tip_radius                             = 8*Units.feet
    propeller.number_of_blades                       = 8
    propeller.hub_radius                             = 0.8
    propeller.cruise.design_freestream_velocity      = 343*0.76
    propeller.cruise.design_angular_velocity         = 2700. * Units.rpm 
    propeller.cruise.design_Cl                       = 0.5 
    propeller.cruise.design_altitude                 = 35000. * Units.feet 
    propeller.cruise.design_thrust                   = 10000 * Units.lbf   
    propeller.clockwise_rotation                     = False
    propeller.variable_pitch                         = True  
    propeller.origin                                 = [[2.,2.5,0.95]] 
    #airfoil                                          = RCAIDE.Components.Airfoils.Airfoil()    
    #airfoil.coordinate_file                          = '../../Vehicles/Airfoils/NACA_4412.txt'
    #airfoil.polar_files                              = ['../../Vehicles/Airfoils/Polars/NACA_4412_polar_Re_50000.txt' ,
                                                        #'../../Vehicles/Airfoils/Polars/NACA_4412_polar_Re_100000.txt' ,
                                                        #'../../Vehicles/Airfoils/Polars/NACA_4412_polar_Re_200000.txt' ,
                                                        #'../../Vehicles/Airfoils/Polars/NACA_4412_polar_Re_500000.txt' ,
                                                        #'../../Vehicles/Airfoils/Polars/NACA_4412_polar_Re_1000000.txt' ] 
    #propeller.append_airfoil(airfoil)              
    propeller.airfoil_polar_stations                 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] 
    propeller                                        = design_propeller(propeller,number_of_stations=10)    
    plot_3d_rotor(propeller, number_of_airfoil_points = 11)

    return 
def second_order_analysis():  

    altitude           = np.linspace(0,36000,10)*Units.feet 
    mach               = np.ones_like(altitude)*0.001 # np.linspace(0.1,0.7,100)

    motor_work = None 
    thrust_baseline , overall_efficiency_baseline  = JT9D_7_turbofan_engine(altitude,mach,motor_work) 

    fig = plt.figure('Turbofan')
    axis = fig.add_subplot(1,2,1)
    axis.set_ylabel('altitude')
    axis.set_xlabel('Thrust (lbf)')
    axis.plot(thrust_baseline/Units.lbf,altitude/Units.feet, color = 'black', label = 'baseline')    
    axis.grid()   

    axis2 = fig.add_subplot(1,2,2)
    axis2.set_ylabel('altitude')
    axis2.set_xlabel(r'overall efficency')
    axis2.plot(overall_efficiency_baseline*100 ,altitude/Units.feet, color = 'black', label = 'baseline')  
    axis2.grid()  


    markers = ['s','x','P','v','o','^','D','*','1']
    motor_work  = np.linspace(0,0.99,5)
    colors  = cm.inferno(np.linspace(0.2,0.8,len(motor_work)))  
    for i in range(len(motor_work)):
        thrust_motor , overall_efficiency_motor  = JT9D_7_turbofan_engine(altitude,mach,motor_work[i])    
        motor_tag = str(round(motor_work[i]*100,0)) +  '% $W_m$' 
        axis.plot(thrust_motor/Units.lbf,altitude/Units.feet,linestyle = '--', color =  colors[i], label = motor_tag  )     
        axis2.plot(overall_efficiency_motor*100,altitude/Units.feet,linestyle = '--', color =  colors[i] )    


        regenerator = RCAIDE.Energy.Thermal_Management.Turbine.Regenerator() 
        regenerator.combustor.temperature_ratio = 0.975
        regenerator.nozzle.temperature_ratio    = 0.975 
        motor_regen_tag =str(round(motor_work[i]*100,0)) +  '% $W_m$' +  r', $\tau_x$ = ' + str(regenerator.combustor.temperature_ratio)
        thrust_motor_regen , overall_efficiency_motor_regen  = JT9D_7_turbofan_engine(altitude,mach,motor_work[i],regenerator)   
        axis.plot(thrust_motor_regen/Units.lbf,altitude/Units.feet,linestyle = '-',  color =  colors[i], marker = markers[i], label = motor_regen_tag  )     
        axis2.plot(overall_efficiency_motor_regen*100,altitude/Units.feet, linestyle = '-',color =  colors[i], marker = markers[i]  )  

        #shaft_power_offtake = RCAIDE.Energy.Converters.Shaft_Power_Off_Take()
        #shaft_power_offtake.power_draw = 1000000 
        #regenerator = RCAIDE.Energy.Thermal_Management.Turbine.Regenerator() 
        #regenerator.combustor.temperature_ratio = 1.0
        #regenerator.nozzle.temperature_ratio    = 1.0 
        #motor_regen_tag ='$W_m$ = ' + str(motor_work[i]) +  r', $\tau_x$ = ' + str(regenerator.combustor.temperature_ratio) + ' 1MW Generator'
        #thrust_motor_regen_gen , overall_efficiency_motor_regen_gen  = JT9D_7_turbofan_engine(altitude,mach,motor_work[i],regenerator,shaft_power_offtake)   
        #axis3.plot(thrust_motor_regen_gen/Units.lbf,altitude/Units.feet,linestyle = '-',  color =  colors[i], marker = markers[i], label = motor_regen_tag  )     
        #axis4.plot(overall_efficiency_motor_regen_gen*100,altitude/Units.feet, linestyle = '-',color =  colors[i], marker = markers[i], label = motor_regen_tag )  


    axis.legend()           
    #axis2.legend()               
    fig.tight_layout()      
    return 


def first_order_analysis():

    work_ratio, eta_th             = brayton_cycle()
    work_ratio_motor, eta_th_motor = brayton_cycle_with_motor()
    work_ratio_regen, eta_th_regen = brayton_cycle_with_regenerator()
    work_ratio_regen_motor, eta_th_regen_motor =brayton_cycle_with_regenerator_and_motor()

    W_motor_additional = np.linspace(0,1,10)  

    fig = plt.figure()
    axis = fig.add_subplot(1,1,1)
    axis.set_xlabel('additional motor work')
    axis.set_ylabel('thermal efficiency work')
    axis.plot(W_motor_additional,eta_th, label = 'brayton' ) 
    axis.plot(W_motor_additional,eta_th_motor, label = 'brayton w. motor' ) 
    axis.plot(W_motor_additional,eta_th_regen, label = 'brayton regen' ) 
    axis.plot(W_motor_additional,eta_th_regen_motor, label = 'brayton regen w. motor' ) 
    axis.legend()
    axis.grid()

    fig.tight_layout()
    return 

def brayton_cycle(): 

    eta_comp    = 0.98
    eta_turb    = 0.98
    C_p         = 1.0
    gamma       = 1.4
    P_2_div_P_1 = 5 # compresson ratio back to free stream 
    P_3_div_P_4 = 5 # expansion ratio back to free stream 

    T_1         = kelvin(25) # compressor inlet temperature 
    T_3         = kelvin(850) # Turbine Inlet Temperature   

    T_2         = T_1*(P_2_div_P_1)**((gamma-1)/gamma)
    T_4         = T_3*(1/P_3_div_P_4)**((gamma-1)/gamma)
    W_comp      = (C_p/eta_comp)*(T_2-T_1)
    W_turb      = eta_turb*C_p*(T_3-T_4)
    work_ratio  = W_comp/W_turb  

    T_2    = W_comp/C_p + T_1
    q_in   = C_p*(T_3-T_2)   
    eta_th      = (W_turb -W_comp) / q_in

    return work_ratio*np.ones(10), eta_th*np.ones(10) 

def brayton_cycle_with_motor():

    eta_comp    = 0.98
    eta_turb    = 0.98
    C_p         = 1.0
    gamma       = 1.4
    P_2_div_P_1 = 5 # compresson ratio back to free stream 
    P_3_div_P_4 = 5 # expansion ratio back to free stream 

    T_1         = kelvin(25) # compressor inlet temperature 
    T_3         = kelvin(850) # Turbine Inlet Temperature   

    T_2         = T_1*(P_2_div_P_1)**((gamma-1)/gamma)
    T_4         = T_3*(1/P_3_div_P_4)**((gamma-1)/gamma)
    W_comp      = (C_p/eta_comp)*(T_2-T_1)
    W_turb      = eta_turb*C_p*(T_3-T_4)

    W_motor     = np.linspace(0,W_comp*0.9,10)  
    work_ratio  = (W_comp + W_motor)/W_turb  

    T_2    = W_comp/C_p + T_1
    q_in   = C_p*(T_3-T_2)

    eta_th      = (W_turb - (W_comp - W_motor)) / q_in
    return work_ratio, eta_th 

def  brayton_cycle_with_regenerator_and_motor(): 

    eta_comp    = 0.98
    eta_turb    = 0.98
    eta_regen   = 0.8
    C_p         = 1.0
    gamma       = 1.4
    P_2_div_P_1 = 5 # compresson ratio back to free stream 
    P_4_div_P_5 = 5 # expansion ratio back to free stream 

    T_1         = kelvin(25) # compressor inlet temperature 
    T_4         = kelvin(850) # Turbine Inlet Temperature   

    T_2         = T_1*(P_2_div_P_1)**((gamma-1)/gamma)
    T_5         = T_4*(1/P_4_div_P_5)**((gamma-1)/gamma)
    W_comp      = (C_p/eta_comp)*(T_2-T_1)
    W_turb      = eta_turb*C_p*(T_4-T_5)  
    W_motor     = np.linspace(0,W_comp*0.9,10) 
    work_ratio  = (W_comp + W_motor)/W_turb  
    T_2         = W_comp/C_p + T_1

    T_3         = (eta_regen*(T_5-T_2)) + T_2
    q_in        = C_p*(T_4-T_3)   

    eta_th  = (W_turb - (W_comp - W_motor)) / q_in

    return work_ratio , eta_th 


def  brayton_cycle_with_regenerator(): 

    eta_comp    = 0.98
    eta_turb    = 0.98
    eta_regen   = 0.8
    C_p         = 1.0
    gamma       = 1.4
    P_2_div_P_1 = 5 # compresson ratio back to free stream 
    P_4_div_P_5 = 5 # expansion ratio back to free stream 

    T_1         = kelvin(25) # compressor inlet temperature 
    T_4         = kelvin(850) # Turbine Inlet Temperature   

    T_2         = T_1*(P_2_div_P_1)**((gamma-1)/gamma)
    T_5         = T_4*(1/P_4_div_P_5)**((gamma-1)/gamma)
    W_comp      = (C_p/eta_comp)*(T_2-T_1)
    W_turb      = eta_turb*C_p*(T_4-T_5)
    work_ratio  = W_comp/W_turb  

    T_2    = W_comp/C_p + T_1

    T_3     = (eta_regen*(T_5-T_2)) + T_2
    q_in    = C_p*(T_4-T_3)   
    eta_th  = (W_turb -W_comp) / q_in

    return work_ratio*np.ones(10), eta_th*np.ones(10) 

def kelvin(T):

    return T + 273

def JT9D_7_turbofan_engine(altitude,mach,motor_work = None, regenerator = None, shaft_power_offtake = None ): 

    turbofan                                    = RCAIDE.Energy.Converters.Turbofan() 
    turbofan.tag                                = 'turbofan'
    turbofan.origin                             = [[13.72, 4.86,-1.1]] 
    turbofan.engine_length                      = 2.71     
    turbofan.bypass_ratio                       = 4.8 # checked  
    turbofan.design_altitude                    = 35000.0*Units.ft # checked  
    turbofan.design_mach_number                 = 0.78    # checked  
    turbofan.design_thrust                      = 10000*Units.lbf # checked   

    # fan                
    fan                                         = RCAIDE.Energy.Converters.Fan()   
    fan.tag                                     = 'fan'
    fan.polytropic_efficiency                   = 0.93
    fan.pressure_ratio                          = 1.56 
    turbofan.fan                                = fan        

    # working fluid                   
    turbofan.working_fluid                      = RCAIDE.Attributes.Gases.Air() 
    ram                                         = RCAIDE.Energy.Converters.Ram()
    ram.tag                                     = 'ram' 
    turbofan.ram                                = ram 

    # inlet nozzle          
    inlet_nozzle                                = RCAIDE.Energy.Converters.Compression_Nozzle()
    inlet_nozzle.tag                            = 'inlet nozzle'
    inlet_nozzle.polytropic_efficiency          = 0.98
    inlet_nozzle.pressure_ratio                 = 0.98 
    turbofan.inlet_nozzle                       = inlet_nozzle 

    # low pressure compressor    
    low_pressure_compressor                       = RCAIDE.Energy.Converters.Compressor()    
    low_pressure_compressor.tag                   = 'lpc'
    low_pressure_compressor.polytropic_efficiency = 0.91
    low_pressure_compressor.pressure_ratio        = 2.3 
    turbofan.low_pressure_compressor              = low_pressure_compressor

    # high pressure compressor  
    high_pressure_compressor                       = RCAIDE.Energy.Converters.Compressor()    
    high_pressure_compressor.tag                   = 'hpc'
    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 10.2
    if motor_work is not None:
        high_pressure_compressor.inputs.external_drive.work_done = motor_work
    turbofan.high_pressure_compressor              = high_pressure_compressor 

    if regenerator is not None:
        turbofan.regenerator = regenerator 

    if shaft_power_offtake is not None:
        turbofan.generator_flag    = True 
        turbofan.shaft_power_offtake = shaft_power_offtake

    # combustor  
    combustor                                      = RCAIDE.Energy.Converters.Combustor()   
    combustor.tag                                  = 'Comb'
    combustor.efficiency                           = 0.99 
    combustor.alphac                               = 1.0     
    combustor.turbine_inlet_temperature            = 1273.889
    combustor.pressure_ratio                       = 0.95
    combustor.fuel_data                            = RCAIDE.Attributes.Propellants.Jet_A()  
    turbofan.combustor                             = combustor

    # low pressure turbine  
    low_pressure_turbine                           = RCAIDE.Energy.Converters.Turbine()   
    low_pressure_turbine.tag                       ='lpt'
    low_pressure_turbine.mechanical_efficiency     = 0.99
    low_pressure_turbine.polytropic_efficiency     = 0.93 
    turbofan.low_pressure_turbine                  = low_pressure_turbine

    # high pressure turbine     
    high_pressure_turbine                          = RCAIDE.Energy.Converters.Turbine()   
    high_pressure_turbine.tag                      ='hpt'
    high_pressure_turbine.mechanical_efficiency    = 0.99
    high_pressure_turbine.polytropic_efficiency    = 0.93 
    turbofan.high_pressure_turbine                 = high_pressure_turbine 

    # core nozzle
    core_nozzle                                    = RCAIDE.Energy.Converters.Expansion_Nozzle()   
    core_nozzle.tag                                = 'core nozzle'
    core_nozzle.polytropic_efficiency              = 0.95
    core_nozzle.pressure_ratio                     = 0.99  
    turbofan.core_nozzle                           = core_nozzle

    # fan nozzle             
    fan_nozzle                                     = RCAIDE.Energy.Converters.Expansion_Nozzle()   
    fan_nozzle.tag                                 = 'fan nozzle'
    fan_nozzle.polytropic_efficiency               = 0.95
    fan_nozzle.pressure_ratio                      = 0.99 
    turbofan.fan_nozzle                            = fan_nozzle 

    # design turbofan
    design_turbofan(turbofan)  

    # ------------------------------------------------------------------------------------------------------------------------
    # Run engine 
    # ------------------------------------------------------------------------------------------------------------------------ 

    thrust             = np.zeros_like(altitude )
    overall_efficiency = np.zeros_like(altitude )
    for i in range(len(altitude)):

        atmosphere_sls = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
        atmo_data      = atmosphere_sls.compute_values(altitude[i],0.0)
        planet         = RCAIDE.Attributes.Planets.Earth()

        p   = atmo_data.pressure          
        T   = atmo_data.temperature       
        rho = atmo_data.density          
        a   = atmo_data.speed_of_sound    
        mu  = atmo_data.dynamic_viscosity      

        # setup conditions
        conditions_sls = RCAIDE.Analyses.Mission.Common.Results()            

        # freestream conditions    
        conditions_sls.freestream.altitude                    = np.atleast_2d(altitude[i])
        conditions_sls.freestream.mach_number                 = np.atleast_2d(mach[i])
        conditions_sls.freestream.pressure                    = np.atleast_2d(p)
        conditions_sls.freestream.temperature                 = np.atleast_2d(T)
        conditions_sls.freestream.density                     = np.atleast_2d(rho)
        conditions_sls.freestream.dynamic_viscosity           = np.atleast_2d(mu)
        conditions_sls.freestream.gravity                     = np.atleast_2d(planet.sea_level_gravity)
        conditions_sls.freestream.isentropic_expansion_factor = np.atleast_2d(turbofan.working_fluid.compute_gamma(T,p))
        conditions_sls.freestream.Cp                          = np.atleast_2d(turbofan.working_fluid.compute_cp(T,p))
        conditions_sls.freestream.R                           = np.atleast_2d(turbofan.working_fluid.gas_specific_constant)
        conditions_sls.freestream.speed_of_sound              = np.atleast_2d(a)
        conditions_sls.freestream.velocity                    = np.atleast_2d(a*mach[i])  

        net          = RCAIDE.Energy.Networks.Turbofan_Engine()  
        fuel_line    = RCAIDE.Energy.Distributors.Fuel_Line()  
        fuel_line.turbofans.append(turbofan)       
        net.fuel_lines.append(fuel_line)   

        conditions_sls.energy[fuel_line.tag]                                = RCAIDE.Analyses.Mission.Common.Conditions()         
        conditions_sls.noise[fuel_line.tag]                                 = RCAIDE.Analyses.Mission.Common.Conditions()      
        sorted_propulsors                                                   = compute_unique_propulsor_groups(fuel_line)
        fuel_line_results                                                   = conditions_sls.energy[fuel_line.tag] 
        fuel_line_results[turbofan.propulsor_group]                         = RCAIDE.Analyses.Mission.Common.Conditions()
        fuel_line_results[turbofan.propulsor_group].turbofan                = RCAIDE.Analyses.Mission.Common.Conditions() 
        fuel_line_results[turbofan.propulsor_group].unique_turbofan_tags    = sorted_propulsors.unique_turbofan_tags 
        fuel_line_results.N_turbofans                                       = sorted_propulsors.N_turbofans
        noise_results                                                       = conditions_sls.noise[fuel_line.tag]
        noise_results[turbofan.propulsor_group]                             = RCAIDE.Analyses.Mission.Common.Conditions() 
        noise_results[turbofan.propulsor_group].turbofan                    = RCAIDE.Analyses.Mission.Common.Conditions() 
        fuel_line_results[turbofan.propulsor_group].turbofan.throttle       = np.array([[1.0]])  

        T , P, mdot,sfc,eta_prop,eta = compute_propulsor_performance(0,fuel_line,turbofan.propulsor_group,fuel_line.turbofans,1,conditions_sls)   

        thrust[i]             = np.linalg.norm(T)
        overall_efficiency[i] = eta 

    return thrust , overall_efficiency
# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def save_results(results,filename): 
    pickle_file  =  filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return   

# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results  



if __name__ == '__main__': 
    main()
    plt.show()