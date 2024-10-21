# NMC_Battery_Cell

#----------------------------------------------------------------------
#   Imports
# --------------------------------------------------------------------- 
import sys 
import RCAIDE  
from RCAIDE.Core import Units, Data 
from RCAIDE.Methods.Power.Battery.Sizing import initialize_from_mass  
from RCAIDE.Components.Energy.Storages.Batteries import Battery
from RCAIDE.Core import Units 
from RCAIDE.Methods.Power.Battery.Sizing import initialize_from_circuit_configuration
from RCAIDE.Core import Data 
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------
def full_setup(c_rate,uncertainties,cycles,days): 
     
    # vehicle data
    vehicle  = vehicle_setup(c_rate,uncertainties)
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(configs_analyses,vehicle,c_rate,cycles,days)
    missions_analyses = missions_setup(mission)

    analyses = RCAIDE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    analyses.missions = missions_analyses


    return vehicle, analyses


# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------
def vehicle_setup(c_rate,uncertainties):

    vehicle                       = RCAIDE.Vehicle() 
    vehicle.tag                   = 'battery'   
    vehicle.reference_area        = 1

    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------    
    # mass properties
    vehicle.mass_properties.takeoff         = 0.048 * Units.kg 
    vehicle.mass_properties.max_takeoff     = 0.048 * Units.kg 
    
    # basic parameters
    vehicle.reference_area      = 1.   
    
    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------   
    wing                         = RCAIDE.Components.Wings.Wing()
    wing.tag                     = 'main_wing' 
    wing.areas.reference         = 1.
    wing.spans.projected         = 1.
    wing.aspect_ratio            = 1.
    wing.symmetric               = True
    wing.thickness_to_chord      = 0.12
    wing.taper                   = 1.
    wing.dynamic_pressure_ratio  = 1.
    wing.chords.mean_aerodynamic = 1.
    wing.chords.root             = 1.
    wing.chords.tip              = 1.
    wing.origin                  = [[0.0,0.0,0.0]] # meters
    wing.aerodynamic_center      = [0.0,0.0,0.0] # meters
    
    # add to vehicle
    vehicle.append_component(wing)
     

    net                           = RCAIDE.Components.Energy.Networks.Battery_Cell_Cycler()
    net.tag                       ='battery_cell'   
    net.dischage_model_fidelity   ='NMC'

    # Battery     
    bat                                      = RCAIDE.Components.Energy.Storages.Batteries.Constant_Mass.Lithium_Ion_LiNiMnCoO2_18650() 
         
    # define variable uncertainties  
    bat.uncertainties.internal_resistance    = uncertainties.R0 
    bat.cell.nominal_voltage                 = uncertainties.V0  
    bat.cell.nominal_capacity                = uncertainties.E      
    bat.specific_heat_capacity               = uncertainties.Cp                           
    bat.cell.specific_heat_capacity          = uncertainties.Cp    
    bat.charging_voltage                     = bat.cell.nominal_voltage   
    bat.convective_heat_transfer_coefficient = 7.17
    net.voltage                              = bat.cell.nominal_voltage 
    initialize_from_circuit_configuration(bat)  
    operating_current                        = bat.cell.nominal_capacity*c_rate
    bat.charging_current                     = operating_current
    net.battery                              = bat   
    
    # append battery mass to vehicle 
    vehicle.mass_properties.takeoff = bat.mass_properties.mass 

    avionics                      = RCAIDE.Components.Energy.Peripherals.Avionics()
    avionics.current              = operating_current
    net.avionics                  = avionics  

    vehicle.append_component(net)

    return vehicle

def analyses_setup(configs):

    analyses = RCAIDE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):   
    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Basic Geometry Relations
    sizing = RCAIDE.Analyses.Sizing.Sizing()
    sizing.features.vehicle = vehicle
    analyses.append(sizing)

    # ------------------------------------------------------------------
    #  Weights
    weights = RCAIDE.Analyses.Weights.Weights_eVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = RCAIDE.Analyses.Aerodynamics.Fidelity_Zero() 
    aerodynamics.vehicle = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.0000
    analyses.append(aerodynamics)  

    # ------------------------------------------------------------------	
    #  Stability Analysis	
    stability = RCAIDE.Analyses.Stability.Fidelity_Zero()    	
    stability.geometry = vehicle	
    analyses.append(stability) 

    # ------------------------------------------------------------------
    #  Energy
    energy= RCAIDE.Analyses.Energy.Energy()
    energy.network = vehicle.networks 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses    



def configs_setup(vehicle): 
    configs         = RCAIDE.Components.Configs.Config.Container()  
    base_config     = RCAIDE.Components.Configs.Config(vehicle)
    base_config.tag = 'base' 
    configs.append(base_config)   
    return configs

def mission_setup(analyses,vehicle,c_rate,cycles_per_day,days):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission     = RCAIDE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission'  
       
    # unpack Segments module
    Segments     = RCAIDE.Analyses.Mission.Segments

    # base segment
    base_segment                                                              = Segments.Segment() 
    base_segment.process.initialize.initialize_battery                        = RCAIDE.Methods.Missions.Segments.Common.Energy.initialize_battery 
    base_segment.process.finalize.post_process.update_battery_state_of_health = RCAIDE.Methods.Missions.Segments.Common.Energy.update_battery_state_of_health    
    
    if days == 0 and cycles_per_day == 0: 
        bat                                                      = vehicle.networks.battery_cell.battery    
        base_segment.max_energy                                  = bat.max_energy
        base_segment.charging_SOC_cutoff                         = bat.cell.charging_SOC_cutoff 
        base_segment.charging_current                            = bat.charging_current
        base_segment.charging_voltage                            = bat.charging_voltage  
        base_segment.state.numerics.number_control_points        = 30
        discharge_time                                           = 0.95 * (1/c_rate) * Units.hrs
          
        # Discharge Segment 
        segment                                                  = Segments.Ground.Battery_Charge_Discharge(base_segment) 
        segment.analyses.extend(analyses.base)            
        segment.tag                                              = 'NMC_Discharge'  
        segment.time                                             = discharge_time 
        segment.battery_energy                                   = bat.max_energy * 1.
        segment = vehicle.networks.battery_cell.add_unknowns_and_residuals_to_segment(segment,initial_battery_cell_temperature = 295 )    
        mission.append_segment(segment)          
    else:
        bat                                                      = vehicle.networks.battery_cell.battery    
        base_segment.max_energy                                  = bat.max_energy
        base_segment.charging_SOC_cutoff                         = bat.cell.charging_SOC_cutoff 
        base_segment.charging_current                            = bat.charging_current
        base_segment.charging_voltage                            = bat.charging_voltage   
        discharge_time                                           = 0.95 * (1/c_rate) * Units.hrs
    
        temp_dev                            = 20
        base_segment.temperature_deviation  = temp_dev
    
        for day in range(days):   
            for cycle in range(cycles_per_day):  
                segment     = Segments.Ground.Battery_Charge_Discharge(base_segment) 
                segment.tag = 'Bat_Discharge_Day_' + str (day) + '_Cycle_' + str (cycle)
                segment.analyses.extend(analyses.base)               
                segment.time                                        = discharge_time
                if  (day == 0) and (cycle== 0):
                    segment.battery_energy                               = bat.max_energy 
                    segment.initial_battery_resistance_growth_factor     = 1
                    segment.initial_battery_capacity_fade_factor         = 1 
                segment = vehicle.networks.battery_cell.add_unknowns_and_residuals_to_segment(segment,initial_battery_cell_temperature = 273 + temp_dev)    
                mission.append_segment(segment) 

                # Charge Model 
                segment     =  Segments.Ground.Battery_Charge_Discharge(base_segment)  
                segment.tag = 'Bat_Charge_Day_' + str (day) + '_Cycle_' + str (cycle)
                segment.analyses.extend(analyses.base)        
                segment.battery_discharge                           = False  
                if cycle  == cycles_per_day:  
                    segment.increment_battery_cycle_day=True                 
                segment.charging_SOC_cutoff                         =1.0
                segment = vehicle.networks.battery_cell.add_unknowns_and_residuals_to_segment(segment,initial_battery_cell_temperature = 273 + temp_dev)    
                mission.append_segment(segment)           
         
    return mission 

def missions_setup(base_mission):

    # the mission container
    missions = RCAIDE.Analyses.Mission.Mission.Container()

    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------
    missions.base = base_mission

    # done!
    return missions   
