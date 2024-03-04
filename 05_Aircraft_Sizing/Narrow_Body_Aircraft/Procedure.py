# Procedure.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import numpy as np

import RCAIDE
from RCAIDE.Core import Units, Data
from RCAIDE.Analyses.Process import Process 
from RCAIDE.Methods.Energy.Propulsors.Turbofan_Propulsor   import design_turbofan,compute_nacelle_geometry
from RCAIDE.Optimization.Common             import write_optimization_outputs

# ----------------------------------------------------------------------        
#   Setup
# ----------------------------------------------------------------------   

def setup():
    
    # ------------------------------------------------------------------
    #   Analysis Procedure
    # ------------------------------------------------------------------ 
    
    # size the base config
    procedure = Process()
    procedure.resize_aircraft = resize_aircraft
    
    # find the weights
    procedure.weights = weight
    
    # performance studies
    procedure.missions                   = Process()
    procedure.missions.design_mission    = design_mission

    # post process the results
    procedure.post_process = post_process
        
    return procedure

# ----------------------------------------------------------------------        
#   Target Range Function
# ----------------------------------------------------------------------     

# ----------------------------------------------------------------------        
#   Design Mission
# ----------------------------------------------------------------------    
def design_mission(nexus):
    
    mission              = nexus.missions.base

    mission.design_range = 1500.*Units.nautical_miles        
    # Given a total target range, compute the cruise distance
    
    segments     = mission.segments
    climb_1      = segments['climb_1']
    climb_2      = segments['climb_2']
    climb_3      = segments['climb_3']
    climb_4      = segments['climb_4']
    climb_5      = segments['climb_5']
     
    descent_1    = segments['descent_1']
    descent_2    = segments['descent_2']
    descent_3    = segments['descent_3']

    x_climb_1    = climb_1.altitude_end/np.tan(np.arcsin(climb_1.climb_rate/climb_1.air_speed))
    x_climb_2    = (climb_2.altitude_end-climb_1.altitude_end)/np.tan(np.arcsin(climb_2.climb_rate/climb_2.air_speed))
    x_climb_3    = (climb_3.altitude_end-climb_2.altitude_end)/np.tan(np.arcsin(climb_3.climb_rate/climb_3.air_speed))
    x_climb_4    = (climb_4.altitude_end-climb_3.altitude_end)/np.tan(np.arcsin(climb_4.climb_rate/climb_4.air_speed))
    x_climb_5    = (climb_5.altitude_end-climb_4.altitude_end)/np.tan(np.arcsin(climb_5.climb_rate/climb_5.air_speed))
    x_descent_1  = (climb_5.altitude_end-descent_1.altitude_end)/np.tan(np.arcsin(descent_1.descent_rate/descent_1.air_speed))
    x_descent_2  = (descent_1.altitude_end-descent_2.altitude_end)/np.tan(np.arcsin(descent_2.descent_rate/descent_2.air_speed))
    x_descent_3  = (descent_2.altitude_end-descent_3.altitude_end)/np.tan(np.arcsin(descent_3.descent_rate/descent_3.air_speed)) 
    segments['cruise'].distance = mission.design_range-(x_climb_1+x_climb_2+x_climb_3+x_climb_4+x_climb_5+x_descent_1+x_descent_2+x_descent_3)
     
    results      = nexus.results
    
    # run mission 
    results.base = mission.evaluate()
    
    return nexus

# ----------------------------------------------------------------------        
#   Design Mission
# ----------------------------------------------------------------------    
def design_mission(nexus):
    mission              = nexus.missions.base

    results      = nexus.results
    
    # run mission 
    results.base = mission.evaluate()
    
    return nexus


# ----------------------------------------------------------------------        
#   Sizing
# ----------------------------------------------------------------------    

def resize_aircraft(nexus):
    #by this point, the aircraft paramters are updated 
    
    configs=nexus.vehicle_configurations
    base=configs.base
    
    # find conditions
    air_speed   = nexus.missions.base.segments['cruise'].air_speed 
    altitude    = nexus.missions.base.segments['climb_5'].altitude_end
    atmosphere  = RCAIDE.Analyses.Atmospheric.US_Standard_1976()
    
    freestream  = atmosphere.compute_values(altitude)
    freestream0 = atmosphere.compute_values(6000.*Units.ft)  #cabin altitude
    
    diff_pressure         = np.max(freestream0.pressure-freestream.pressure,0)
    fuselage              = base.fuselages['fuselage']
    fuselage.differential_pressure = diff_pressure 
    
    #now size engine
    mach_number        = air_speed/freestream.speed_of_sound
    
    #now add to freestream data object
    freestream.velocity    = air_speed
    freestream.mach_number = mach_number
    freestream.gravity     = 9.81
    
    conditions             = RCAIDE.Analyses.Mission.Common.Results()   #assign conditions in form for propulsor sizing
    conditions.freestream  = freestream
    
    for config in configs:
        config.wings.horizontal_stabilizer.areas.reference = (26.0/92.0)*config.wings.main_wing.areas.reference
            
        for wing in config.wings:
            
            wing = RCAIDE.Methods.Geometry.Two_Dimensional.Planform.wing_planform(wing)
            wing.areas.exposed  = 0.8 * wing.areas.wetted
            wing.areas.affected = 0.6 * wing.areas.reference
            
        fuselage                         = config.fuselages['fuselage']
        fuselage.differential_pressure = diff_pressure 
        
        for network in config.networks:
            for fuel_line in network.fuel_lines:
                for turbofan in fuel_line.propulsors: 
                    turbofan.design_mach_number = mach_number
                    turbofan.design_altitude    = altitude
                    design_turbofan(turbofan)
                    compute_nacelle_geometry(turbofan,turbofan.nacelle)

    return nexus

# ----------------------------------------------------------------------        
#   Weights
# ----------------------------------------------------------------------    

def weight(nexus):
    vehicle=nexus.vehicle_configurations.base

    # weight analysis
    weights = nexus.analyses.base.weights.evaluate(method="New SUAVE")
    weights = nexus.analyses.cruise.weights.evaluate(method="New SUAVE")
    vehicle.mass_properties.breakdown = weights
    weights = nexus.analyses.landing.weights.evaluate(method="New SUAVE")
    weights = nexus.analyses.takeoff.weights.evaluate(method="New SUAVE")
    weights = nexus.analyses.short_field_takeoff.weights.evaluate(method="New SUAVE")

    return nexus

# ----------------------------------------------------------------------
#   Finalizing Function
# ----------------------------------------------------------------------    

def finalize(nexus):
    
    nexus.analyses.finalize()   
    
    return nexus         

# ----------------------------------------------------------------------
#   Post Process Results to give back to the optimizer
# ----------------------------------------------------------------------   

def post_process(nexus):
    
    # Unpack data
    vehicle                           = nexus.vehicle_configurations.base
    results                           = nexus.results
    summary                           = nexus.summary
    nexus.total_number_of_iterations +=1
    
    # ----------------------------------------------------------------------       
    # BASE MISSION POST PROCESSING 
    # ----------------------------------------------------------------------   
    # power and throttle in base mission
    MTOW         = vehicle.mass_properties.max_takeoff
    S_ref        = vehicle.reference_area
    max_throttle = 0 
    max_power    = 0
    for i in range(len(results.base.segments)):  
        for network in results.base.segments[i].analyses.energy.networks: 
            fuel_lines      = network.fuel_lines 
            for fuel_line in fuel_lines:
                for propulsor in fuel_line.propulsors:  
                    max_segment_throttle = np.max(results.base.segments[i].conditions.energy[fuel_line.tag][propulsor.tag].throttle[:,0])
                    max_segment_power    = np.max(results.base.segments[i].conditions.energy.power[:,0])     
                    if max_segment_throttle > max_throttle:
                        max_throttle = max_segment_throttle   
                    if max_segment_power > max_power:
                        max_power = max_segment_power  
            
    summary.max_throttle = max_throttle
    
    # Fuel margin and base fuel calculations
    zero_fuel_weight         = vehicle.mass_properties.breakdown.zero_fuel_weight[0]
    design_landing_weight    = results.base.segments[-1].conditions.weights.total_mass[-1,0]
    max_zero_fuel_margin     = (design_landing_weight- zero_fuel_weight)/zero_fuel_weight
    base_mission_fuelburn    = vehicle.mass_properties.takeoff - results.base.segments['descent_3'].conditions.weights.total_mass[-1,0]
    summary.max_zero_fuel_margin  = max_zero_fuel_margin
    summary.base_mission_fuelburn = base_mission_fuelburn
    
    # wing loading and power loading 
    wing_loading          = MTOW*9.81/S_ref
    power_loading         = MTOW*9.81/max_power
    summary.wing_loading  = wing_loading
    summary.power_loading = power_loading
    

    # ----------------------------------------------------------------------       
    # PRINT OPTIMIZATION RESULTS
    # ----------------------------------------------------------------------   
    print('\n\nMTOW (kg)              : ',MTOW)
    print('Wing Area  (m^2)       : ',S_ref)
    print('Max zero-fuel margin   : ',max_zero_fuel_margin)
    print('fuel burn (kg)         : ',base_mission_fuelburn)
    print('Wing Loading  (kg/m^2) : ',wing_loading)
    print('Power Loading (kg/W)   : ',power_loading)
    
    # when you run want to output results to a file 
    #filename = 'results.txt' 
    #write_optimization_outputs(nexus, filename)    
    return nexus    
