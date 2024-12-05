#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE 
from RCAIDE.Framework.Core import Units, Data   
from  RCAIDE.Framework.Analyses.Geodesics.Geodesics import Calculate_Distance  
from RCAIDE.Library.Plots import *       
from RCAIDE import  load 
from RCAIDE import  save
from matplotlib import  pyplot as  plt
import  pickle

# python imports 
import os 
import pickle
import sys 
import pandas as pd
import numpy as  np 

local_path_1 =  os.path.split(sys.path[0])[0] 
local_path_2 =  os.path.split(os.path.split(sys.path[0])[0])[0]  

sys.path.append( os.path.join(local_path_1, 'Flight_Path_Functions'))
sys.path.append( os.path.join(local_path_1, 'Post_Processing_Functions'))
sys.path.append(os.path.join(local_path_2, 'Aircraft' + os.path.sep + 'Hexacopter'))
sys.path.append(os.path.join(local_path_2, 'Aircraft' + os.path.sep + 'Tiltrotor'))
sys.path.append(os.path.join(local_path_2, 'Aircraft' + os.path.sep + 'Tilt_Stopped_Rotor'))
from Hexacopter                       import vehicle_setup as  HC_vehicle_setup
from Hexacopter                       import configs_setup as  HC_configs_setup 
from Tiltrotor                        import vehicle_setup as  TR_vehicle_setup
from Tiltrotor                        import configs_setup as  TR_configs_setup 
from Tilt_Stopped_Rotor_V_Tail        import vehicle_setup as  TSR_vehicle_setup
from Tilt_Stopped_Rotor_V_Tail        import configs_setup as  TSR_configs_setup 
from Aircraft_Noise_Emissions         import approach_departure_flight_simulation_results , post_process_approach_departure_noise_data
# ----------------------------------------------------------------------------------------------------------------------
#  Main 
# ----------------------------------------------------------------------------------------------------------------------  
def main():            
     
    # ----------------------------------------------------------------------------------------------------------------------
    #  SIMULATION SETTINGS 
    # ----------------------------------------------------------------------------------------------------------------------  
    aircraft_codes                   = ['HC', 'TSR','TR']  
    cruise_altitude                  = 1000*Units.feet 
    radius_Vert1                     = 3600* Units.feet # circular pattern radius around vertiport 1
    radius_Vert2                     = 3600* Units.feet # circular pattern radius around vertiport 2
    dep_heading                      = 0   * Units.degree # Heading [degrees] of the departure from vertiport 1
    app_heading                      = 180  * Units.degree # Heading [degrees] of the approach to vertiport 2  
    number_of_cpts                   = 2  
    dep_sector                       = 180  * Units.degree      
    app_sector                       = 180  * Units.degree      
    path_heading                     = 270 *  Units.degree
    level_cruise_distance            = 10 * Units.miles
 
    # ----------------------------------------------------------------------------------------------------------------------
    #  RUN MISSION TO GET NOISE 
    # ----------------------------------------------------------------------------------------------------------------------
    raw_noise_filenames =  []
    analyzed_filenames  =  [] 
    
    for i in  range(3):
        aircraft_code      =  aircraft_codes[i]
        results            = run_noise_mission(number_of_cpts,aircraft_code, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, level_cruise_distance,cruise_altitude) # Run noise simulation  
        raw_noise_filename =  aircraft_code  +  '_Hemisphere_Raw_Data'
        raw_noise_filenames.append(raw_noise_filename)
        save(results, raw_noise_filename ,pickle_format=True)  
    
        
        # ----------------------------------------------------------------------------------------------------------------------
        #  DATA ANALYSIS 
        # ----------------------------------------------------------------------------------------------------------------------
        raw_noise_filename =  raw_noise_filenames[i]
        results =  load(raw_noise_filename, pickle_format=True)
        analyzed_results  =  approach_departure_flight_simulation_results(results) 
        analyzed_filename =  aircraft_code  +  '_Hemisphere_Analyzed_Noise'
        analyzed_filenames.append(analyzed_filename)
        save(analyzed_results, analyzed_filename + '.res')
        
        F =  Data(filename_list=analyzed_filenames)
        analyzed_filename_name = '_Hemisphere_Analyzed_Noise' 
        save(F, analyzed_filename_name + '.res')                  
              
    return


def run_noise_mission(number_of_cpts,aircraft_code, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, level_cruise_distance,cruise_altitude):

    if aircraft_code == 'HC':
        vehicle  = HC_vehicle_setup(redesign_rotors = False)     
        configs  = HC_configs_setup(vehicle)
    elif aircraft_code == 'TSR':    
        vehicle  = TSR_vehicle_setup(redesign_rotors = False)     
        configs  = TSR_configs_setup(vehicle)
    elif aircraft_code == 'TR':
        vehicle  = TR_vehicle_setup(redesign_rotors = False)     
        configs  = TR_configs_setup(vehicle)

    # vehicle analyses
    analyses = noise_analyses_setup(configs)

    # mission analyses

    if aircraft_code == 'HC':    
        mission  = HC_noise_mission_setup(number_of_cpts, analyses, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, level_cruise_distance,cruise_altitude)
    elif aircraft_code == 'TSR':
        mission  = TSR_noise_mission_setup(number_of_cpts, analyses, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, level_cruise_distance,cruise_altitude)
    elif aircraft_code == 'TR':
        mission  = TR_noise_mission_setup(number_of_cpts, analyses, radius_Vert1, radius_Vert2, dep_heading, app_heading, dep_sector, app_sector, path_heading, level_cruise_distance,cruise_altitude)
 
    missions = missions_setup(mission) 
     
    results = missions.base_mission.evaluate() 
    return results

def noise_analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses 

# ------------------------------------------------------------------
# Base Analysis
# ------------------------------------------------------------------
def base_analysis(vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = RCAIDE.Framework.Analyses.Vehicle() 
    
    # ------------------------------------------------------------------
    #  Weights
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics         = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle = vehicle 
    analyses.append(aerodynamics)
     
    # ------------------------------------------------------------------
    #  Stability Analysis'HC_mission_LA_ONT_BUR_1000ft'
    stability         = RCAIDE.Framework.Analyses.Stability.Vortex_Lattice_Method() 
    stability.vehicle = vehicle 
    analyses.append(stability)    

    #  Noise Analysis   
    noise = RCAIDE.Framework.Analyses.Noise.Frequency_Domain_Buildup()   
    noise.vehicle = vehicle    
    noise.settings.noise_hemisphere_phi_angles   = np.linspace(0,np.pi / 2,12)
    noise.settings.noise_hemisphere_theta_angles = np.linspace(-1 * np.pi, 1*np.pi,24)
    analyses.append(noise)
 
    # ------------------------------------------------------------------
    #  Energy
    energy          = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = RCAIDE.Framework.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses
 

# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def HC_noise_mission_setup(number_of_cpts, analyses, radius_Vert1=3600*Units.ft, radius_Vert2=3600*Units.ft, dep_heading=200*Units.degrees, app_heading=90*Units.degrees, dep_sector=90*Units.degrees, app_sector=90*Units.degrees, path_heading = 100, level_cruise_distance=10*Units.miles,cruise_altitude=1000*Units.ft): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment()
    base_segment.state.numerics.number_of_control_points    = number_of_cpts
    
    
    # ------------------------------------------------------------------
    #   Mission Constants
    # ------------------------------------------------------------------
    
    pattern_speed    = 55  * Units['mph'] #CHANGE FOR EACH AIRCRAFT 
    cruise_speed     = 75. * Units['mph'] #CHANGE FOR EACH AIRCRAFT
    transition_speed = 35. * Units['mph'] #CHANGE FOR EACH AIRCRAFT 
             
    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                                  = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                              = "Cruise"  
    segment.analyses.extend(analyses.forward_flight)   
    segment.initial_battery_state_of_charge                = 1.0
    segment.altitude                         = cruise_altitude      
    segment.air_speed                        = cruise_speed    
    segment.distance                         = level_cruise_distance
    segment.true_course                      = path_heading

    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True        
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['rotor_propulsor_1','rotor_propulsor_2','rotor_propulsor_3',
                                                                         'rotor_propulsor_4','rotor_propulsor_5','rotor_propulsor_6']]
    segment.assigned_control_variables.body_angle.active             = True     
    mission.append_segment(segment)      
    
                
    
    return mission


# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def TSR_noise_mission_setup(number_of_cpts, analyses, radius_Vert1=3600*Units.ft, radius_Vert2=3600*Units.ft, dep_heading=200*Units.degrees, app_heading=90*Units.degrees, dep_sector=90*Units.degrees, app_sector=90*Units.degrees, path_heading = 100, level_cruise_distance=10*Units.miles,cruise_altitude=1000*Units.ft): 
  
     
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission        = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag    = 'baseline_mission' 
    
    # unpack Segments module
    Segments       = RCAIDE.Framework.Mission.Segments 
    base_segment   = Segments.Segment()
    base_segment.state.numerics.number_of_control_points    = number_of_cpts
    

    
    # ------------------------------------------------------------------
    #   Mission Constants
    # ------------------------------------------------------------------
    
    pattern_speed    = 90 * Units.kts     #CHANGE FOR EACH AIRCRAFT 
    cruise_speed     = 110. * Units['mph'] #CHANGE FOR EACH AIRCRAFT 
     
    #------------------------------------------------------------------------------------------------------------------------------------  
    # Cruise 
    #------------------------------------------------------------------------------------------------------------------------------------  
    segment                                                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                                                      = "Cruise"  
    segment.analyses.extend( analyses.forward_flight )                    
    segment.initial_battery_state_of_charge                = 1.0
    segment.altitude                                                 = 1000.0 * Units.ft  
    segment.air_speed                                                = cruise_speed
    segment.distance                                                 = level_cruise_distance
    segment.true_course                                              = path_heading
                                                                     
    # define flight dynamics to model                                
    segment.flight_dynamics.force_x                                  = True  
    segment.flight_dynamics.force_z                                  = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['prop_rotor_propulsor_1','prop_rotor_propulsor_2','prop_rotor_propulsor_3',
                                                                         'prop_rotor_propulsor_4','prop_rotor_propulsor_5','prop_rotor_propulsor_6'] ] 
    segment.assigned_control_variables.body_angle.active             = True                
         
    mission.append_segment(segment)   
   
    return mission

# ------------------------------------------------------------------
#   Baseline Mission Setup
# ------------------------------------------------------------------
def TR_noise_mission_setup(number_of_cpts, analyses, radius_Vert1=3600*Units.ft, radius_Vert2=3600*Units.ft, dep_heading=200*Units.degrees, app_heading=90*Units.degrees, dep_sector=90*Units.degrees, app_sector=90*Units.degrees, path_heading = 100, level_cruise_distance=10*Units.miles,cruise_altitude=1000*Units.ft): 
    
    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'

    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments  
    base_segment = Segments.Segment()
    base_segment.state.numerics.number_of_control_points    = number_of_cpts
    

    # ------------------------------------------------------------------
    #   Mission Constants
    # ------------------------------------------------------------------
    
    pattern_speed    = 90 * Units.kts 
    cruise_speed     = 125.  * Units['mph']   
    transition_speed = 100 * Units['mph']    
    
 

    # ------------------------------------------------------------------
    #   First Cruise Segment: Constant Acceleration, Constant Altitude
    # ------------------------------------------------------------------ 
    segment                          = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag                      = "Cruise"  
    segment.analyses.extend(analyses.cruise) 
    segment.initial_battery_state_of_charge                = 1.0
    segment.altitude                 = cruise_altitude
    segment.air_speed                = cruise_speed
    segment.distance                 = level_cruise_distance
    segment.true_course              = path_heading
    
    # define flight dynamics to model 
    segment.flight_dynamics.force_x                       = True  
    segment.flight_dynamics.force_z                       = True     
    
    # define flight controls 
    segment.assigned_control_variables.throttle.active               = True           
    segment.assigned_control_variables.throttle.assigned_propulsors  = [['prop_rotor_propulsor_1','prop_rotor_propulsor_2','prop_rotor_propulsor_3',
                                                                         'prop_rotor_propulsor_4','prop_rotor_propulsor_5','prop_rotor_propulsor_6']]  
    segment.assigned_control_variables.body_angle.active             = True
    mission.append_segment(segment)     
   
    
    return mission


def missions_setup(mission): 
 
    missions         = RCAIDE.Framework.Mission.Missions()
    
    # base mission 
    mission.tag  = 'base_mission'
    missions.append(mission)
 
    return missions 
 
if __name__ == '__main__': 
    main()     