
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE
from RCAIDE.Framework.Core import Units,Data  
from RCAIDE.Library.Plots  import *  
import numpy as np  
# ----------------------------------------------------------------------------------------------------------------------
#  REGRESSION
# ----------------------------------------------------------------------------------------------------------------------  
def compute_control_surface_derivatives(vehicle):           
        
    elevator_span_fractions = np.array([0.2, 0.2, 0.6])
    CM_delta_e    =  np.zeros(4)
    CLift_delta_e =  np.zeros(4)
    for control_surface in vehicle.main_wing.control_surfaces:
        if type(control_surface) == RCAIDE.Library.Components.Wings.Control_Surfaces.Elevator:
            for i in  range(len(elevator_span_fractions)):
                # modify control surface
                control_surface.span_fraction_start = 0.1 
                control_surface.span_fraction_end   = elevator_span_fractions[i]
                
                # evaluate aircraft
                results =  evalaute_aircraft(vehicle)
                
                # store properties
                CM_delta_e[i]     = results.CM_delta_e_i 
                CLift_delta_e[i]  = results.CLift_delta_e_i
                
                
    derivatives =  Data(
        CM_delta_e = CM_delta_e, 
        CLift_delta_e = CLift_delta_e, 
    )
                
    return derivatives
                 
                
def evalaute_aircraft(vehicle): 
        
    # Set up configs
    configs  = configs_setup(vehicle)

    # vehicle analyses
    analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(analyses)
    missions = missions_setup(mission) 
     
    results = missions.base_mission.evaluate()
    
          
    return results
 
def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

       # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses        = RCAIDE.Framework.Analyses.Vehicle() 

    # ------------------------------------------------------------------
    #  Weights
    # ------------------------------------------------------------------
    weights         = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    # ------------------------------------------------------------------
    aerodynamics                                      = RCAIDE.Framework.Analyses.Aerodynamics.Vortex_Lattice_Method() 
    aerodynamics.vehicle                              = vehicle  
    aerodynamics.settings.use_surrogate               = False 
    aerodynamics.settings.trim_aircraft               = False 
    analyses.append(aerodynamics)
    
    # ------------------------------------------------------------------
    #  Energy
    # ------------------------------------------------------------------
    energy     = RCAIDE.Framework.Analyses.Energy.Energy()
    energy.vehicle = vehicle  
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    # ------------------------------------------------------------------
    planet     = RCAIDE.Framework.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    # ------------------------------------------------------------------
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses

  
def configs_setup(vehicle): 
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------ 
    configs = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config                                                       = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag                                                   = 'base'     
    configs.append(base_config) 
 
 

    return configs
 

# ----------------------------------------------------------------------
#   Define the Mission
# ---------------------------------------------------------------------- 
def mission_setup(analyses): 
    
    '''
    This sets up the nominal cruise of the aircraft
    '''
     
    mission = RCAIDE.Framework.Mission.Sequential_Segments()
    mission.tag = 'mission'
  
    # unpack Segments module
    Segments = RCAIDE.Framework.Mission.Segments

    #   Cruise Segment: constant Speed, constant altitude 
    segment                           = Segments.Untrimmed.Untrimmed()
    segment.analyses.extend( analyses )   
    segment.tag                       = "cruise"
    segment.angle_of_attack           = 1E-20
    segment.air_speed                 = 120 * Units.mph

    segment.flight_dynamics.force_x   = True    
    segment.flight_dynamics.force_z   = True    
    segment.flight_dynamics.force_y   = True     
    segment.flight_dynamics.moment_y  = True 
    segment.flight_dynamics.moment_x  = True
    segment.flight_dynamics.moment_z  = True
    
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