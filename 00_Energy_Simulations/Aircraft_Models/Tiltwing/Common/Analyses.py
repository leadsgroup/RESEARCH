'''
# Analyses.py
# 
# Created: May 2019, M Clarke
#          Sep 2020, M.Clarke 

'''

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import MARC  

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ---------------------------------------------------------------------- 
def analyses_setup(configs,run_noise_analysis_flag,use_topology_flag,microphone_terrain_data,airport_geospacial_data):

    analyses = MARC.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config,run_noise_analysis_flag,use_topology_flag,microphone_terrain_data,airport_geospacial_data)
        analyses[tag] = analysis

    return analyses

# ------------------------------------------------------------------
# Base Analysis
# ------------------------------------------------------------------
def base_analysis(vehicle,run_noise_analysis_flag,use_topology_flag,microphone_terrain_data,airport_geospacial_data):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = MARC.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Basic Geometry Relations
    sizing = MARC.Analyses.Sizing.Sizing()
    sizing.features.vehicle = vehicle
    analyses.append(sizing)

    # ------------------------------------------------------------------
    #  Weights
    weights = MARC.Analyses.Weights.Weights_eVTOL()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = MARC.Analyses.Aerodynamics.Fidelity_Zero()
    aerodynamics.geometry = vehicle  
    aerodynamics.settings.model_fuselage = True 
    aerodynamics.settings.number_spanwise_vortices           = 25
    aerodynamics.settings.number_chordwise_vortices          = 5    
    analyses.append(aerodynamics)  
        
    if run_noise_analysis_flag:  
        # ------------------------------------------------------------------
        #  Noise Analysis
        noise = MARC.Analyses.Noise.Fidelity_One()   
        noise.geometry = vehicle

        # ------------------------------------------------------------------
        #  Noise Analysis
        noise = MARC.Analyses.Noise.Fidelity_One()   
        noise.geometry = vehicle  
        noise.settings.mean_sea_level_altitude           = False 
        noise.settings.ground_microphone_x_resolution    = microphone_terrain_data.ground_microphone_x_resolution           
        noise.settings.ground_microphone_y_resolution    = microphone_terrain_data.ground_microphone_y_resolution             
        noise.settings.ground_microphone_min_x           = microphone_terrain_data.ground_microphone_min_x                 
        noise.settings.ground_microphone_max_x           = microphone_terrain_data.ground_microphone_max_x                 
        noise.settings.ground_microphone_min_y           = microphone_terrain_data.ground_microphone_min_y                 
        noise.settings.ground_microphone_max_y           = microphone_terrain_data.ground_microphone_max_y    
        noise.settings.ground_microphone_x_stencil       = microphone_terrain_data.ground_microphone_x_stencil             
        noise.settings.ground_microphone_y_stencil       = microphone_terrain_data.ground_microphone_y_stencil        
        
        if use_topology_flag:                             
            noise.settings.ground_microphone_locations   = microphone_terrain_data.cartesian_microphone_locations 
            noise.settings.aircraft_departure_location   = airport_geospacial_data.departure_location   
            noise.settings.aircraft_destination_location = airport_geospacial_data.destination_location   
        
        analyses.append(noise)                                                       
                                                                              
    # ------------------------------------------------------------------
    #  Energy
    energy= MARC.Analyses.Energy.Energy()
    energy.network = vehicle.networks
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = MARC.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = MARC.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    return analyses   
