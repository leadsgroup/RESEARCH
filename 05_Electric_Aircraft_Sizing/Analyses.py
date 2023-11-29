# Analyses.py
# 
# Created:  Mar 2023, M. Clarke

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import MARC  

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ---------------------------------------------------------------------- 
def setup(configs):

    analyses = MARC.Analyses.Analysis.Container()

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
     
    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = MARC.Analyses.Stability.Fidelity_Zero()
    stability.geometry = vehicle
    analyses.append(stability)       
                                                                              
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