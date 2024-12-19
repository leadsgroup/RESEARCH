# Analyses.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import RCAIDE
from RCAIDE.Framework.Core import Units

import numpy as np 

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ---------------------------------------------------------------------- 
def analyses_setup(configs):

    analyses = RCAIDE.Framework.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config, configs)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle, configs):

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
    planet     = RCAIDE.Framework.Analyses.Planets.Earth()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    # ------------------------------------------------------------------
    atmosphere = RCAIDE.Framework.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    # done!
    return analyses

 