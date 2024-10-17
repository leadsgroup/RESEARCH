# Missions.py

#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import RCAIDE

import sys
sys.path.append('Common')  

# ---------------------------------------------------------------------
#   Define the Configurations
# --------------------------------------------------------------------- 
def configs_setup(vehicle, tms_operation, f_idx):

    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------

    configs         = RCAIDE.Library.Components.Configs.Config.Container() 
    base_config     = RCAIDE.Library.Components.Configs.Config(vehicle)
    base_config.tag = 'base'  
    configs.append(base_config)

    return configs