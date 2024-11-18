# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
# RCAIDE imports 
import RCAIDE

# python imports 
import sys
#sys.path.append('../../../Aircraft/Cessna')
#import Cessna_172_mission_simulation
#sys.path.append('../../../Aircraft/Navion')
#import Navion_regression
sys.path.append('../../../Aircraft/C-5a')
import Lockheed_C5a

def main():
    
    #vehicle = Cessna_172_mission_simulation.vehicle_setup()
    
    ## ------------------------------------------------------------------
    ##   Weight Breakdown 
    ## ------------------------------------------------------------------  
    #weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_General_Aviation()
    #weight_analysis.vehicle                       = vehicle
    #weight_analysis.method                        = 'Raymer'
    #results                                       = weight_analysis.evaluate() 
    #print("Operating empty weight estimate for the aircraft: "+str(results))
    
    ####################################################################
    
    #vehicle = Navion_regression.vehicle_setup()
    
    ## ------------------------------------------------------------------
    ##   Weight Breakdown 
    ## ------------------------------------------------------------------  
    #weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_General_Aviation()
    #weight_analysis.vehicle                       = vehicle
    #weight_analysis.method                        = 'Raymer'
    #results                                       = weight_analysis.evaluate() 
    #print("Operating empty weight estimate for the aircraft: "+str(results))    
    
    ####################################################################
    
    vehicle = Lockheed_C5a.vehicle_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    weight_analysis.method                        = 'Raymer'
    results                                       = weight_analysis.evaluate() 
    print("Operating empty weight estimate for the aircraft: "+str(results))     
    
if __name__ == '__main__': 
    main()   