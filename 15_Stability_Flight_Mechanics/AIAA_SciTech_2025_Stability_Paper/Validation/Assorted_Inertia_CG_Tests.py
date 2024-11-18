# Vehicle.py 

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------     

import RCAIDE
from RCAIDE.Framework.Core                              import Units   
from RCAIDE.Library.Methods.Propulsors.Converters.Rotor import design_propeller
from RCAIDE.Library.Methods.Geometry.Planform           import segment_properties
from RCAIDE.Library.Plots                               import *
import numpy as np
import os

# ----------------------------------------------------------------------------------------------------------------------
#  IMPORT
# ---------------------------------------------------------------------------------------------------------------------- 
# RCAIDE imports 
import RCAIDE
from RCAIDE.Framework.Core import Units
from RCAIDE.Library.Plots                 import *     
from RCAIDE.Library.Methods.Weights.Correlation_Buildups import Common
from RCAIDE.Library.Methods.Weights.Moment_of_Inertia.compute_aircraft_moment_of_inertia import compute_aircraft_moment_of_inertia
from RCAIDE.Library.Methods.Weights.Moment_of_Inertia.compute_aircraft_moment_of_inertia import compute_cuboid_moment_of_inertia
from RCAIDE.Library.Methods.Weights.Center_of_Gravity     import compute_vehicle_center_of_gravity


# python imports 
import numpy as np  
import matplotlib.pyplot as plt

from  Cessna_172 import  vehicle_setup as  C172_vehicle_setup
from Twin_Otter import  vehicle_setup as  twin_otter_setup
from Boeing_737_800 import vehicle_setup as  B738_setup
from Lockheed_C5a import vehicle_setup as C5a_setup
from Joby import  vehicle_setup as  Joby_setup
from Stopped_Rotor import vehicle_setup as stopped_rotor_setup
from Boeing_777_200 import  vehicle_setup as  B777_setup
from  ATR_72 import  vehicle_setup as  ATR72_setup
from  CRJ_700_Mission_Simulation import  vehicle_setup as  CRJ7_setup

def main():

    #print("\n ############################\n")
    #C172_Analysis()
    #print("\n ############################\n")
    #Twin_Otter_Analysis()
    #print("\n ############################\n")
    #B738_Analysis()
    #print("\n ############################\n")    
    #C5a_Analysis()
    print("\n ############################\n")
    Joby_Analysis()
    #print("\n ############################\n")
    #Stop_Rotor_Analysis()
    #print("\n ############################\n")
    #B777_Analysis()
    #print("\n ############################\n")
    #ATR72_Analysis()
    #print("\n ############################\n")
    #CRJ7_Analysis()
    print("\n ############################\n")    
    
    return
    
def C172_Analysis():
    vehicle = C172_vehicle_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_General_Aviation()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Operating empty weight estimate for C172: "+str(results))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location      = vehicle.mass_properties.center_of_gravity
    CG_location      =  [[2.2, 0, 0]]
    print("CG location: " + str(CG_location))
    
    # ------------------------------------------------------------------
    #   Operating Aircraft MOI
    # ------------------------------------------------------------------    
    MOI, mass = compute_aircraft_moment_of_inertia(vehicle, CG_location)
  
    print(MOI)
    sft2     = 1.355817 # 1 slug*ft^2 to 1.355817 kg*m^2
    Cessna_true = np.array([[948.0 , 0, 0], [0, 1346.0, 0], [0, 0, 1967.0]]) * sft2
    error    = (MOI - Cessna_true) / Cessna_true * 100
    print(error)
    print("Percent of empty mass used for Cessna inertia tensor calcs: "+str((mass)/results['empty']*100)+"%")
    
def Twin_Otter_Analysis():
    vehicle = twin_otter_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Operating empty weight estimate for DHC-6 Twin Otter: "+str(results))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location      = vehicle.mass_properties.center_of_gravity
    CG_location      =  [[5.95, 0, 0.95]]
    print("CG location: " + str(CG_location))
    
    # ------------------------------------------------------------------
    #   Operating Aircraft MOI
    # ------------------------------------------------------------------    
    MOI, mass = compute_aircraft_moment_of_inertia(vehicle, CG_location)
  
    print(MOI)
    sft2     = 1.355817 # 1 slug*ft^2 to 1.355817 kg*m^2
    MOI_true = np.array([[15680.0 , 0, 620.], [0, 22106.0, 0], [620., 0, 33149.0]]) * sft2
    error    = (MOI - MOI_true) / MOI_true * 100
    print(error)
    print("Percent of empty mass used for Twin Otter inertia tensor calcs: "+str((mass)/results['empty']*100)+"%")

def B738_Analysis():
    vehicle = B738_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Operating empty weight estimate for Boeing 737: " + str(results["operating_empty"]))
    print("Empty weight estimate for Boeing 737: " + str(results["empty"]))
    print("Max takeoff weight estimate for Boeing 737: " + str(results["max_takeoff"]))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location      = vehicle.mass_properties.center_of_gravity
    print("CG location: " + str(CG_location))
    
    # ------------------------------------------------------------------
    #   Operating Aircraft MOI
    # ------------------------------------------------------------------    
    MOI, mass = compute_aircraft_moment_of_inertia(vehicle, CG_location)
  
    print(MOI)
    #sft2     = 1.355817 # 1 slug*ft^2 to 1.355817 kg*m^2
    #Cessna_true = np.array([[948.0 , 0, 0], [0, 1346.0, 0], [0, 0, 1967.0]]) * sft2
    #error    = (MOI - Cessna_true) / Cessna_true * 100
    #print(error)
    print("Percent of empty mass used for Boeing 737 inertia tensor calcs: "+str((mass)/results['empty']*100)+"%")
    
def C5a_Analysis():
    vehicle = C5a_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Operating empty weight estimate for C-5a: " + str(results["operating_empty"]))
    print("Empty weight estimate for C-5a: " + str(results["empty"]))
    print("Max takeoff weight estimate for C-5a: " + str(results["max_takeoff"]))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location      = vehicle.mass_properties.center_of_gravity
    print("CG location: " + str(CG_location))
    
    # ------------------------------------------------------------------
    #   Operating Aircraft MOI
    # ------------------------------------------------------------------    
    MOI, mass = compute_aircraft_moment_of_inertia(vehicle, CG_location)
  
    print(MOI)
    sft2     = 1.355817 # 1 slug*ft^2 to 1.355817 kg*m^2
    MOI_true = np.array([[19.1e6 , 0, 2.5e6], [0, 31.3e6, 0], [2.5e6, 0, 47.0e6]]) * sft2
    error    = (MOI - MOI_true) / MOI_true * 100
    print(error)
    print("Percent of empty mass used for C-5a inertia tensor calcs: "+str((mass)/results['empty']*100)+"%")

    
def Joby_Analysis():
    vehicle = Joby_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Empty weight estimate for Joby: " + str(results["empty"]))
    print("Total weight estimate for Joby: " + str(results["total"]))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location      = [[1.8, 0, 0.8]]
    print("CG location: " + str(CG_location))
    
    # ------------------------------------------------------------------
    #   Operating Aircraft MOI
    # ------------------------------------------------------------------    
    MOI, mass = compute_aircraft_moment_of_inertia(vehicle, CG_location)
  
    print(MOI)
    #sft2     = 1.355817 # 1 slug*ft^2 to 1.355817 kg*m^2
    #Cessna_true = np.array([[948.0 , 0, 0], [0, 1346.0, 0], [0, 0, 1967.0]]) * sft2
    #error    = (MOI - Cessna_true) / Cessna_true * 100
    #print(error)
    print("Percent of empty mass used for Joby inertia tensor calcs: "+str((mass)/results['empty']*100)+"%")
    
def Stop_Rotor_Analysis():
    vehicle = stopped_rotor_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_EVTOL()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Empty weight estimate for stop rotor: " + str(results["empty"]))
    print("Total weight estimate for stop rotor: " + str(results["total"]))
    
    # ------------------------------------------------------------------
    #   CG Location
    # ------------------------------------------------------------------    
    compute_vehicle_center_of_gravity(vehicle) 
    CG_location      = [[1.8, 0, 0.8]]
    print("CG location: " + str(CG_location))
    
    # ------------------------------------------------------------------
    #   Operating Aircraft MOI
    # ------------------------------------------------------------------    
    MOI, mass = compute_aircraft_moment_of_inertia(vehicle, CG_location)
  
    print(MOI)
    #sft2     = 1.355817 # 1 slug*ft^2 to 1.355817 kg*m^2
    #Cessna_true = np.array([[948.0 , 0, 0], [0, 1346.0, 0], [0, 0, 1967.0]]) * sft2
    #error    = (MOI - Cessna_true) / Cessna_true * 100
    #print(error)
    print("Percent of empty mass used for stop rotor inertia tensor calcs: "+str((mass)/results['empty']*100)+"%")

def B777_Analysis():
    vehicle = B777_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Empty weight estimate for B777-200: " + str(results["empty"]))
    print("Operating empty weight estimate for B777-200: " + str(results["operating_empty"]))
    print("Max takeoff weight estimate for B777-200: " + str(results["max_takeoff"]))

def ATR72_Analysis():
    vehicle = ATR72_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Empty weight estimate for ATR72: " + str(results["empty"]))
    print("Operating empty weight estimate for ATR72: " + str(results["operating_empty"]))
    print("Max takeoff weight estimate for ATR72: " + str(results["max_takeoff"]))

def CRJ7_Analysis():
    vehicle = CRJ7_setup()
    
    # ------------------------------------------------------------------
    #   Weight Breakdown 
    # ------------------------------------------------------------------  
    weight_analysis                               = RCAIDE.Framework.Analyses.Weights.Weights_Transport()
    weight_analysis.vehicle                       = vehicle
    results                                       = weight_analysis.evaluate() 
    print("Empty weight estimate for CRJ-700: " + str(results["empty"]))
    print("Operating empty weight estimate for CRJ-700: " + str(results["operating_empty"]))
    print("Max takeoff weight estimate for CRJ-700: " + str(results["max_takeoff"]))


if __name__ == '__main__': 
    main()