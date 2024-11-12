
from RCAIDE.Library.Plots import *       
from RCAIDE import  load 
from RCAIDE import  save

import time  
import numpy as np
import pylab as plt 
import sys 
import os 
import pickle
import  pandas as pd 

sys.path.append(os.path.join( os.path.split(os.path.split(sys.path[0])[0])[0], 'Aircraft' + os.sep + 'Stopped_Rotor'))
from Stopped_Rotor  import vehicle_setup, configs_setup,analyses_setup , mission_setup, missions_setup, plot_results
 

def main():
     
    # vehicle data
    new_geometry    = True
    redesign_rotors =  False 
    if new_geometry :
        vehicle  = vehicle_setup(redesign_rotors)
        save_aircraft_geometry(vehicle , 'Stopped_Rotor')
    else: 
        vehicle = load_aircraft_geometry('Stopped_Rotor') 

    # Set up configs
    configs  = configs_setup(vehicle)

    # vehicle analyses
    analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(analyses)
    missions = missions_setup(mission) 
     
    results = missions.base_mission.evaluate() 
     
    # plot the results 
    plot_results(results) 

    ## plot vehicle 
    #plot_3d_vehicle(vehicle, 
                    #min_x_axis_limit            = -5,
                    #max_x_axis_limit            = 15,
                    #min_y_axis_limit            = -10,
                    #max_y_axis_limit            = 10,
                    #min_z_axis_limit            = -10,
                    #max_z_axis_limit            = 10,
                    #show_figure                 = False 
                    #)  
    
    
    return


def save_aircraft_geometry(geometry,filename): 
    pickle_file  = filename + '.pkl'
    with open(pickle_file, 'wb') as file:
        pickle.dump(geometry, file) 
    return 


def load_aircraft_geometry(filename):  
    load_file = filename + '.pkl' 
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results


def load_rotor(filename):
    rotor =  load(filename)
    return rotor

def save_rotor(rotor, filename):
    save(rotor, filename)
    return

if __name__ == '__main__': 
    main()    
    plt.show()