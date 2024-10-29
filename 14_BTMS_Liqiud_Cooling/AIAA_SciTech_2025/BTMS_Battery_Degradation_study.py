
#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import pickle
import time  
import pylab as plt
import os
import sys



sys.path.append(os.path.join(sys.path[0], 'Common'))
import Vehicle
import Analyses 
import Missions
import Plots
import Create_Excel
 

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():  
    # start simulation clock
    ti                         = time.time()
    RUN_NEW_MODEL_FLAG         = True
    
    # -------------------------------------------------------------------------------------------    
    # SET UP SIMULATION PARAMETERS   
    # -------------------------------------------------------------------------------------------  
    days_per_group             = 1#5               # total number of days simulated
    flights_per_day            = 1#5               # number of flights per day
    day_group                  = [2,3]#,3,4,5,6,7,8,9,10,11]
    plot_mission               = False            # plot mission flag  
    resize_aircraft            = False
   
    if RUN_NEW_MODEL_FLAG:    
        for g_idx, group  in  enumerate(day_group):
            print(' ***********  Simulating Group ' + str(group) + ' ***********  ')
            # -------------------------------------------------------------------------------------------    
            # SET UP VEHICLE
            # -------------------------------------------------------------------------------------------  
            vehicle = Vehicle.vehicle_setup(resize_aircraft,'e_twin_otter_vehicle')
            configs = Vehicle.configs_setup(vehicle)
            analyses = Analyses.analyses_setup(configs)
            charge_throughput = load_charge_throughput()
            
           
            # -------------------------------------------------------------------------------------------    
            # SET UP MISSION PROFILE  
            # -------------------------------------------------------------------------------------------    
            base_mission      = Missions.repeated_flight_operation_setup(configs,analyses,day_group,g_idx,group,days_per_group, flights_per_day, charge_throughput)
            missions_analyses = Missions.missions_setup(base_mission)
       
        
            ## -------------------------------------------------------------------------------------------    
            # RUN SIMULATION !!
            # -------------------------------------------------------------------------------------------
            results     = missions_analyses.base.evaluate()
            charge_throughput[str(group)] = results.segments[-1].conditions.energy.bus.battery_modules.lithium_ion_nmc.cell.charge_throughput[-1]
            save_charge_throughput(charge_throughput)
            
            # -------------------------------------------------------------------------------------------    
            # SAVE RESULTS
            # -------------------------------------------------------------------------------------------
            filename          = 'e_Twin_Otter_'
            save_results(results,filename,group)
            create_excel(filename,group)
            
            if plot_mission: 
                Plots.plot_results(results,save_figure_flag = False)       
                        
    tf = time.time() 
    print ('time taken: '+ str(round(((tf-ti)/60),3)) + ' mins')
    
    
    return 
# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def create_excel(filename,group):
    results = load_results(filename,group)
    Create_Excel.write_data(results,filename,group)
    return


# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def save_results(results,filename,group):
   #  Pickle Backup Files
    pickle_file = os.path.join("Raw_Data", filename + 'group_number' + str(group) + '.pkl')
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return

# Function to load the existing pickle file, if it exists
def load_charge_throughput(filename='charge_throughput.pkl'):
    if os.path.exists(filename):
        with open(filename, 'rb') as f:
            return pickle.load(f)
    else:
        return {}  # Return an empty dictionary if the file does not exist
    
# Function to save the updated dictionary to the pickle file
def save_charge_throughput(data, filename='charge_throughput.pkl'):
    with open(filename, 'wb') as f:
        pickle.dump(data, f)


# ------------------------------------------------------------------
#   Load Results
# ------------------------------------------------------------------   
def load_results(filename, group):  
    load_file =  os.path.join("Raw_Data", filename + 'group_number' + str(group) + '.pkl')
    with open(load_file, 'rb') as file:
        results = pickle.load(file) 
    return results

  
    

def show_notification():
    os.system('osascript -e \'display notification "The simulation has completed successfully." with title "Simulation Complete"\'')


if __name__ == '__main__':
    main()
    show_notification()
    plt.show()