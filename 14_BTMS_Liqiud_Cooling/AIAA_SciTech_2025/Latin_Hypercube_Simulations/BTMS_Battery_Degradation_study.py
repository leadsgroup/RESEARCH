# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
import argparse
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
#   Parse Command-Line Arguments
# ----------------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(description="Degradation Simulator for Electric Aircraft Batteries")
    parser.add_argument("HAS_power", type=float, help="Power for the Heat-Activated System (HAS)")
    parser.add_argument("HEX_power", type=float, help="Power for the Heat Exchanger (HEX)")
    parser.add_argument("RES_dimensions", type=float, help="Dimensions for the reservoir (RES)")
    parser.add_argument("--storage_dir", type=str, required=True, help="Directory to store results")

    return parser.parse_args()

# ----------------------------------------------------------------------
#   Degradation Simulator
# ----------------------------------------------------------------------
def degradation_simulator(HAS_power, HEX_power, RES_dimensions, storage_dir):  
    # Start simulation clock
    ti = time.time()
    RUN_NEW_MODEL_FLAG = True
    
    # -------------------------------------------------------------------------------------------    
    # SET UP SIMULATION PARAMETERS   
    # -------------------------------------------------------------------------------------------  
    days_per_group = 5  # total number of days simulated
    flights_per_day = 6  # number of flights per day
    day_group = list(range(1, 31)) 
    plot_mission = True  # plot mission flag  
   
    if RUN_NEW_MODEL_FLAG:    
        for g_idx, group in enumerate(day_group):
            print(f' ***********  Simulating Group {group} ***********  ')
            
            # -------------------------------------------------------------------------------------------    
            # SET UP VEHICLE
            # -------------------------------------------------------------------------------------------  
            resize_aircraft = group == 1
            vehicle = Vehicle.vehicle_setup(resize_aircraft, 'e_twin_otter_vehicle')
            configs = Vehicle.configs_setup(vehicle)
            analyses = Analyses.analyses_setup(configs)
            charge_throughput, cycle_day, resistance_growth, capacity_fade = load_charge_throughput(storage_dir)
            
            # -------------------------------------------------------------------------------------------    
            # SET UP MISSION PROFILE  
            # -------------------------------------------------------------------------------------------    
            base_mission = Missions.repeated_flight_operation_setup(
                configs, analyses, day_group, g_idx, group, days_per_group, flights_per_day,
                charge_throughput, cycle_day, resistance_growth, capacity_fade
            )
            missions_analyses = Missions.missions_setup(base_mission)
        
            # -------------------------------------------------------------------------------------------    
            # RUN SIMULATION !!
            # -------------------------------------------------------------------------------------------
            results = missions_analyses.base.evaluate()
            charge_throughput[str(group)] = results.segments[-1].conditions.energy.bus.battery_modules.lithium_ion_nmc.cell.charge_throughput[-1]
            resistance_growth[str(group)] = results.segments[-1].conditions.energy.bus.battery_modules.lithium_ion_nmc.cell.resistance_growth_factor
            capacity_fade[str(group)] = results.segments[-1].conditions.energy.bus.battery_modules.lithium_ion_nmc.cell.capacity_fade_factor
            cycle_day[str(group)] = results.segments[-1].conditions.energy.bus.battery_modules.lithium_ion_nmc.cell.cycle_in_day
            save_charge_throughput(charge_throughput, cycle_day, resistance_growth, capacity_fade, storage_dir)
            
            # -------------------------------------------------------------------------------------------    
            # SAVE RESULTS
            # -------------------------------------------------------------------------------------------
            filename = f'e_Twin_Otter_nmc_case_{group}_'
            save_results(results, filename, group, storage_dir)
            create_excel(filename, group, storage_dir)
            
            if plot_mission: 
                Plots.plot_results(results, save_figure_flag=False)       
                        
    tf = time.time() 
    print(f'time taken: {round(((tf - ti) / 60), 3)} mins')
    return 

# ----------------------------------------------------------------------
#   Save Results
# ----------------------------------------------------------------------
def create_excel(filename, group, storage_dir):
    results = load_results(filename, group, storage_dir)
    Create_Excel.write_data(results, filename, group)
    return

def save_results(results, filename, group, storage_dir):
    save_dir = storage_dir
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)  # Create the directory if it doesn't exist
    pickle_file = os.path.join(save_dir, f"{filename}group_number_{group}.pkl")
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return

# ----------------------------------------------------------------------
#   Load Results
# ----------------------------------------------------------------------
def load_charge_throughput(storage_dir, filename='previous_day_data'):
    load_dir = storage_dir
    file_path = os.path.join(load_dir, f"{filename}.pkl")
    if os.path.exists(file_path):
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
            return data['charge_throughput'], data['cycle_day'], data['resistance_growth'], data['capacity_fade']
    else:   
        return {}, {}, {}, {}  # Return empty dictionaries if the file does not exist

def save_charge_throughput(charge_throughput, cycle_day, resistance_growth, capacity_fade, storage_dir, filename='previous_day_data'):
    save_dir = storage_dir
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)  # Create the directory if it doesn't exist
    file_path = os.path.join(save_dir, f"{filename}.pkl")
    data = {
        'charge_throughput': charge_throughput, 
        'cycle_day': cycle_day, 
        'resistance_growth': resistance_growth, 
        'capacity_fade': capacity_fade
    }
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)

def load_results(filename, group, storage_dir):  
    load_dir = storage_dir
    load_file = os.path.join(load_dir, f"{filename}group_number_{group}.pkl")
    if os.path.exists(load_file):
        with open(load_file, 'rb') as file:
            results = pickle.load(file) 
        return results
    else:
        raise FileNotFoundError(f"File {load_file} not found.")

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_arguments()

    # Extract arguments
    HAS_power = args.HAS_power
    HEX_power = args.HEX_power
    RES_dimensions = args.RES_dimensions
    storage_dir = args.storage_dir

    # Print to verify the inputs (optional)
    print(f"HAS_power: {HAS_power}, HEX_power: {HEX_power}, RES_dimensions: {RES_dimensions}")
    print(f"Storage directory: {storage_dir}")

    # Run the degradation simulator
    degradation_simulator(HAS_power, HEX_power, RES_dimensions, storage_dir)