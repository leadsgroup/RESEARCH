#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import pickle

def generate_lhs_samples(variable_limits, num_samples, variable_names,lhs_samples):
    
    # Generate Latin Hypercube Samples
    num_variables = len(variable_limits)
    lhs_samples = np.zeros((num_samples, num_variables))
    for i, (low, high) in enumerate(variable_limits):
        intervals = np.linspace(low, high, num_samples + 1)
        samples = np.random.uniform(intervals[:-1], intervals[1:])
        np.random.shuffle(samples)
        lhs_samples[:, i] = samples
    
    save_results(lhs_samples,"LHS Samples")
    
    return lhs_samples

def save_results(results, filename, storage_dir):
    save_dir = storage_dir
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)  # Create the directory if it doesn't exist
    pickle_file = os.path.join(save_dir, f"{filename}.pkl")
    with open(pickle_file, 'wb') as file:
        pickle.dump(results, file) 
    return
