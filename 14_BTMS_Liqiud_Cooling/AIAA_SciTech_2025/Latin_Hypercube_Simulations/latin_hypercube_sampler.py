#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_lhs_samples(variable_limits, num_samples, variable_names=None):
    
    # Generate Latin Hypercube Samples
    num_variables = len(variable_limits)
    lhs_samples = np.zeros((num_samples, num_variables))
    for i, (low, high) in enumerate(variable_limits):
        intervals = np.linspace(low, high, num_samples + 1)
        samples = np.random.uniform(intervals[:-1], intervals[1:])
        np.random.shuffle(samples)
        lhs_samples[:, i] = samples
    
    # Visualize if there are three variables
    if num_variables == 3:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(lhs_samples[:, 0], lhs_samples[:, 1], lhs_samples[:, 2], 
                              c='blue', alpha=0.7, edgecolors='k')
        
        # Add labels
        if variable_names and len(variable_names) == 3:
            ax.set_xlabel(variable_names[0])
            ax.set_ylabel(variable_names[1])
            ax.set_zlabel(variable_names[2])
        else:
            ax.set_xlabel('Variable 1')
            ax.set_ylabel('Variable 2')
            ax.set_zlabel('Variable 3')
        
        ax.set_title('3D Latin Hypercube Sampling Visualization')
    
    return lhs_samples
