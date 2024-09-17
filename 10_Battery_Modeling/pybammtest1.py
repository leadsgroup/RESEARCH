import pybamm
import numpy as np
import matplotlib.pyplot as plt

def main():
    
    model = pybamm.lithium_ion.DFN()  # Doyle-Fuller-Newman model
    sim = pybamm.Simulation(model)
    
    sim.solve([0, 3600])  # solve for 1 hour
    sim.plot()
    
    return

if __name__ == '__main__': 
    main() 
    plt.show()