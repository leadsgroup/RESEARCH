import cantera as ct
import numpy as np
import pandas as pd
import math
import subprocess

##############################################################
# Input Parameters
##############################################################

reaction_mechanism = 'JetFuelSurrogate.yaml'
fuel = 'N-C12H26'         

T_0 = 2000  #[K]
p_0 = 25*ct.one_atm 
mass_flow_rate = 20  #[kg/s]
composition_0 = 'N-C12H26:1, O2:1.0, N2:3.76'
length = 0.2 #this length is double that of the previously considered legnth
area_0 = 5e-3 #initial cross-sectional area [m**2]
n_steps = 1000

##############################################################
# Chain of Reactors
##############################################################
# The plug flow reactor is represented by a linear chain of zero-dimensional reactors. The gas at the inlet to the first one has the specified inlet composition, and for all others the inlet composition is fixed at the composition of the reactor immediately upstream. Since in a PFR model there is no diffusion, the upstream reactors are not affected by any downstream reactors, and therefore the problem may be solved by simply marching from the first to last reactor, integrating each one to steady state.

# import the gas model and set the initial conditions
gas = ct.Solution(reaction_mechanism)
gas.TPY = T_0, p_0, composition_0
gas2 = gas
dz = length / n_steps
r_vol = area_0 * dz

# create a new reactor
r2 = ct.IdealGasReactor(gas)
r2.volume = r_vol

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
upstream = ct.Reservoir(gas, name='upstream')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
downstream = ct.Reservoir(gas2, name='downstream')

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
m = ct.MassFlowController(upstream, r2, mdot=mass_flow_rate)

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
v = ct.PressureController(r2, downstream, master=m, K=1e-5)
#v = ct.PressureController(r2, downstream, master=m, K=0)

sim2 = ct.ReactorNet([r2])
# define time, space, and other information vectors
z = (np.arange(n_steps) + 1) * dz
u = np.zeros_like(z)  # velocity in each reactor
t_r = np.zeros_like(z)  # residence time in each reactor
t2 = np.zeros_like(z)
T = np.zeros_like(z)
states = ct.SolutionArray(r2.thermo)
x_CO2 = np.zeros_like(z)
x_CO = np.zeros_like(z)
x_H2O = np.zeros_like(z)
y_CO2 = np.zeros_like(z)
y_CO = np.zeros_like(z)
y_H2O = np.zeros_like(z)
EI_CO2 = np.zeros_like(z)
EI_CO = np.zeros_like(z)
EI_H2O = np.zeros_like(z)
mdot_fuel = np.zeros_like(z)

for n in range(n_steps):
    # Set the state of the reservoir to match that of the previous reactor
    gas2.TDY = r2.thermo.TDY
    upstream.syncState() 
    sim2.reinitialize()
    sim2.advance_to_steady_state() 
    t_r[n] = r2.mass / mass_flow_rate  # residence time in this reactor
    t2[n] = np.sum(t_r)
    states.append(r2.thermo.state)
    x_CO2 = states.X[:, gas.species_index('CO2')]
    x_CO = states.X[:, gas.species_index('CO')]
    x_H2O = states.X[:, gas.species_index('H2O')]
    y_CO2 = states.Y[:, gas.species_index('CO2')]
    y_CO = states.Y[:, gas.species_index('CO')]
    y_H2O = states.Y[:, gas.species_index('H2O')]
    u[n] = mass_flow_rate / area_0 / r2.thermo.density
    mdot_fuel[n] = mass_flow_rate * gas.molecular_weights[gas.species_index('N-C12H26')] / gas.density
    
    # CO2
    mdot_CO2_PFR = r2.thermo['CO2'].Y * r2.thermo.density * u[n] * area_0
    EI_CO2[n] = (mdot_CO2_PFR) / mdot_fuel[n] if mdot_fuel[n] > 0 else 0
    
    # CO
    mdot_CO_PFR = r2.thermo['CO'].Y * r2.thermo.density * u[n] * area_0
    EI_CO[n] = (mdot_CO_PFR) / mdot_fuel[n] if mdot_fuel[n] > 0 else 0

    # H2O
    mdot_H2O_PFR = r2.thermo['H2O'].Y * r2.thermo.density * u[n] * area_0
    EI_H2O[n] = (mdot_H2O_PFR) / mdot_fuel[n] if mdot_fuel[n] > 0 else 0      


##############################################################
# Results in matplotlib
##############################################################

import matplotlib.pyplot as plt

plt.figure()
plt.plot(z, x_CO2, label = 'CO2')
plt.plot(z, x_CO, label='CO')
plt.plot(z, x_H2O, label='H2O')
plt.xlabel('$z$ [m]')
plt.ylabel('$X$ [-]')
plt.legend(loc=0)

plt.figure()
plt.plot(z, y_CO2,label = 'CO2')
plt.plot(z, y_CO,  label='CO')
plt.plot(z, y_H2O, label='H2O')
plt.xlabel('$z$ [m]')
plt.ylabel('$Y$ [-]')
plt.legend(loc=0)

plt.figure()
plt.plot(z, EI_CO2, label='EI CO2')
plt.plot(z, EI_CO, label='EI CO')
plt.plot(z, EI_H2O, label='EI H2O')
plt.xlabel('$z$ [m]')
plt.ylabel('$EI$ [-]')
plt.legend(loc=0)

plt.show()