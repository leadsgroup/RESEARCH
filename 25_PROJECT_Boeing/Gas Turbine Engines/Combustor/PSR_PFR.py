
#"""
#Created on Wed Jun  5 23:01:04 2019

#@author: ravis
#"""


#import cantera as ct
#import numpy as np
#import pandas as pd
#import math
#import subprocess

###############################################################
## Input Parameters
###############################################################

## input file containing the reaction mechanism

##reaction_mechanism = 'gri30.xml'

##input_cti = 'POLIMI_C1C3_HT_NOX_1412.cti'
##output_yaml = 'POLIMI_C1C3_HT_NOX_1412.yaml'

### Call the cti2yaml conversion tool
##subprocess.run(['cti2yaml', input_cti, output_yaml])

#reaction_mechanism = 'gri30.yaml'

#operation = 'cruise'  # cruise or takeoff
#fuel = 'CH4'          # LH2 or LNG


#T_0 = 1400  #[K]
#p_0 = 13.93 * 100000 #[Pa]
#mass_flow_rate = 0.6589695  #[kg/s]
#composition_0 = 'CH4:1.0, O2:1.0, N2:3.76'

##composition from 2000 reactors and the first 0.001 run
##composition_0 = 'HE:0,AR:0,N2:0.749291,NO:1.52998e-05,N2O:1.39006e-06,O2:0.135015,NO2:5.84178e-08,HNO:3.0538e-10,HNNO:1.21234e-13,HONO:8.94188e-10,HNO2:4.06267e-11,HONO2:6.01211e-12,H2:3.11582e-06,N2H2:2.69884e-15,H2O:0.0512336,H2O2:2.77577e-07,NH3:2.14563e-12,N2H4:8.81081e-22,CO:5.95068e-05,CO2:0.0632162,N:5.97991e-10,O:8.91149e-05,NO3:2.81154e-13,H:2.31571e-07,NH:5.68016e-11,NNH:3.05065e-10,OH:0.00106853,HO2:6.212e-06,NH2:1.25931e-11,H2NO:2.05585e-12,N2H3:1.34003e-18,HCO:1.48393e-12'



##length = 0.17022  # PFR axial length [m] 
##Since the above length results in a very low residence time, it was suggested to increase the
##length in order to compensate for the area correction that is absent in this code, which may be the
##reason for the very low residence time
#length = 0.17022 #this length is double that of the previously considered legnth
#area_0 = 4.514e-3 #initial cross-sectional area [m**2]


## Resolution: The PFR will be simulatHE:0,AR:0,N2:0.749291545,NO:1.44618E-05,N2O:1.21131E-06,O2:0.135460734,NO2:1.60039E-07,HNO:5.02724E-11,HNNO:0,HONO:3.93049E-09,HNO2:2.03227E-11,HONO2:6.12118E-12,H2:2.79905E-07,N2H2:0,H2O:0.05171935,H2O2:1.28107E-07,NH3:1.86471E-10,N2H4:0,CO:2.09032E-05,CO2:0.063276845,N:1.12728E-12,O:7.2018E-06,NO3:7.45695E-13,H:7.68897E-09,NH:1.32867E-12,NNH:5.51473E-12,OH:0.000205489,HO2:1.67872E-06,NH2:2.57024E-11,H2NO:9.66349E-13,N2H3:0,HCO:7.89355E-12ed by 'n_steps' time steps or by a chain
## of 'n_steps' stirred reactors.
#n_steps = 1000



###############################################################
## Chain of Reactors
###############################################################
## The plug flow reactor is represented by a linear chain of zero-dimensional reactors. The gas at the inlet to the first one has the specified inlet composition, and for all others the inlet composition is fixed at the composition of the reactor immediately upstream. Since in a PFR model there is no diffusion, the upstream reactors are not affected by any downstream reactors, and therefore the problem may be solved by simply marching from the first to last reactor, integrating each one to steady state.

## import the gas model and set the initial conditions
#gas = ct.Solution(reaction_mechanism)
#gas.TPY = T_0, p_0, composition_0
#gas2 = gas
#dz = length / n_steps
#r_vol = area_0 * dz

## create a new reactor
#r2 = ct.IdealGasReactor(gas)
#r2.volume = r_vol


## create a reservoir to represent the reactor immediately upstream. Note
## that the gas object is set already to the state of the upstream reactor
#upstream = ct.Reservoir(gas, name='upstream')

## create a reservoir for the reactor to exhaust into. The composition of
## this reservoir is irrelevant.
#downstream = ct.Reservoir(gas2, name='downstream')

## The mass flow rate into the reactor will be fixed by using a
## MassFlowController object.
#m = ct.MassFlowController(upstream, r2, mdot=mass_flow_rate)

## We need an outlet to the downstream reservoir. This will determine the
## pressure in the reactor. The value of K will only affect the transient
## pressure difference.
#v = ct.PressureController(r2, downstream, master=m, K=1e-5)
##v = ct.PressureController(r2, downstream, master=m, K=0)

#sim2 = ct.ReactorNet([r2])
## define time, space, and other information vectors
#z = (np.arange(n_steps) + 1) * dz
#u = np.zeros_like(z)  # velocity in each reactor
#t_r = np.zeros_like(z)  # residence time in each reactor
#t2 = np.zeros_like(z)
#T = np.zeros_like(z)
#area = np.zeros_like(z)
#Vol = np.zeros_like(z)
#states = ct.SolutionArray(r2.thermo)
#y_NO_PFR = np.zeros(100)
#y_HE = np.zeros(100)
#y_AR = np.zeros(100)
#y_N2 = np.zeros(100)
#y_NO_PFR = np.zeros(100)
#y_N2O = np.zeros(100)
#y_O2 = np.zeros(100)
#y_NO2_PFR = np.zeros(100)
#y_HNO = np.zeros(100)
#y_HNNO = np.zeros(100)
#y_HONO = np.zeros(100)
#y_HNO2 = np.zeros(100)
#y_HONO2 = np.zeros(100)
#y_H2 = np.zeros(100)
#y_N2H2 = np.zeros(100)
#y_H2O = np.zeros(100)
#y_H2O2 = np.zeros(100)
#y_NH3 = np.zeros(100)
#y_N2H4 = np.zeros(100)
#y_CO_PFR = np.zeros(100)
#y_CO2_PFR = np.zeros(100)
#y_N = np.zeros(100)
#y_O = np.zeros(100)
#y_NO3 = np.zeros(100)
#y_H = np.zeros(100)
#y_NH = np.zeros(100)
#y_NNH = np.zeros(100)
#y_OH = np.zeros(100)
#y_HO2 = np.zeros(100)
#y_NH2 = np.zeros(100)
#y_H2NO = np.zeros(100)
#y_N2H3 = np.zeros(100)
#y_HCO = np.zeros(100)
#d_0 = np.zeros(100)
#mass = np.zeros_like(z)
#r_rate = np.zeros_like(z)
#r = gas.reaction(100)
#major = gas['CO', 'CO2', 'NO', 'NO2']
#wdot_major = np.zeros_like(z)
#udot_major = np.zeros_like(z)
#vdot_major = np.zeros_like(z)
#for n in range(n_steps):
    ## Set the state of the reservoir to match that of the previous reactor
    #gas2.TDY = r2.thermo.TDY
    #upstream.syncState() #Set the state of the Reactor to match that of the associated ThermoPhase object. 
    ##After calling syncState(), call ReactorNet.reinitialize() before further integration.
    ## integrate the reactor forward in time until steady state is reached
    #sim2.reinitialize()
    #sim2.advance_to_steady_state() #Advance the reactor network in time until steady state is reached.
##The steady state is defined by requiring that the state of the system only changes below a certain
##threshold. The residual is computed using feature scaling
    #Vol[n] = r2.volume
    #t_r[n] = r2.mass / mass_flow_rate  # residence time in this reactor
    #t2[n] = np.sum(t_r)
    #states.append(r2.thermo.state)
    #d_0 = states.DP
    #mass[n] = r2.mass
    #y_AR = states.Y[:, gas.species_index('AR')]
    #y_N2 = states.Y[:, gas.species_index('N2')]
    #y_NO_PFR = states.Y[:, gas.species_index('NO')]
    #y_N2O = states.Y[:, gas.species_index('N2O')]
    #y_O2 = states.Y[:, gas.species_index('O2')]
    #y_NO2_PFR = states.Y[:, gas.species_index('NO2')]
    #y_HNO = states.Y[:, gas.species_index('HNO')]
    #y_H2 = states.Y[:, gas.species_index('H2')]
    #y_CO_PFR = states.Y[:, gas.species_index('CO')]
    #y_CO2_PFR = states.Y[:, gas.species_index('CO2')]
    #y_N = states.Y[:, gas.species_index('N')]
    #y_O = states.Y[:, gas.species_index('O')]
    #y_H = states.Y[:, gas.species_index('H')]
    #y_NH = states.Y[:, gas.species_index('NH')]
    #y_NNH = states.Y[:, gas.species_index('NNH')]
    #y_OH = states.Y[:, gas.species_index('OH')]
    #y_HO2 = states.Y[:, gas.species_index('HO2')]
    #y_NH2 = states.Y[:, gas.species_index('NH2')]
    #rf=gas2.forward_rates_of_progress
    #rr=gas2.reverse_rates_of_progress
    #df = pd.DataFrame(r_rate)
    #rc = gas.forward_rate_constants
    #rd = gas.reverse_rate_constants
    #wdot_major = major.creation_rates
    #udot_major = major.destruction_rates
    #vdot_major = major.net_production_rates
    

###############################################################
## Results in matplotlib
###############################################################

#import matplotlib.pyplot as plt

#plt.figure()
#plt.plot(z, states.T, label='Reactor Chain')
#plt.xlabel('$z$ [m]')
#plt.ylabel('$T$ [K]')
#axes = plt.gca()
#axes.set_ylim([(1300),(T_0 + 100)])
#plt.legend(loc=0)

#plt.figure()
#plt.plot(z, states.P, label='Reactor Chain')
#plt.xlabel('$z$ [m]')
#plt.ylabel('$p$ [Pa]')
#axes = plt.gca()

#plt.figure()
#plt.plot(z, states.Y[:, gas.species_index('CO')], label = 'CO')
#plt.plot(z, states.Y[:, gas.species_index('NO')], label='NO')
#plt.plot(z, states.Y[:, gas.species_index('NO2')], label='NO2')
#plt.plot(z, states.Y[:, gas.species_index('CO2')], label='CO2')
#plt.plot(z, states.Y[:, gas.species_index('HCO')], label='HCO')
#plt.plot(z, states.Y[:, gas.species_index('C2H2')], label='C2H2')
#plt.xlabel('$z$ [m]')
#plt.ylabel('$Y_{H_2}$ [-]')
#plt.legend(loc=0)

#plt.figure()
#plt.plot(z, states.Y[:, gas.species_index('CO2')], label='Reactor Chain')
#plt.xlabel('$z$ [m]')
#plt.ylabel('$Y_{CO_2}$ [-]')
#plt.legend(loc=0)

#plt.figure()
#plt.plot(z, states.Y[:, gas.species_index('CO')], label='Reactor Chain')
#plt.xlabel('$z$ [m]')
#plt.ylabel('$Y_{CO}$ [-]')
#plt.legend(loc=0)

#plt.figure()
#plt.plot(z, states.Y[:, gas.species_index('HCO')], label='Reactor Chain')
#plt.xlabel('$z$ [m]')
#plt.ylabel('$Y_{HCO}$ [-]')
#plt.legend(loc=0)

#plt.figure()
#plt.plot(z, states.Y[:, gas.species_index('NO')], label='Reactor Chain')
#plt.xlabel('$z$ [m]')
#plt.ylabel('$Y_{NO}$ [-]')
#plt.legend(loc=0)

#plt.figure()
#plt.plot(z, states.Y[:, gas.species_index('NO2')], label='Reactor Chain')
#plt.xlabel('$z$ [ms]')
#plt.ylabel('$Y_{NO2}$ [-]')
#plt.legend(loc=0)
#plt.show()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
Perfectly Stirred Reactor (PSR) + Plug Flow Reactor (PFR)
=========================================================

This code solves a PSR + PFR problem.
The PFR is computed as a chain of reactors.
Methane, Ethane and Propane emissions are computed.

"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------------------------
# Methane
# -------------------------------------------------------------------------------------------

#######################################################################
# Perfectly Stirred Reactor (PSR)
#######################################################################

# Use reaction mechanism GRI-Mech 3.0. For 0-D simulations, no transport model is necessary.
Methane_PSR = ct.Solution('gri30.yaml')

# Create a Reservoir for the inlet, set to a methane/air mixture at a specified equivalence ratio
equiv_ratio = 0.5  # lean combustion
Methane_PSR.TP = 300.0, ct.one_atm
Methane_PSR.set_equivalence_ratio(equiv_ratio, 'CH4:1.0', 'O2:1.0, N2:3.76')
inlet_Methane_PSR = ct.Reservoir(Methane_PSR)

# Create the combustor, and fill it initially with a mixture consisting of the equilibrium products of the inlet mixture.
Methane_PSR.equilibrate('HP')
combustor_Methane_PSR = ct.IdealGasReactor(Methane_PSR)
combustor_Methane_PSR.volume = 1.0

# Create a reservoir for the exhaust
exhaust_Methane_PSR = ct.Reservoir(Methane_PSR)

# Use a variable mass flow rate to keep the residence time in the reactor constant.
def mdot(t):
    return combustor_Methane_PSR.mass / residence_time_PSR

inlet_mfc_Methane_PSR = ct.MassFlowController(inlet_Methane_PSR, combustor_Methane_PSR, mdot=mdot)

# A PressureController has a baseline mass flow rate matching the 'primary' MassFlowController.
outlet_mfc_Methane_PSR = ct.PressureController(combustor_Methane_PSR, exhaust_Methane_PSR, primary=inlet_mfc_Methane_PSR, K=0.01)

# The simulation only contains one reactor
sim_Methane_PSR = ct.ReactorNet([combustor_Methane_PSR])

# Initialize total mass emission dictionary
total_emissions_Methane_PSR = {species: 0.0 for species in Methane_PSR.species_names}

# Run a loop over decreasing residence times, until the reactor is extinguished.
states_Methane_PSR = ct.SolutionArray(Methane_PSR, extra=['tres', 'EI_CO2_Methane_PSR', 'EI_CO_Methane_PSR', 'EI_NO2_Methane_PSR', 'EI_NO_Methane_PSR'])
#residence_time_PFR_n[n] = combustor_Methane_PFR.mass / mass_flow_rate_PFR
residence_time_PSR = 0.1  # starting residence time

while combustor_Methane_PSR.T > 500:
    sim_Methane_PSR.initial_time = 0.0  # reset the integrator
    sim_Methane_PSR.advance_to_steady_state()

    # Ensure combustor.T is a float
    temperature_Methane_PSR = combustor_Methane_PSR.T if np.isscalar(combustor_Methane_PSR.T) else combustor_Methane_PSR.T.item()

    # Compute mass flow rates
    # mdot_fuel = mdot(0)* Methane_PSR.molecular_weights[Methane_PSR.species_index('CH4')] / Methane_PSR.density
    mdot_fuel = mdot(0)* Methane_PSR.Y[Methane_PSR.species_index('CH4')]
    for species in Methane_PSR.species_names:
        species_index_Methane_PSR = Methane_PSR.species_index(species)
        #mdot_species_Methane_PSR = combustor_Methane_PSR.thermo[species].Y * combustor_Methane_PSR.thermo.density * combustor_Methane_PSR.volume / residence_time_PSR
        mdot_species_Methane_PSR = combustor_Methane_PSR.thermo[species].Y * mdot_fuel
        total_emissions_Methane_PSR[species] += mdot_species_Methane_PSR * residence_time_PSR

    # Compute emission index for CO2
    mdot_CO2 = total_emissions_Methane_PSR['CO2']/residence_time_PSR
    EI_CO2_Methane_PSR = (mdot_CO2) / mdot_fuel if mdot_fuel > 0 else 0
    mdot_CO = total_emissions_Methane_PSR['CO']/residence_time_PSR
    EI_CO_Methane_PSR = (mdot_CO) / mdot_fuel if mdot_fuel > 0 else 0    
    mdot_NO2 = total_emissions_Methane_PSR['NO2']/residence_time_PSR
    EI_NO2_Methane_PSR = (mdot_NO2) / mdot_fuel if mdot_fuel > 0 else 0 
    mdot_NO = total_emissions_Methane_PSR['NO']/residence_time_PSR
    EI_NO_Methane_PSR = (mdot_NO) / mdot_fuel if mdot_fuel > 0 else 0      
    states_Methane_PSR.append(combustor_Methane_PSR.thermo.state, tres=residence_time_PSR, EI_CO2_Methane_PSR=EI_CO2_Methane_PSR, EI_CO_Methane_PSR=EI_CO_Methane_PSR, EI_NO2_Methane_PSR=EI_NO2_Methane_PSR, EI_NO_Methane_PSR=EI_NO_Methane_PSR)
    residence_time_PSR *= 0.9  # decrease the residence time for the next iteration

# Print total emissions for each species
print("\nMethane PSR emissions for each species (in kg):")
for species, total_emission in total_emissions_Methane_PSR.items():
    # Ensure total_emission is a float
    total_emission_scalar = total_emission if np.isscalar(total_emission) else total_emission.item()
    print(f"{species}: {total_emission_scalar:.6e} kg")

# Calculate outflow velocity from PSR
mass_flow_rate_out = mdot(0)  # mass flow rate out of the combustor
density_out = combustor_Methane_PSR.thermo.density  # density of the gas in the combustor
area_out = 0.02  # Assuming the area is 1 m^2 for simplification
outflow_velocity = mass_flow_rate_out / (density_out * area_out)

# Extract the final state of the methane combustor
final_state = combustor_Methane_PSR.thermo.state
final_temperature = combustor_Methane_PSR.T
final_pressure = combustor_Methane_PSR.thermo.P
final_composition_X = combustor_Methane_PSR.thermo.X
final_composition_Y = combustor_Methane_PSR.thermo.Y

#######################################################################
# Plug Flow Reactor (PFR)
#######################################################################

T_0 = final_temperature  # inlet temperature [K]
P_0 = final_pressure  # constant pressure [Pa]
#T_0 = 1400  #[K]
#P_0 = 13.93 * 100000 #[Pa]
X_0 = final_composition_X
length = 0.15  # *approximate* PFR length [m]
u_0 = outflow_velocity  # inflow velocity [m/s]
PFR_area = area_out  # cross-sectional area [m**2]

# Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
# of 'n_steps' stirred reactors.
n_steps = 1000

# The plug flow reactor is represented by a linear chain of zero-dimensional
# reactors. The gas at the inlet to the first one has the specified inlet
# composition, and for all others the inlet composition is fixed at the
# composition of the reactor immediately upstream. Since in a PFR model there
# is no diffusion, the upstream reactors are not affected by any downstream
# reactors, and therefore the problem may be solved by simply marching from
# the first to last reactor, integrating each one to steady state.

# import the gas model and set the initial conditions
Methane_PFR = ct.Solution('gri30.yaml')
Methane_PFR.TPX = T_0, P_0, X_0 
mass_flow_rate_PFR = u_0 * Methane_PFR.density * PFR_area
dz = length / n_steps
PFR_vol = PFR_area * dz

# create a new reactor
combustor_Methane_PFR = ct.IdealGasReactor(Methane_PFR)
combustor_Methane_PFR.volume = PFR_vol

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
inlet_Methane_PFR = ct.Reservoir(Methane_PFR,name='inlet')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
exhaust_Methane_PFR = ct.Reservoir(Methane_PFR,name='exhaust')

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
inlet_mfc_Methane_PFR = ct.MassFlowController(inlet_Methane_PFR, combustor_Methane_PFR, mdot=mass_flow_rate_PFR)

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
outlet_mfc_Methane_PFR = ct.PressureController(combustor_Methane_PFR, exhaust_Methane_PFR, primary=inlet_mfc_Methane_PFR, K=1e-5)

sim_Methane_PFR = ct.ReactorNet([combustor_Methane_PFR])

# Initialize total mass emission dictionary
total_emissions_Methane_PFR = {species: 0.0 for species in Methane_PFR.species_names}

# define time, space, and other information vectors
PFR_z = (np.arange(n_steps) + 1) * dz
residence_time_PFR_n = np.zeros_like(PFR_z)  # residence time in each reactor
u_PFR = np.zeros_like(PFR_z)
residence_time_PFR = np.zeros_like(PFR_z)
mdot_fuel_PFR = np.zeros_like(PFR_z)
states_Methane_PFR = ct.SolutionArray(combustor_Methane_PFR.thermo)
x_NO2_PFR = np.zeros_like(PFR_z)
x_NO_PFR = np.zeros_like(PFR_z)
x_CO2_PFR = np.zeros_like(PFR_z)
x_CO_PFR = np.zeros_like(PFR_z)
y_NO2_PFR = np.zeros_like(PFR_z)
y_NO_PFR = np.zeros_like(PFR_z)
y_CO2_PFR = np.zeros_like(PFR_z)
y_CO_PFR = np.zeros_like(PFR_z)
EI_CO2_PFR = np.zeros_like(PFR_z)
EI_CO_PFR = np.zeros_like(PFR_z)
EI_NO2_PFR = np.zeros_like(PFR_z)
EI_NO_PFR = np.zeros_like(PFR_z)

# iterate through the PFR cells
for n in range(n_steps):
    # Set the state of the reservoir to match that of the previous reactor
    Methane_PFR.TDY = combustor_Methane_PFR.thermo.TDY
    inlet_Methane_PFR.syncState()
    # integrate the reactor forward in time until steady state is reached
    sim_Methane_PFR.reinitialize()
    sim_Methane_PFR.advance_to_steady_state()
    # compute velocity and transform into time
    u_PFR[n] = mass_flow_rate_PFR / PFR_area / combustor_Methane_PFR.thermo.density
    residence_time_PFR_n[n] = combustor_Methane_PFR.mass / mass_flow_rate_PFR  # residence time in this reactor
    residence_time_PFR[n] = np.sum(residence_time_PFR_n)
    # write output data
    states_Methane_PFR.append(combustor_Methane_PFR.thermo.state)
    x_NO2_PFR = states_Methane_PFR.X[:, Methane_PFR.species_index('NO2')]
    x_NO_PFR = states_Methane_PFR.X[:, Methane_PFR.species_index('NO')]
    x_CO2_PFR = states_Methane_PFR.X[:, Methane_PFR.species_index('CO2')]
    x_CO_PFR = states_Methane_PFR.X[:, Methane_PFR.species_index('CO')]    
    y_NO2_PFR = states_Methane_PFR.Y[:, Methane_PFR.species_index('NO2')]
    y_NO_PFR = states_Methane_PFR.Y[:, Methane_PFR.species_index('NO')]
    y_CO2_PFR = states_Methane_PFR.Y[:, Methane_PFR.species_index('CO2')]
    y_CO_PFR = states_Methane_PFR.Y[:, Methane_PFR.species_index('CO')] 
    
    # Compute emission indices for CO2, CO, NO2, and NO
    mdot_fuel_PFR[n] = mass_flow_rate_PFR * Methane_PFR.molecular_weights[Methane_PFR.species_index('CH4')] / Methane_PFR.density
    
    # NO2
    mdot_NO2_PFR = combustor_Methane_PFR.thermo['NO2'].Y * combustor_Methane_PFR.thermo.density * u_PFR[n] * PFR_area
    EI_NO2_PFR[n] = (mdot_NO2_PFR) / mdot_fuel_PFR[n] if mdot_fuel_PFR[n] > 0 else 0
    
    # NO
    mdot_NO_PFR = combustor_Methane_PFR.thermo['NO'].Y * combustor_Methane_PFR.thermo.density * u_PFR[n] * PFR_area
    EI_NO_PFR[n] = (mdot_NO_PFR) / mdot_fuel_PFR[n] if mdot_fuel_PFR[n] > 0 else 0
    
    # CO2
    mdot_CO2_PFR = combustor_Methane_PFR.thermo['CO2'].Y * combustor_Methane_PFR.thermo.density * u_PFR[n] * PFR_area
    EI_CO2_PFR[n] = (mdot_CO2_PFR) / mdot_fuel_PFR[n] if mdot_fuel_PFR[n] > 0 else 0
    
    # CO
    mdot_CO_PFR = combustor_Methane_PFR.thermo['CO'].Y * combustor_Methane_PFR.thermo.density * u_PFR[n] * PFR_area
    EI_CO_PFR[n] = (mdot_CO_PFR) / mdot_fuel_PFR[n] if mdot_fuel_PFR[n] > 0 else 0
        

## Print total emissions for each species from the PFR
#print("\nMethane PFR emissions for each species (in kg):")
#for species_PFR, total_emission_PFR in total_emissions_Methane_PFR.items():
    #total_emission_scalar_PFR = np.sum(total_emission_PFR) if isinstance(total_emission_PFR, np.ndarray) else total_emission_PFR
    #print(f"{species_PFR}: {total_emission_scalar_PFR:.6e} kg")
#print("\nMethane emissions for each species (in kg):")
#for species, total_emission in total_emissions_Methane_PSR.items():
    #for species_PFR, total_emission_PFR in total_emissions_Methane_PFR.items():
        #total_emission_scalar = np.sum(total_emission) if isinstance(total_emission, np.ndarray) else total_emission
        #total_emission_scalar_PFR = np.sum(total_emission_PFR) if isinstance(total_emission_PFR, np.ndarray) else total_emission_PFR
    #print(f"{species}: {total_emission_scalar + total_emission_scalar_PFR:.6e} kg")

# Plot results
f, ax1 = plt.subplots(2, 2, figsize=(16, 12))
f.suptitle('CH4 PSR Emissions')
ax1[0, 0].plot(states_Methane_PSR.tres, states_Methane_PSR.EI_CO2_Methane_PSR, '.-', color='C0')
ax1[0, 0].set_xlabel('residence time [s]')
ax1[0, 0].set_title('Emission Index CO2', color='C0')
ax1[0, 0].set_ylabel('EI [kg/kg]')
ax1[0, 1].plot(states_Methane_PSR.tres, states_Methane_PSR.EI_CO_Methane_PSR, '.-', color='C1')
ax1[0, 1].set_xlabel('residence time [s]')
ax1[0, 1].set_title('Emission Index CO', color='C1')
ax1[0, 1].set_ylabel('EI [kg/kg]')
ax1[1, 0].plot(states_Methane_PSR.tres, states_Methane_PSR.EI_NO2_Methane_PSR, '.-', color='C2')
ax1[1, 0].set_xlabel('residence time [s]')
ax1[1, 0].set_title('Emission Index NO2', color='C2')
ax1[1, 0].set_ylabel('EI [kg/kg]')
ax1[1, 1].plot(states_Methane_PSR.tres, states_Methane_PSR.EI_NO_Methane_PSR, '.-', color='C3')
ax1[1, 1].set_xlabel('residence time [s]')
ax1[1, 1].set_title('Emission Index NO', color='C3')
ax1[1, 1].set_ylabel('EI [kg/kg]')

#Plotting emission indices
plt.figure()
plt.plot(PFR_z, EI_CO2_PFR, label='EI CO')
plt.plot(PFR_z, EI_CO_PFR, label='EI CO')
plt.plot(PFR_z, EI_NO2_PFR, label='EI NO2')
plt.plot(PFR_z, EI_NO_PFR, label='EI NO')
plt.xlabel('Distance (m)')
plt.ylabel('Emission Index (g/kg fuel)')
plt.legend()
plt.title('Emission Indices along the PFR')

plt.figure()
plt.plot(PFR_z, states_Methane_PFR.T, label='Reactor Chain')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)

plt.figure()
plt.plot(PFR_z, states_Methane_PFR.P, label='Reactor Chain')
plt.xlabel('$z$ [m]')
plt.ylabel('$p$ [Pa]')

plt.figure()
plt.plot(PFR_z, states_Methane_PFR.X[:, Methane_PFR.species_index('CO')], label = 'CO')
plt.plot(PFR_z, states_Methane_PFR.X[:, Methane_PFR.species_index('NO')], label='NO')
plt.plot(PFR_z, states_Methane_PFR.X[:, Methane_PFR.species_index('NO2')], label='NO2')
plt.plot(PFR_z, states_Methane_PFR.X[:, Methane_PFR.species_index('CO2')], label='CO2')
plt.xlabel('$z$ [m]')
plt.ylabel('$Molar Fractions (X)$ [-]')
plt.legend(loc=0)

plt.figure()
plt.plot(PFR_z, states_Methane_PFR.Y[:, Methane_PFR.species_index('CO')], label = 'CO')
plt.plot(PFR_z, states_Methane_PFR.Y[:, Methane_PFR.species_index('NO')], label='NO')
plt.plot(PFR_z, states_Methane_PFR.Y[:, Methane_PFR.species_index('NO2')], label='NO2')
plt.plot(PFR_z, states_Methane_PFR.Y[:, Methane_PFR.species_index('CO2')], label='CO2')
plt.xlabel('$z$ [m]')
plt.ylabel('$Mass Fractions (Y)$ [-]')
plt.legend(loc=0)

plt.show()
