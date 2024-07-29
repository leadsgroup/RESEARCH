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
residence_time_PSR = 0.1  # starting residence time

while combustor_Methane_PSR.T > 500:
    sim_Methane_PSR.initial_time = 0.0  # reset the integrator
    sim_Methane_PSR.advance_to_steady_state()

    # Ensure combustor.T is a float
    temperature_Methane_PSR = combustor_Methane_PSR.T if np.isscalar(combustor_Methane_PSR.T) else combustor_Methane_PSR.T.item()

    # Compute mass flow rates
    mdot_fuel = mdot(0)* Methane_PSR.molecular_weights[Methane_PSR.species_index('CH4')] / Methane_PSR.density
    for species in Methane_PSR.species_names:
        species_index_Methane_PSR = Methane_PSR.species_index(species)
        mdot_species_Methane_PSR = combustor_Methane_PSR.thermo[species].Y * combustor_Methane_PSR.thermo.density * combustor_Methane_PSR.volume / residence_time_PSR
        total_emissions_Methane_PSR[species] += mdot_species_Methane_PSR * residence_time_PSR

    # Compute emission index for CO2
    mdot_CO2 = total_emissions_Methane_PSR['CO2']/residence_time_PSR
    EI_CO2_Methane_PSR = (mdot_CO2*1000) / mdot_fuel if mdot_fuel > 0 else 0
    mdot_CO = total_emissions_Methane_PSR['CO']/residence_time_PSR
    EI_CO_Methane_PSR = (mdot_CO*1000) / mdot_fuel if mdot_fuel > 0 else 0    
    mdot_NO2 = total_emissions_Methane_PSR['NO2']/residence_time_PSR
    EI_NO2_Methane_PSR = (mdot_NO2*1000) / mdot_fuel if mdot_fuel > 0 else 0 
    mdot_NO = total_emissions_Methane_PSR['NO']/residence_time_PSR
    EI_NO_Methane_PSR = (mdot_NO*1000) / mdot_fuel if mdot_fuel > 0 else 0      
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
area_out = 1  # Assuming the area is 1 m^2 for simplification
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
X_0 = final_composition_X
length = 1.5  # *approximate* PFR length [m]
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

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
inlet_Methane_PFR = ct.Reservoir(Methane_PFR)

# create a new reactor
combustor_Methane_PFR = ct.IdealGasReactor(Methane_PFR)
combustor_Methane_PFR.volume = PFR_vol

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
exhaust_Methane_PFR = ct.Reservoir(Methane_PFR)

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
states_Methane_PFR = ct.SolutionArray(combustor_Methane_PFR.thermo)

# Initialize lists for emission indices
EI_CO2_PFR = []
EI_CO_PFR = []
EI_NO2_PFR = []
EI_NO_PFR = []

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
    
# Compute emission indices
mdot_fuel_PFR = mass_flow_rate_PFR * Methane_PFR.molecular_weights[Methane_PFR.species_index('CH4')] / Methane_PFR.density
for species in Methane_PFR.species_names:
    species_index_Methane_PFR = Methane_PFR.species_index(species)
    mdot_species_Methane_PFR = combustor_Methane_PFR.thermo[species].Y * combustor_Methane_PFR.thermo.density * combustor_Methane_PFR.volume / residence_time_PFR
    total_emissions_Methane_PFR[species] += mdot_species_Methane_PFR * residence_time_PFR    
    
mdot_CO2_PFR = total_emissions_Methane_PFR['CO2']
EI_CO2_Methane_PFR = (mdot_CO2_PFR*1000) / mdot_fuel_PFR if mdot_fuel_PFR > 0 else 0

mdot_CO_PFR = total_emissions_Methane_PFR['CO']
EI_CO2_Methane_PFR = (mdot_CO_PFR*1000) / mdot_fuel_PFR if mdot_fuel_PFR > 0 else 0

states_Methane_PFR.append(combustor_Methane_PFR.thermo.state)
#mdot_CO2_PFR = mass_flow_rate_PFR * Methane_PFR['CO2'].Y * Methane_PFR.molecular_weights[Methane_PFR.species_index('CO2')] / Methane_PFR.density
    #EI_CO2_PFR.append((mdot_CO2_PFR * 1000) / mdot_fuel_PFR if mdot_fuel_PFR > 0 else 0)
    
    #mdot_CO_PFR = mass_flow_rate_PFR * Methane_PFR['CO'].Y * Methane_PFR.molecular_weights[Methane_PFR.species_index('CO')] / Methane_PFR.density
    #EI_CO_PFR.append((mdot_CO_PFR * 1000) / mdot_fuel_PFR if mdot_fuel_PFR > 0 else 0)
    
    #mdot_NO2_PFR = mass_flow_rate_PFR * Methane_PFR['NO2'].Y * Methane_PFR.molecular_weights[Methane_PFR.species_index('NO2')] / Methane_PFR.density
    #EI_NO2_PFR.append((mdot_NO2_PFR * 1000) / mdot_fuel_PFR if mdot_fuel_PFR > 0 else 0)
    
    #mdot_NO_PFR = mass_flow_rate_PFR * Methane_PFR['NO'].Y * Methane_PFR.molecular_weights[Methane_PFR.species_index('NO')] / Methane_PFR.density
    #EI_NO_PFR.append((mdot_NO_PFR * 1000) / mdot_fuel_PFR if mdot_fuel_PFR > 0 else 0)    

## Print total emissions for each species
#print("\nMethane total emissions for each species (in kg):")
#for species, total_emission in total_emissions_Methane_PFR.items():
    ## Ensure total_emission is a float
    #total_emission_scalar = total_emission if np.isscalar(total_emission) else total_emission.item()
    #print(f"{species}: {total_emission_scalar:.6e} kg")

# -------------------------------------------------------------------------------------------------------------------------------------------

# Plot results
f, ax1 = plt.subplots(2, 2, figsize=(16, 12))
f.suptitle('CH4 PSR Emissions')
ax1[0, 0].plot(states_Methane_PSR.tres, states_Methane_PSR.EI_CO2_Methane_PSR, '.-', color='C0')
ax1[0, 0].set_xlabel('residence time [s]')
ax1[0, 0].set_title('Emission Index CO2', color='C0')
ax1[0, 0].set_ylabel('EI [g/kg]')
ax1[0, 1].plot(states_Methane_PSR.tres, states_Methane_PSR.EI_CO_Methane_PSR, '.-', color='C1')
ax1[0, 1].set_xlabel('residence time [s]')
ax1[0, 1].set_title('Emission Index CO', color='C1')
ax1[0, 1].set_ylabel('EI [g/kg]')
ax1[1, 0].plot(states_Methane_PSR.tres, states_Methane_PSR.EI_NO2_Methane_PSR, '.-', color='C2')
ax1[1, 0].set_xlabel('residence time [s]')
ax1[1, 0].set_title('Emission Index NO2', color='C2')
ax1[1, 0].set_ylabel('EI [g/kg]')
ax1[1, 1].plot(states_Methane_PSR.tres, states_Methane_PSR.EI_NO_Methane_PSR, '.-', color='C3')
ax1[1, 1].set_xlabel('residence time [s]')
ax1[1, 1].set_title('Emission Index NO', color='C3')
ax1[1, 1].set_ylabel('EI [g/kg]')

## Plotting emission indices
#plt.figure()
#plt.plot(PFR_z, EI_CO2_PFR, label='EI CO2')
#plt.plot(PFR_z, EI_CO_PFR, label='EI CO')
#plt.plot(PFR_z, EI_NO2_PFR, label='EI NO2')
#plt.plot(PFR_z, EI_NO_PFR, label='EI NO')
#plt.xlabel('Distance (m)')
#plt.ylabel('Emission Index (g/kg fuel)')
#plt.legend()
#plt.title('Emission Indices along the PFR')

#plt.figure()
#plt.plot(PFR_z, states_Methane_PFR.T, label='Reactor Chain')
#plt.xlabel('$z$ [m]')
#plt.ylabel('$T$ [K]')
#plt.legend(loc=0)
#plt.figure()

#plt.plot(residence_time_PFR, states_Methane_PFR.X[:, Methane_PFR.species_index('H2')], label='Reactor Chain')
#plt.xlabel('$t$ [s]')
#plt.ylabel('$X_{H_2}$ [-]')
#plt.legend(loc=0)

plt.show()
