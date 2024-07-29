"""
Combustor residence time
========================

Calculate steady-state solutions for a combustor, modeled as a single well-stirred
reactor, for different residence times.

We are interested in the steady-state burning solution. This example explores
the effect of changing the residence time on completeness of reaction (through
the burned gas temperature) and on the total heat release rate.

Demonstrates the use of a `MassFlowController` where the mass flow rate function
depends on variables other than time by capturing these variables from the
enclosing scope. Also shows the use of a `PressureController` to create a constant
pressure reactor with a fixed volume.

Requires: cantera >= 3.0, matplotlib >= 2.0

.. tags:: Python, combustion, reactor network, well-stirred reactor, plotting
"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# -------------------------------------------------------------------------------------------
# Methane
# -------------------------------------------------------------------------------------------

# Use reaction mechanism GRI-Mech 3.0. For 0-D simulations, no transport model is necessary.
Methane = ct.Solution('gri30.yaml')

# Create a Reservoir for the inlet, set to a methane/air mixture at a specified equivalence ratio
equiv_ratio = 0.5  # lean combustion
Methane.TP = 300.0, ct.one_atm
Methane.set_equivalence_ratio(equiv_ratio, 'CH4:1.0', 'O2:1.0, N2:3.76')
inlet_Methane = ct.Reservoir(Methane)

# Create the combustor, and fill it initially with a mixture consisting of the equilibrium products of the inlet mixture.
Methane.equilibrate('HP')
combustor_Methane = ct.IdealGasReactor(Methane)
combustor_Methane.volume = 1.0

# Create a reservoir for the exhaust
exhaust_Methane = ct.Reservoir(Methane)

# Use a variable mass flow rate to keep the residence time in the reactor constant.
def mdot(t):
    return combustor_Methane.mass / residence_time

inlet_mfc_Methane = ct.MassFlowController(inlet_Methane, combustor_Methane, mdot=mdot)

# A PressureController has a baseline mass flow rate matching the 'primary' MassFlowController.
outlet_mfc_Methane = ct.PressureController(combustor_Methane, exhaust_Methane, primary=inlet_mfc_Methane, K=0.01)

# The simulation only contains one reactor
sim_Methane = ct.ReactorNet([combustor_Methane])

# Initialize total mass emission dictionary
total_emissions_Methane = {species: 0.0 for species in Methane.species_names}

# Run a loop over decreasing residence times, until the reactor is extinguished.
states_Methane = ct.SolutionArray(Methane, extra=['tres', 'EI_CO2_Methane', 'EI_CO_Methane', 'EI_NO2_Methane', 'EI_NO_Methane'])
residence_time = 0.1  # starting residence time

while combustor_Methane.T > 500:
    sim_Methane.initial_time = 0.0  # reset the integrator
    sim_Methane.advance_to_steady_state()

    # Ensure combustor.T is a float
    temperature_Methane = combustor_Methane.T if np.isscalar(combustor_Methane.T) else combustor_Methane.T.item()
    # print('tres = {:.2e}; T = {:.1f}'.format(residence_time, temperature_Methane))

    # Compute mass flow rates
    mdot_fuel = mdot(0)* Methane.molecular_weights[Methane.species_index('CH4')] / Methane.density
    for species in Methane.species_names:
        species_index_Methane = Methane.species_index(species)
        mdot_species_Methane = combustor_Methane.thermo[species].Y * combustor_Methane.thermo.density * combustor_Methane.volume / residence_time
        total_emissions_Methane[species] += mdot_species_Methane * residence_time

    # Compute emission index for CO2
    mdot_CO2 = total_emissions_Methane['CO2']/residence_time
    EI_CO2_Methane = (mdot_CO2*1000) / mdot_fuel if mdot_fuel > 0 else 0
    mdot_CO = total_emissions_Methane['CO']/residence_time
    EI_CO_Methane = (mdot_CO*1000) / mdot_fuel if mdot_fuel > 0 else 0    
    mdot_NO2 = total_emissions_Methane['NO2']/residence_time
    EI_NO2_Methane = (mdot_NO2*1000) / mdot_fuel if mdot_fuel > 0 else 0 
    mdot_NO = total_emissions_Methane['NO']/residence_time
    EI_NO_Methane = (mdot_NO*1000) / mdot_fuel if mdot_fuel > 0 else 0      
    states_Methane.append(combustor_Methane.thermo.state, tres=residence_time, EI_CO2_Methane=EI_CO2_Methane, EI_CO_Methane=EI_CO_Methane, EI_NO2_Methane=EI_NO2_Methane, EI_NO_Methane=EI_NO_Methane)
    residence_time *= 0.9  # decrease the residence time for the next iteration

# Print total emissions for each species
print("\nMethane total emissions for each species (in kg):")
for species, total_emission in total_emissions_Methane.items():
    # Ensure total_emission is a float
    total_emission_scalar = total_emission if np.isscalar(total_emission) else total_emission.item()
    print(f"{species}: {total_emission_scalar:.6e} kg")

# Plot results
f, ax1 = plt.subplots(2, 2, figsize=(16, 12))
f.suptitle('CH4 Emissions')
ax1[0, 0].plot(states_Methane.tres, states_Methane.EI_CO2_Methane, '.-', color='C0')
ax1[0, 0].set_xlabel('residence time [s]')
ax1[0, 0].set_title('Emission Index CO2', color='C0')
ax1[0, 0].set_ylabel('EI [g/kg]')
ax1[0, 1].plot(states_Methane.tres, states_Methane.EI_CO_Methane, '.-', color='C1')
ax1[0, 1].set_xlabel('residence time [s]')
ax1[0, 1].set_title('Emission Index CO', color='C1')
ax1[0, 1].set_ylabel('EI [g/kg]')
ax1[1, 0].plot(states_Methane.tres, states_Methane.EI_NO2_Methane, '.-', color='C2')
ax1[1, 0].set_xlabel('residence time [s]')
ax1[1, 0].set_title('Emission Index NO2', color='C2')
ax1[1, 0].set_ylabel('EI [g/kg]')
ax1[1, 1].plot(states_Methane.tres, states_Methane.EI_NO_Methane, '.-', color='C3')
ax1[1, 1].set_xlabel('residence time [s]')
ax1[1, 1].set_title('Emission Index NO', color='C3')
ax1[1, 1].set_ylabel('EI [g/kg]')

# -------------------------------------------------------------------------------------------
# Ethane
# -------------------------------------------------------------------------------------------

# Use reaction mechanism GRI-Mech 3.0. For 0-D simulations, no transport model is necessary.
Ethane = ct.Solution('gri30.yaml')

# Create a Reservoir for the inlet, set to a methane/air mixture at a specified equivalence ratio
equiv_ratio = 0.5  # lean combustion
Ethane.TP = 300.0, ct.one_atm
Ethane.set_equivalence_ratio(equiv_ratio, 'C2H6:1.0', 'O2:1.0, N2:3.76')
inlet_Ethane = ct.Reservoir(Ethane)

# Create the combustor, and fill it initially with a mixture consisting of the equilibrium products of the inlet mixture.
Ethane.equilibrate('HP')
combustor_Ethane = ct.IdealGasReactor(Ethane)
combustor_Ethane.volume = 1.0

# Create a reservoir for the exhaust
exhaust_Ethane = ct.Reservoir(Ethane)

# Use a variable mass flow rate to keep the residence time in the reactor constant.
def mdot(t):
    return combustor_Ethane.mass / residence_time

inlet_mfc_Ethane = ct.MassFlowController(inlet_Ethane, combustor_Ethane, mdot=mdot)

# A PressureController has a baseline mass flow rate matching the 'primary' MassFlowController.
outlet_mfc_Ethane = ct.PressureController(combustor_Ethane, exhaust_Ethane, primary=inlet_mfc_Ethane, K=0.01)

# The simulation only contains one reactor
sim_Ethane = ct.ReactorNet([combustor_Ethane])

# Initialize total mass emission dictionary
total_emissions_Ethane = {species: 0.0 for species in Ethane.species_names}

# Run a loop over decreasing residence times, until the reactor is extinguished.
states_Ethane = ct.SolutionArray(Ethane, extra=['tres', 'EI_CO2_Ethane', 'EI_CO_Ethane', 'EI_NO2_Ethane', 'EI_NO_Ethane'])
residence_time = 0.01  # starting residence time
while combustor_Ethane.T > 500:
    sim_Ethane.initial_time = 0.0  # reset the integrator
    sim_Ethane.advance_to_steady_state()

    # Ensure combustor.T is a float
    temperature_Ethane = combustor_Ethane.T if np.isscalar(combustor_Ethane.T) else combustor_Ethane.T.item()
    # print('tres = {:.2e}; T = {:.1f}'.format(residence_time, temperature_Ethane))

    # Compute mass flow rates
    mdot_fuel = mdot(0)* Ethane.molecular_weights[Ethane.species_index('C2H6')] / Ethane.density
    for species in Ethane.species_names:
        species_index_Ethane = Ethane.species_index(species)
        mdot_species_Ethane = combustor_Ethane.thermo[species].Y * combustor_Ethane.thermo.density * combustor_Ethane.volume / residence_time
        total_emissions_Ethane[species] += mdot_species_Ethane * residence_time

    # Compute emission index for CO2
    mdot_CO2 = total_emissions_Ethane['CO2']/residence_time
    EI_CO2_Ethane = (mdot_CO2*1000) / mdot_fuel if mdot_fuel > 0 else 0
    mdot_CO = total_emissions_Ethane['CO']/residence_time
    EI_CO_Ethane = (mdot_CO*1000) / mdot_fuel if mdot_fuel > 0 else 0    
    mdot_NO2 = total_emissions_Ethane['NO2']/residence_time
    EI_NO2_Ethane = (mdot_NO2*1000) / mdot_fuel if mdot_fuel > 0 else 0 
    mdot_NO = total_emissions_Ethane['NO']/residence_time
    EI_NO_Ethane = (mdot_NO*1000) / mdot_fuel if mdot_fuel > 0 else 0      
    states_Ethane.append(combustor_Ethane.thermo.state, tres=residence_time, EI_CO2_Ethane=EI_CO2_Ethane, EI_CO_Ethane=EI_CO_Ethane, EI_NO2_Ethane=EI_NO2_Ethane, EI_NO_Ethane=EI_NO_Ethane)
    residence_time *= 0.9  # decrease the residence time for the next iteration

# Print total emissions for each species
print("\nEthane total emissions for each species (in kg):")
for species, total_emission in total_emissions_Ethane.items():
    # Ensure total_emission is a float
    total_emission_scalar = total_emission if np.isscalar(total_emission) else total_emission.item()
    print(f"{species}: {total_emission_scalar:.6e} kg")





# Plot results
f, ax1 = plt.subplots(2, 2, figsize=(16, 12))
f.suptitle('C2H6 Emissions')
ax1[0, 0].plot(states_Ethane.tres, states_Ethane.EI_CO2_Ethane, '.-', color='C0')
ax1[0, 0].set_xlabel('residence time [s]')
ax1[0, 0].set_title('Emission Index CO2', color='C0')
ax1[0, 0].set_ylabel('EI [g/kg]')
ax1[0, 1].plot(states_Ethane.tres, states_Ethane.EI_CO_Ethane, '.-', color='C1')
ax1[0, 1].set_xlabel('residence time [s]')
ax1[0, 1].set_title('Emission Index CO', color='C1')
ax1[0, 1].set_ylabel('EI [g/kg]')
ax1[1, 0].plot(states_Ethane.tres, states_Ethane.EI_NO2_Ethane, '.-', color='C2')
ax1[1, 0].set_xlabel('residence time [s]')
ax1[1, 0].set_title('Emission Index NO2', color='C2')
ax1[1, 0].set_ylabel('EI [g/kg]')
ax1[1, 1].plot(states_Ethane.tres, states_Ethane.EI_NO_Ethane, '.-', color='C3')
ax1[1, 1].set_xlabel('residence time [s]')
ax1[1, 1].set_title('Emission Index NO', color='C3')
ax1[1, 1].set_ylabel('EI [g/kg]')

# -------------------------------------------------------------------------------------------
# Propane
# -------------------------------------------------------------------------------------------

# Use reaction mechanism GRI-Mech 3.0. For 0-D simulations, no transport model is necessary.
Propane = ct.Solution('gri30.yaml')

# Create a Reservoir for the inlet, set to a methane/air mixture at a specified equivalence ratio
equiv_ratio = 0.5  # lean combustion
Propane.TP = 300.0, ct.one_atm
Propane.set_equivalence_ratio(equiv_ratio, 'C3H8:1.0', 'O2:1.0, N2:3.76')
inlet_Propane = ct.Reservoir(Propane)

# Create the combustor, and fill it initially with a mixture consisting of the equilibrium products of the inlet mixture.
Propane.equilibrate('HP')
combustor_Propane = ct.IdealGasReactor(Propane)
combustor_Propane.volume = 1.0

# Create a reservoir for the exhaust
exhaust_Propane = ct.Reservoir(Propane)

# Use a variable mass flow rate to keep the residence time in the reactor constant.
def mdot(t):
    return combustor_Propane.mass / residence_time

inlet_mfc_Propane = ct.MassFlowController(inlet_Propane, combustor_Propane, mdot=mdot)

# A PressureController has a baseline mass flow rate matching the 'primary' MassFlowController.
outlet_mfc_Propane = ct.PressureController(combustor_Propane, exhaust_Propane, primary=inlet_mfc_Propane, K=0.01)

# The simulation only contains one reactor
sim_Propane = ct.ReactorNet([combustor_Propane])

# Initialize total mass emission dictionary
total_emissions_Propane = {species: 0.0 for species in Propane.species_names}

# Run a loop over decreasing residence times, until the reactor is extinguished.
states_Propane = ct.SolutionArray(Propane, extra=['tres', 'EI_CO2_Propane', 'EI_CO_Propane', 'EI_NO2_Propane', 'EI_NO_Propane'])
residence_time = 0.1  # starting residence time
while combustor_Propane.T > 500:
    sim_Propane.initial_time = 0.0  # reset the integrator
    sim_Propane.advance_to_steady_state()

    # Ensure combustor.T is a float
    temperature_Propane = combustor_Propane.T if np.isscalar(combustor_Propane.T) else combustor_Propane.T.item()
    # print('tres = {:.2e}; T = {:.1f}'.format(residence_time, temperature_Propane))

    # Compute mass flow rates
    mdot_fuel = mdot(0)* Propane.molecular_weights[Propane.species_index('C3H8')] / Propane.density
    for species in Propane.species_names:
        species_index_Propane = Propane.species_index(species)
        mdot_species_Propane = combustor_Propane.thermo[species].Y * combustor_Propane.thermo.density * combustor_Propane.volume / residence_time
        total_emissions_Propane[species] += mdot_species_Propane * residence_time

    # Compute emission index for CO2
    mdot_CO2 = total_emissions_Propane['CO2']/residence_time
    EI_CO2_Propane = (mdot_CO2*1000) / mdot_fuel if mdot_fuel > 0 else 0
    mdot_CO = total_emissions_Propane['CO']/residence_time
    EI_CO_Propane = (mdot_CO*1000) / mdot_fuel if mdot_fuel > 0 else 0    
    mdot_NO2 = total_emissions_Propane['NO2']/residence_time
    EI_NO2_Propane = (mdot_NO2*1000) / mdot_fuel if mdot_fuel > 0 else 0 
    mdot_NO = total_emissions_Propane['NO']/residence_time
    EI_NO_Propane = (mdot_NO*1000) / mdot_fuel if mdot_fuel > 0 else 0      
    states_Propane.append(combustor_Propane.thermo.state, tres=residence_time, EI_CO2_Propane=EI_CO2_Propane, EI_CO_Propane=EI_CO_Propane, EI_NO2_Propane=EI_NO2_Propane, EI_NO_Propane=EI_NO_Propane)
    residence_time *= 0.9  # decrease the residence time for the next iteration

# Print total emissions for each species
print("\nPropane total emissions for each species (in kg):")
for species, total_emission in total_emissions_Propane.items():
    # Ensure total_emission is a float
    total_emission_scalar = total_emission if np.isscalar(total_emission) else total_emission.item()
    print(f"{species}: {total_emission_scalar:.6e} kg")

# Plot results
f, ax1 = plt.subplots(2, 2, figsize=(16, 12))
f.suptitle('C3H8 Emissions')
ax1[0, 0].plot(states_Propane.tres, states_Propane.EI_CO2_Propane, '.-', color='C0')
ax1[0, 0].set_xlabel('residence time [s]')
ax1[0, 0].set_title('Emission Index CO2', color='C0')
ax1[0, 0].set_ylabel('EI [g/kg]')
ax1[0, 1].plot(states_Propane.tres, states_Propane.EI_CO_Propane, '.-', color='C1')
ax1[0, 1].set_xlabel('residence time [s]')
ax1[0, 1].set_title('Emission Index CO', color='C1')
ax1[0, 1].set_ylabel('EI [g/kg]')
ax1[1, 0].plot(states_Propane.tres, states_Propane.EI_NO2_Propane, '.-', color='C2')
ax1[1, 0].set_xlabel('residence time [s]')
ax1[1, 0].set_title('Emission Index NO2', color='C2')
ax1[1, 0].set_ylabel('EI [g/kg]')
ax1[1, 1].plot(states_Propane.tres, states_Propane.EI_NO_Propane, '.-', color='C3')
ax1[1, 1].set_xlabel('residence time [s]')
ax1[1, 1].set_title('Emission Index NO', color='C3')
ax1[1, 1].set_ylabel('EI [g/kg]')

plt.show()