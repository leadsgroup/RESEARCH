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

# Use reaction mechanism GRI-Mech 3.0. For 0-D simulations, no transport model is necessary.
gas = ct.Solution('gri30.yaml')

# Create a Reservoir for the inlet, set to a methane/air mixture at a specified equivalence ratio
equiv_ratio = 0.5  # lean combustion
gas.TP = 300.0, ct.one_atm
gas.set_equivalence_ratio(equiv_ratio, 'CH4:1.0', 'O2:1.0, N2:3.76')
inlet = ct.Reservoir(gas)

# Create the combustor, and fill it initially with a mixture consisting of the equilibrium products of the inlet mixture.
gas.equilibrate('HP')
combustor = ct.IdealGasReactor(gas)
combustor.volume = 1.0

# Create a reservoir for the exhaust
exhaust = ct.Reservoir(gas)

# Use a variable mass flow rate to keep the residence time in the reactor constant.
def mdot(t):
    return combustor.mass / residence_time

inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)

# A PressureController has a baseline mass flow rate matching the 'primary' MassFlowController.
outlet_mfc = ct.PressureController(combustor, exhaust, primary=inlet_mfc, K=0.01)

# The simulation only contains one reactor
sim = ct.ReactorNet([combustor])

# Initialize total mass emission dictionary
total_emissions = {species: 0.0 for species in gas.species_names}

# Run a loop over decreasing residence times, until the reactor is extinguished.
states = ct.SolutionArray(gas, extra=['tres', 'EI_CO2', 'EI_CO', 'EI_NO2', 'EI_NO'])
residence_time = 0.1  # starting residence time
while combustor.T > 500:
    sim.initial_time = 0.0  # reset the integrator
    sim.advance_to_steady_state()

    # Ensure combustor.T is a float
    temperature = combustor.T if np.isscalar(combustor.T) else combustor.T.item()
    print('tres = {:.2e}; T = {:.1f}'.format(residence_time, temperature))

    # Compute mass flow rates
    mdot_fuel = mdot(0)# * gas.molecular_weights[gas.species_index('CH4')] / gas.density
    for species in gas.species_names:
        species_index = gas.species_index(species)
        mdot_species = combustor.thermo[species].Y * combustor.thermo.density * combustor.volume / residence_time
        total_emissions[species] += mdot_species * residence_time

    # Compute emission index for CO2
    mdot_CO2 = total_emissions['CO2']
    EI_CO2 = mdot_CO2 / mdot_fuel if mdot_fuel > 0 else 0
    mdot_CO = total_emissions['CO']
    EI_CO = mdot_CO / mdot_fuel if mdot_fuel > 0 else 0    
    mdot_NO2 = total_emissions['NO2']
    EI_NO2 = mdot_NO2 / mdot_fuel if mdot_fuel > 0 else 0 
    mdot_NO = total_emissions['NO']
    EI_NO = mdot_NO / mdot_fuel if mdot_fuel > 0 else 0      
    states.append(combustor.thermo.state, tres=residence_time, EI_CO2=EI_CO2, EI_CO=EI_CO, EI_NO2=EI_NO2, EI_NO=EI_NO)
    residence_time *= 0.9  # decrease the residence time for the next iteration

    

# Print total emissions for each species
print("\nTotal emissions for each species (in kg):")
for species, total_emission in total_emissions.items():
    # Ensure total_emission is a float
    total_emission_scalar = total_emission if np.isscalar(total_emission) else total_emission.item()
    print(f"{species}: {total_emission_scalar:.6e} kg")

# Plot results
f, ax1 = plt.subplots(1, 1)
ax1.plot(states.tres, states.EI_CO2, '.-', color='C0')
ax1.set_xlabel('residence time [s]')
ax1.set_ylabel('Emission Index CO2 [kg/kg fuel]', color='C0')
f.tight_layout()

f, ax1 = plt.subplots(1, 1)
ax1.plot(states.tres, states.EI_CO, '.-', color='C0')
ax1.set_xlabel('residence time [s]')
ax1.set_ylabel('Emission Index CO [kg/kg fuel]', color='C0')
f.tight_layout()

f, ax1 = plt.subplots(1, 1)
ax1.plot(states.tres, states.EI_NO2, '.-', color='C0')
ax1.set_xlabel('residence time [s]')
ax1.set_ylabel('Emission Index NO2 [kg/kg fuel]', color='C0')
f.tight_layout()

f, ax1 = plt.subplots(1, 1)
ax1.plot(states.tres, states.EI_NO, '.-', color='C0')
ax1.set_xlabel('residence time [s]')
ax1.set_ylabel('Emission Index NO [kg/kg fuel]', color='C0')
f.tight_layout()

plt.show()
