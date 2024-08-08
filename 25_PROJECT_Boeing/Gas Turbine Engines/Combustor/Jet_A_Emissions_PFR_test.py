import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

#######################################################################
# Plug Flow Reactor (PFR)
#######################################################################

JetA_PFR = ct.Solution('JetFuelSurrogate.yaml')

T_0 = 2000.0  # inlet temperature [K]
P_0 = 25 * ct.one_atm  # constant pressure [Pa]
length = 0.6  # *approximate* PFR length [m]
area = 0.05  # cross-sectional area [m**2]
equiv_ratio = 0.5  # lean combustion
JetA_PFR.TP = T_0, P_0
fuel = 'N-C12H26:0.6, A1CH3:0.2, A1:0.2'
oxidizer = 'O2:0.2095, N2:0.7809, AR:0.0093, CO2:0.0003'
JetA_PFR.set_equivalence_ratio(equiv_ratio, fuel, oxidizer)
u_0 = 15  # inflow velocity [m/s]
mass_flow_rate_PFR = u_0 * JetA_PFR.density * area
n_steps = 100

dz = length / n_steps
PFR_vol = area * dz

# The plug flow reactor is represented by a linear chain of zero-dimensional
# reactors. The gas at the inlet to the first one has the specified inlet
# composition, and for all others the inlet composition is fixed at the
# composition of the reactor immediately upstream. Since in a PFR model there
# is no diffusion, the upstream reactors are not affected by any downstream
# reactors, and therefore the problem may be solved by simply marching from
# the first to last reactor, integrating each one to steady state.

# create a new reactor
combustor_JetA_PFR = ct.IdealGasReactor(JetA_PFR)
combustor_JetA_PFR.volume = PFR_vol

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
inlet_JetA_PFR = ct.Reservoir(JetA_PFR, name='inlet')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
exhaust_JetA_PFR = ct.Reservoir(JetA_PFR, name='exhaust')

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
inlet_mfc_JetA_PFR = ct.MassFlowController(inlet_JetA_PFR, combustor_JetA_PFR, mdot=mass_flow_rate_PFR)

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
outlet_mfc_JetA_PFR = ct.PressureController(combustor_JetA_PFR, exhaust_JetA_PFR, primary=inlet_mfc_JetA_PFR, K=1e-5)

sim_JetA_PFR = ct.ReactorNet([combustor_JetA_PFR])

# Set the relative and absolute tolerances for the integrator
sim_JetA_PFR.rtol = 1e-9
sim_JetA_PFR.atol = 1e-21

# Initialize total mass emission dictionary
total_emissions_JetA_PFR = {species: 0.0 for species in JetA_PFR.species_names}

# define time, space, and other information vectors
PFR_z = (np.arange(n_steps) + 1) * dz
residence_time_PFR_n = np.zeros_like(PFR_z)
residence_time_PFR = np.zeros_like(PFR_z)
u_PFR = np.zeros_like(PFR_z)
X_CO2 = np.zeros_like(PFR_z)
X_CO = np.zeros_like(PFR_z)
X_H2O = np.zeros_like(PFR_z)
Y_CO2 = np.zeros_like(PFR_z)
Y_CO = np.zeros_like(PFR_z)
Y_H2O = np.zeros_like(PFR_z)
M_CO2 = np.zeros_like(PFR_z)
M_CO = np.zeros_like(PFR_z)
M_H2O = np.zeros_like(PFR_z)
EI_CO2_JetA_PFR = np.zeros_like(PFR_z)
EI_CO_JetA_PFR = np.zeros_like(PFR_z)
EI_H2O_JetA_PFR = np.zeros_like(PFR_z)
mass = 0
mass_air = 0
states_JetA_PFR = ct.SolutionArray(JetA_PFR, extra=['mass', 'mass_air'])

# iterate through the PFR cells
for n in range(n_steps):
    # Set the state of the reservoir to match that of the previous reactor
    JetA_PFR.TDY = combustor_JetA_PFR.thermo.TDY
    inlet_JetA_PFR.syncState()
    # integrate the reactor forward in time until steady state is reached
    sim_JetA_PFR.reinitialize()
    try:
        sim_JetA_PFR.advance_to_steady_state()
    except ct.CanteraError as e:
        print(f"Error at step {n}: {e}")
        break
    
    # compute velocity and transform into time
    u_PFR[n] = mass_flow_rate_PFR / area / combustor_JetA_PFR.thermo.density
    residence_time_PFR_n[n] = combustor_JetA_PFR.mass / mass_flow_rate_PFR  # residence time in this reactor
    residence_time_PFR[n] = np.sum(residence_time_PFR_n[:n+1])
    
    # Extract mole fractions
    X_CO2[n] = combustor_JetA_PFR.thermo['CO2'].X[0]
    X_CO[n] = combustor_JetA_PFR.thermo['CO'].X[0]
    X_H2O[n] = combustor_JetA_PFR.thermo['H2O'].X[0]
    
    # Extract mass fractions
    Y_CO2[n] = combustor_JetA_PFR.thermo['CO2'].Y[0]
    Y_CO[n] = combustor_JetA_PFR.thermo['CO'].Y[0]
    Y_H2O[n] = combustor_JetA_PFR.thermo['H2O'].Y[0]
    
    M_CO2[n] = combustor_JetA_PFR.thermo['CO2'].Y[0] * combustor_JetA_PFR.mass
    M_CO[n] = combustor_JetA_PFR.thermo['CO'].Y[0] * combustor_JetA_PFR.mass
    M_H2O[n] = combustor_JetA_PFR.thermo['H2O'].Y[0] * combustor_JetA_PFR.mass
    
    for species in JetA_PFR.species_names:
        total_emissions_JetA_PFR[species] += combustor_JetA_PFR.thermo[species].Y * combustor_JetA_PFR.mass
    
    mass_air += (combustor_JetA_PFR.thermo['O2'].Y + combustor_JetA_PFR.thermo['N2'].Y + combustor_JetA_PFR.thermo['AR'].Y) * combustor_JetA_PFR.mass
    
    mass += combustor_JetA_PFR.mass
    
    EI_CO2_JetA_PFR[n] = total_emissions_JetA_PFR['CO2'] / (mass - mass_air)
    EI_CO_JetA_PFR[n] = total_emissions_JetA_PFR['CO'] / (mass - mass_air)
    EI_H2O_JetA_PFR[n] = total_emissions_JetA_PFR['H2O'] / (mass - mass_air)
    
    # write output data
    states_JetA_PFR.append(combustor_JetA_PFR.thermo.state, mass=mass, mass_air=mass_air)

    
# Plot results
f, ax1 = plt.subplots(3, 1, figsize=(16, 12))
f.suptitle('JetA Surrogate PFR Emissions')
subtitle = f'Equivalence ratio: {equiv_ratio}, Temperature: {JetA_PFR.T:.1f} K, Pressure: {JetA_PFR.P/ct.one_atm:.1f} atm, Length: {length} m, Area: {area} m^2'
plt.figtext(0.5, 0.925, subtitle, ha='center', fontsize=12)
ax1[0].plot(residence_time_PFR, EI_CO2_JetA_PFR, '.-', color='C0')
ax1[0].axhline(y=3.16, color='r', linestyle='--')
ax1[0].annotate('Typical EI value: 3.16', xy=(0.5, 3.16), xytext=(0.5, 3.18), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
ax1[0].set_title('Emission Index CO2', color='C0')
ax1[0].set_ylabel('EI [kg/kg]')
ax1[1].plot(residence_time_PFR, EI_CO_JetA_PFR, '.-', color='C1')
ax1[1].set_title('Emission Index CO', color='C1')
ax1[1].set_ylabel('EI [kg/kg]')
ax1[2].plot(residence_time_PFR, EI_H2O_JetA_PFR, '.-', color='C2')
ax1[2].set_xlabel('residence time [s]')
ax1[2].axhline(y=1.34, color='r', linestyle='--')
ax1[2].annotate('Typical EI value: 1.34', xy=(0.5, 1.34), xytext=(0.5, 1.36), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
ax1[2].set_title('Emission Index H2O', color='C2')
ax1[2].set_ylabel('EI [kg/kg]')

plt.figure(figsize=(16, 12))
plt.plot(residence_time_PFR, X_CO2, '.-', label='CO2 Mole Fraction', color='C0')
plt.plot(residence_time_PFR, X_CO, '.-', label='CO Mole Fraction', color='C1')
plt.plot(residence_time_PFR, X_H2O, '.-', label='H2O Mole Fraction', color='C2')
plt.xlabel('Residence Time [s]')
plt.ylabel('Mole Fraction')
plt.title('Mole Fraction of CO2, CO, H2O, NO2, NO and soot vs. Residence Time')
plt.legend()
plt.grid(True)

plt.figure(figsize=(16, 12))
plt.plot(residence_time_PFR, Y_CO2, '.-', label='CO2 Mass Fraction', color='C0')
plt.plot(residence_time_PFR, Y_CO, '.-', label='CO Mass Fraction', color='C1')
plt.plot(residence_time_PFR, Y_H2O, '.-', label='H2O Mass Fraction', color='C2')
plt.xlabel('Residence Time [s]')
plt.ylabel('Mass Fraction')
plt.title('Mass Fraction of CO2, CO, H2O, NO2, NO and soot vs. Residence Time')
plt.legend()
plt.grid(True)

plt.figure(figsize=(16, 12))
plt.plot(residence_time_PFR, M_CO2, '.-', label='CO2 Mass', color='C0')
plt.plot(residence_time_PFR, M_CO, '.-', label='CO Mass', color='C1')
plt.plot(residence_time_PFR, M_H2O, '.-', label='H2O Mass', color='C2')
plt.xlabel('Residence Time [s]')
plt.ylabel('Mass[kg]')
plt.title('Mass of CO2, CO, H2O, NO2, NO and soot vs. Residence Time')
plt.legend()
plt.grid(True)

plt.show()