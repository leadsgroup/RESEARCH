import cantera           as ct
import numpy             as np
import matplotlib.pyplot as plt

#######################################################################
# Plug Flow Reactor (PFR): modeled as a series of infinitesimal PSR
#######################################################################

# The plug flow reactor is represented by a linear chain of zero-dimensional
# reactors. The gas at the inlet to the first one has the specified inlet
# composition, and for all others the inlet composition is fixed at the
# composition of the reactor immediately upstream. Since in a PFR model there
# is no diffusion, the upstream reactors are not affected by any downstream
# reactors, and therefore the problem may be solved by simply marching from
# the first to last reactor, integrating each one to steady state.

# import kinetic mechanism for Jet A + air
JetA_PFR = ct.Solution('JetFuelSurrogate.yaml')             # [1]    Simpler kinetic mechanism, faster simulation
#JetA_PFR = ct.Solution('chem.yaml')                        # [2]    More accurate kinetic mechanism, slower simulation
                                                            
# set initial conditions                                    
T_0    = 700                                                # [K]    combustor inlet temperature 
P_0    = 20 * ct.one_atm                                    # [atm]  combustor pressure 
length = 0.6                                                # [m]    combustor length 
area   = 0.15                                               # [m**2] combustor cross-sectional area
u_0    = 8                                                  # [m/s]  combustor inlet velocity 
equiv_ratio = 1                                             # [-]    equivalence ratio (0.5 lean combustion, 1 stoichiometric combustion)
n_steps = 100                                               # [-]    number of infinitesimal PSR that compose the PFR

# define the fuel and oxidizer composition by species 
# and molar fraction    
fuel     = 'N-C12H26:0.6, A1CH3:0.2, A1:0.2'                                                                                                              # [1] Simpler kinetic mechanism, faster simulation    
#fuel = 'NC10H22:0.16449, NC12H26:0.34308, NC16H34:0.10335, IC8H18:0.08630, NC7H14:0.07945, C6H5C2H5: 0.07348, C6H5C4H9: 0.05812, C10H7CH3: 0.10972'      # [2] More accurate kinetic mechanism, slower simulation    
oxidizer = 'O2:0.2095, N2:0.7809, AR:0.0093, CO2:0.0003'                                                                                                  # air composition

# compute additional combustor parameters 
mass_flow_rate_PFR = u_0 * JetA_PFR.density * area          # [kg/s] mass flow rate of Jet A + air inside the combustor 
dz = length / n_steps                                       # [m]    infinitesimal length of a single PSR
PFR_vol = area * dz                                         # [m**3] combustor volume

# set Jet A + air mixture properties
JetA_PFR.TP = T_0, P_0                                      # set temperature and pressure 
JetA_PFR.set_equivalence_ratio(equiv_ratio, fuel, oxidizer) # set equivalence ratio, fuel and oxidizer

# create the PFR reactor
combustor_JetA_PFR = ct.IdealGasReactor(JetA_PFR)           # create the ideal gas reactor 
combustor_JetA_PFR.volume = PFR_vol                         # set the combustor volume

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
inlet_JetA_PFR = ct.Reservoir(JetA_PFR, name='inlet')       # set the inlet resevoir

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
exhaust_JetA_PFR = ct.Reservoir(JetA_PFR, name='exhaust')   # set the exhaust resevoir

# The mass flow rate into the reactor is fixed by using a
# MassFlowController object.
inlet_mfc_JetA_PFR = ct.MassFlowController(inlet_JetA_PFR, combustor_JetA_PFR, mdot=mass_flow_rate_PFR)                # set the combustor inlet mass flow rate

# Set the outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
outlet_mfc_JetA_PFR = ct.PressureController(combustor_JetA_PFR, exhaust_JetA_PFR, primary=inlet_mfc_JetA_PFR, K=1e-5)  # set the outlet to the downstream resevoir

sim_JetA_PFR = ct.ReactorNet([combustor_JetA_PFR])                                                                     # define the simulation reactor network composed by the combustor

# Set the relative and absolute tolerances for the integrator
sim_JetA_PFR.rtol = 1e-9 
sim_JetA_PFR.atol = 1e-21

# Initialize total mass emission dictionary
total_emissions_JetA_PFR = {species: 0.0 for species in JetA_PFR.species_names}                                        # initialize the total si

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
states_JetA_PFR = ct.SolutionArray(combustor_JetA_PFR.thermo, extra=['mass', 'mass_air'])

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
    
    M_CO2[n] += combustor_JetA_PFR.thermo['CO2'].Y[0] * combustor_JetA_PFR.mass
    M_CO[n] += combustor_JetA_PFR.thermo['CO'].Y[0] * combustor_JetA_PFR.mass
    M_H2O[n] += combustor_JetA_PFR.thermo['H2O'].Y[0] * combustor_JetA_PFR.mass
    
    for species in JetA_PFR.species_names:
        total_emissions_JetA_PFR[species] += combustor_JetA_PFR.thermo[species].Y * combustor_JetA_PFR.mass
    
    mass_air += (equiv_ratio*14.5*combustor_JetA_PFR.mass)/(1 + equiv_ratio*14.5)
    
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
ax1[0].plot(PFR_z, EI_CO2_JetA_PFR, '.-', color='C0')
ax1[0].axhline(y=3.16, color='r', linestyle='--')
ax1[0].annotate('Typical EI value: 3.16', xy=(0.5, 3.16), xytext=(0.5, 3.18), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
ax1[0].set_title('Emission Index CO2', color='C0')
ax1[0].set_ylabel('EI [kg/kg]')
ax1[1].plot(PFR_z, EI_CO_JetA_PFR, '.-', color='C1')
ax1[1].set_title('Emission Index CO', color='C1')
ax1[1].set_ylabel('EI [kg/kg]')
ax1[2].plot(PFR_z, EI_H2O_JetA_PFR, '.-', color='C2')
ax1[2].set_xlabel('Combustor length [m]')
ax1[2].axhline(y=1.34, color='r', linestyle='--')
ax1[2].annotate('Typical EI value: 1.34', xy=(0.5, 1.34), xytext=(0.5, 1.36), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
ax1[2].set_title('Emission Index H2O', color='C2')
ax1[2].set_ylabel('EI [kg/kg]')

plt.figure(figsize=(16, 12))
plt.plot(PFR_z, X_CO2, '.-', label='CO2 Mole Fraction', color='C0')
plt.plot(PFR_z, X_CO, '.-', label='CO Mole Fraction', color='C1')
plt.plot(PFR_z, X_H2O, '.-', label='H2O Mole Fraction', color='C2')
plt.xlabel('Combustor length [m]')
plt.ylabel('Mole Fraction')
plt.title('Mole Fraction of CO2, CO, H2O, NO2, NO and soot vs. Residence Time')
plt.legend()
plt.grid(True)

plt.figure(figsize=(16, 12))
plt.plot(PFR_z, Y_CO2, '.-', label='CO2 Mass Fraction', color='C0')
plt.plot(PFR_z, Y_CO, '.-', label='CO Mass Fraction', color='C1')
plt.plot(PFR_z, Y_H2O, '.-', label='H2O Mass Fraction', color='C2')
plt.xlabel('Combustor length [m]')
plt.ylabel('Mass Fraction')
plt.title('Mass Fraction of CO2, CO, H2O, NO2, NO and soot vs. Residence Time')
plt.legend()
plt.grid(True)

plt.figure(figsize=(16, 12))
plt.plot(PFR_z, M_CO2, '.-', label='CO2 Mass', color='C0')
plt.plot(PFR_z, M_CO, '.-', label='CO Mass', color='C1')
plt.plot(PFR_z, M_H2O, '.-', label='H2O Mass', color='C2')
plt.xlabel('Combustor length [m]')
plt.ylabel('Mass[kg]')
plt.title('Mass of CO2, CO, H2O, NO2, NO and soot vs. Residence Time')
plt.legend()
plt.grid(True)

plt.show()