import cantera           as ct
import numpy             as np
import matplotlib.pyplot as plt

#######################################################################
# Plug Flow Reactor (PFR): modeled as a series of infinitesimal PFR
#######################################################################

# The plug flow reactor is represented by a linear chain of zero-dimensional
# reactors. The gas at the inlet to the first one has the specified inlet
# composition, and for all others the inlet composition is fixed at the
# composition of the reactor immediately upstream. Since in a PFR model there
# is no diffusion, the upstream reactors are not affected by any downstream
# reactors, and therefore the problem may be solved by simply marching from
# the first to last reactor, integrating each one to steady state.

# import kinetic mechanism for Jet A + air
#JetA_PFR = ct.Solution('JetFuelSurrogate.yaml')             # [1]    Simpler kinetic mechanism, faster simulation
JetA_PFR = ct.Solution('chem.yaml')                        # [2]    More accurate kinetic mechanism, slower simulation
                                                            
# set initial conditions                                    
T_0    = 700                                                # [K]    combustor inlet temperature 
P_0    = 20 * ct.one_atm                                    # [atm]  combustor pressure 
length = 0.6                                                # [m]    combustor length 
area   = 0.15                                               # [m**2] combustor cross-sectional area
u_0    = 8                                                  # [m/s]  combustor inlet velocity 
equiv_ratio = 1                                             # [-]    equivalence ratio (0.5 lean combustion, 1 stoichiometric combustion)
n_steps = 100                                               # [-]    number of infinitesimal PFR that compose the PFR

# define the fuel and oxidizer composition by species 
# and molar fraction    
#fuel     = 'N-C12H26:0.6, A1CH3:0.2, A1:0.2'                                                                                                              # [1] Simpler kinetic mechanism, faster simulation    
fuel = 'NC10H22:0.16449, NC12H26:0.34308, NC16H34:0.10335, IC8H18:0.08630, NC7H14:0.07945, C6H5C2H5: 0.07348, C6H5C4H9: 0.05812, C10H7CH3: 0.10972'      # [2] More accurate kinetic mechanism, slower simulation    
oxidizer = 'O2:0.2095, N2:0.7809, AR:0.0093, CO2:0.0003'                                                                                                  # air composition

# compute additional combustor parameters 
mass_flow_rate_PFR = u_0 * JetA_PFR.density * area          # [kg/s] mass flow rate of Jet A + air inside the combustor 
dz = length / n_steps                                       # [m]    infinitesimal length of a single PFR
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
sim_JetA_PFR.rtol = 1e-9                                                                                               # set the relative tolerance for the integrator
sim_JetA_PFR.atol = 1e-21                                                                                              # set the absolute tolerance for the integrator

# Initialize total mass emission dictionary
total_emissions_JetA_PFR = {species: 0.0 for species in JetA_PFR.species_names}                                        # initialize the total emissions for each species

# define vectors and variables
PFR_z = (np.arange(n_steps) + 1) * dz                                                                                  # vector of length for the entire combustor
residence_time_PFR_n = np.zeros_like(PFR_z)                                                                            # initialize the residence time in each infinitesimal PFR
residence_time_PFR = np.zeros_like(PFR_z)                                                                              # initialize the residence time for the entire combustor
X_CO2 = np.zeros_like(PFR_z)                                                                                           # initialize the vector of molar ratio for CO2
X_CO = np.zeros_like(PFR_z)                                                                                            # initialize the vector of molar ratio for CO
X_H2O = np.zeros_like(PFR_z)                                                                                           # initialize the vector of molar ratio for H2O
X_NO2 = np.zeros_like(PFR_z)                                                                                           # initialize the vector of molar ratio for NO2
X_NO = np.zeros_like(PFR_z)                                                                                            # initialize the vector of molar ratio for NO
#X_soot = np.zeros_like(PFR_z)                                                                                          # initialize the vector of molar ratio for soot
Y_CO2 = np.zeros_like(PFR_z)                                                                                           # initialize the vector of mass ratio for CO2
Y_CO = np.zeros_like(PFR_z)                                                                                            # initialize the vector of mass ratio for CO
Y_H2O = np.zeros_like(PFR_z)                                                                                           # initialize the vector of mass ratio for H2O
Y_NO2 = np.zeros_like(PFR_z)                                                                                           # initialize the vector of mass ratio for NO2
Y_NO = np.zeros_like(PFR_z)                                                                                            # initialize the vector of mass ratio for NO
#Y_soot = np.zeros_like(PFR_z)                                                                                          # initialize the vector of mass ratio for soot
M_CO2 = np.zeros_like(PFR_z)                                                                                           # initialize the vector of mass for CO2
M_CO = np.zeros_like(PFR_z)                                                                                            # initialize the vector of mass for CO
M_H2O = np.zeros_like(PFR_z)                                                                                           # initialize the vector of mass for H2O
M_NO2 = np.zeros_like(PFR_z)                                                                                           # initialize the vector of mass ratio for NO2
M_NO = np.zeros_like(PFR_z)                                                                                            # initialize the vector of mass ratio for NO
#M_soot = np.zeros_like(PFR_z)                                                                                          # initialize the vector of mass ratio for soot
EI_CO2_JetA_PFR = np.zeros_like(PFR_z)                                                                                 # initialize the vector of EI for CO2
EI_CO_JetA_PFR = np.zeros_like(PFR_z)                                                                                  # initialize the vector of EI for CO
EI_H2O_JetA_PFR = np.zeros_like(PFR_z)                                                                                 # initialize the vector of EI for H2O
EI_NO2_JetA_PFR = np.zeros_like(PFR_z)                                                                                 # initialize the vector of EI for NO2
EI_NO_JetA_PFR = np.zeros_like(PFR_z)                                                                                  # initialize the vector of EI for NO
#EI_soot_JetA_PFR = np.zeros_like(PFR_z)                                                                                # initialize the vector of EI for soot
mass = 0                                                                                                               # initialize the mass of Jet A + air
mass_air = 0                                                                                                           # initialize the mass air
states_JetA_PFR = ct.SolutionArray(combustor_JetA_PFR.thermo, extra=['mass', 'mass_air'])                              # define the variables of the solution

# iterate through the PFR cells
for n in range(n_steps):
    # Set the state of the reservoir to match that of the previous reactor
    JetA_PFR.TDY = combustor_JetA_PFR.thermo.TDY                                                                       # set the same temperature, density and mass ratio
    inlet_JetA_PFR.syncState()                                                                                         # synchronize the internal state
    # integrate the reactor forward in time until steady state is reached
    sim_JetA_PFR.reinitialize()                                                                                        # reinitialize the state 
    try:
        sim_JetA_PFR.advance_to_steady_state()                                                                         # advance the simulation until it reaches a steady-state condition
    except ct.CanteraError as e:
        print(f"Error at step {n}: {e}")
        break
    
    residence_time_PFR_n[n] = combustor_JetA_PFR.mass / mass_flow_rate_PFR                                             # residence time in the infinitesimal PFR
    residence_time_PFR[n] = np.sum(residence_time_PFR_n[:n+1])                                                         # residence time for the entire combustor
    
    # Extract mole fractions
    X_CO2[n] = combustor_JetA_PFR.thermo['CO2'].X[0]                                                                   # vector of molar ratio for CO2
    X_CO[n] = combustor_JetA_PFR.thermo['CO'].X[0]                                                                     # vector of molar ratio for CO
    X_H2O[n] = combustor_JetA_PFR.thermo['H2O'].X[0]                                                                   # vector of molar ratio for H2O
    X_NO2[n] = combustor_JetA_PFR.thermo['NO2'].X[0]                                                                   # vector of molar ratio for NO2
    X_NO[n] = combustor_JetA_PFR.thermo['NO'].X[0]                                                                     # vector of molar ratio for NO
    #X_soot[n] = combustor_JetA_PFR.thermo['CSOLID'].X[0]                                                               # vector of molar ratio for soot
    
    # Extract mass fractions
    Y_CO2[n] = combustor_JetA_PFR.thermo['CO2'].Y[0]                                                                   # vector of mass ratio for CO2
    Y_CO[n] = combustor_JetA_PFR.thermo['CO'].Y[0]                                                                     # vector of mass ratio for CO
    Y_H2O[n] = combustor_JetA_PFR.thermo['H2O'].Y[0]                                                                   # vector of mass ratio for H2O
    Y_NO2[n] = combustor_JetA_PFR.thermo['NO2'].Y[0]                                                                   # vector of molar ratio for NO2
    Y_NO[n] = combustor_JetA_PFR.thermo['NO'].Y[0]                                                                     # vector of molar ratio for NO
    #Y_soot[n] = combustor_JetA_PFR.thermo['CSOLID'].Y[0]                                                               # vector of molar ratio for soot  
    
    if n == 0:
        M_CO2[n] = Y_CO2[n] * combustor_JetA_PFR.mass                                                                  
        M_CO[n] = Y_CO[n] * combustor_JetA_PFR.mass                                                                    
        M_H2O[n] = Y_H2O[n] * combustor_JetA_PFR.mass                                                                  
        M_NO2[n] = Y_NO2[n] * combustor_JetA_PFR.mass                                                                  
        M_NO[n] = Y_NO[n] * combustor_JetA_PFR.mass                                                                    
        #M_soot[n] = Y_soot[n] * combustor_JetA_PFR.mass                                                                
    else:
        M_CO2[n] = M_CO2[n-1] + Y_CO2[n] * combustor_JetA_PFR.mass                                                     # vector of mass for CO2
        M_CO[n] = M_CO[n-1] + Y_CO[n] * combustor_JetA_PFR.mass                                                        # vector of mass for CO
        M_H2O[n] = M_H2O[n-1] + Y_H2O[n] * combustor_JetA_PFR.mass                                                     # vector of mass for H2O
        M_NO2[n] = M_NO2[n-1] + Y_NO2[n] * combustor_JetA_PFR.mass                                                     # vector of mass for NO2
        M_NO[n] = M_NO[n-1] + Y_NO[n] * combustor_JetA_PFR.mass                                                        # vector of mass for NO
        #M_soot[n] = M_soot[n-1] + Y_soot[n] * combustor_JetA_PFR.mass                                                  # vector of mass for soot       
    
    for species in JetA_PFR.species_names:
        total_emissions_JetA_PFR[species] += combustor_JetA_PFR.thermo[species].Y * combustor_JetA_PFR.mass            # total emissions for each species
    
    mass_air += (equiv_ratio*14.5*combustor_JetA_PFR.mass)/(1 + equiv_ratio*14.5)                                      # mass of Jet A + air
    
    mass += combustor_JetA_PFR.mass                                                                                    # mass of air
    
    EI_CO2_JetA_PFR[n] = total_emissions_JetA_PFR['CO2'] / (mass - mass_air)                                           # vector of EI for CO2
    EI_CO_JetA_PFR[n] = total_emissions_JetA_PFR['CO'] / (mass - mass_air)                                             # vector of EI for CO
    EI_H2O_JetA_PFR[n] = total_emissions_JetA_PFR['H2O'] / (mass - mass_air)                                           # vector of EI for H2O
    EI_NO2_JetA_PFR[n] = total_emissions_JetA_PFR['NO2'] / (mass - mass_air)                                           # vector of EI for NO2
    EI_NO_JetA_PFR[n] = total_emissions_JetA_PFR['NO'] / (mass - mass_air)                                             # vector of EI for NO
    #EI_soot_JetA_PFR[n] = total_emissions_JetA_PFR['CSOLID'] / (mass - mass_air)                                       # vector of EI for soot    
    
    states_JetA_PFR.append(combustor_JetA_PFR.thermo.state, mass=mass, mass_air=mass_air)                              # write output data

# Plot results
f, ax1 = plt.subplots(3, 1, figsize=(16, 12))
f.suptitle('JetA Surrogate PFR Emissions')
subtitle = f'Equivalence ratio: {equiv_ratio}, Temperature: {JetA_PFR.T:.1f} K, Pressure: {JetA_PFR.P/ct.one_atm:.1f} atm, Length: {length} m, Area: {area} m^2'
plt.figtext(0.5, 0.925, subtitle, ha='center', fontsize=12)
ax1[0].plot(PFR_z, EI_CO2_JetA_PFR, '.-', color='C0')
ax1[0].axhline(y=3.16, color='r', linestyle='--')
ax1[0].annotate('Typical EI value: 3.16', xy=(0.5, 3.16), xytext=(0.5, 3.14), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
ax1[0].set_title('Emission Index CO2', color='C0')
ax1[0].set_ylabel('EI [kg/kg]')
ax1[1].plot(PFR_z, EI_CO_JetA_PFR, '.-', color='C1')
ax1[1].axhline(y=0.05, color='r', linestyle='--')                                                                                                                 # https://ntrs.nasa.gov/api/citations/19750007129/downloads/19750007129.pdf, https://pdf.sciencedirectassets.com/271798/1-s2.0-S1352231015X00148/1-s2.0-S1352231015301722/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEFQaCXVzLWVhc3QtMSJHMEUCIQCbJLgpzlmLNnfVH4WG%2FuFB%2FqOu8xrYTmZDE%2FFdwY9CRQIgehzVJdCWAPcetTxMSNhAjOt304L1xdH5qkaf5ZtBqVQqswUITRAFGgwwNTkwMDM1NDY4NjUiDFPyX8EpoZkVNqV3lSqQBWYvEJUcWvUdCATtDtCfyxKBy4jWz2XQVw7nsF9S4bmfJqHSb4XTUr5b1SQzBkqZusEjar1Qi65i0CxHA8FSi2T4hdDbv399n2x9e3fLe%2FaI3MBt7ln%2BzuRZTJ7WH18TSE7v1I5rsejXhUloB3jFtkAn%2FWVF2o7uJnrpGPJgsIxMhYMhDk6OZb4FM%2BnhS%2Fb1FE2Q275EvoF2Dvce%2FlJ5srqnH38xoC2CoxxEWmQEDf9dZcQdaANorm2HkZCb2Wfe2o3tFnJruPAgxBOWMrSrt1B01CP3nfahP0R4svM03hcbnG0AFPkZs4ASaPE6cuJgIM0atyIF6QOndgrfRz5Itq5hSixBdjDtwy8AkDPbeDIT%2FTChYcRvJug46fRys56a9gB8w2IIKR5PbGreEDqD76O%2FSFm2%2FthB%2BZPEBGx1VGi7BPtHx%2FCCMxGPQToqx7M0XdDGWB3CLP7ae81xPu%2BtQo5B%2BnUPSVUzCExpO%2Fy6yISXRVeXOGFdbKUUAXFvn6EG0HPnRsPGke7B59RGvcycjRNT77rUD71wsBrrQ4sTmPpjLPrQsgg03z9a8QqWRNuUx4nkGyBCTjaxCCyeH81ZDZcTO1wpOya6AGbrjh61%2Fztmo2D8LiCkhR%2F9%2BuYkSusLq4JhFs%2BPh5%2BUtPA%2BwzZQGmKNCByIsJrvL5dEU8NeSEkwMp2Lc9AcHz8LZFi4FfXwUpnE9N%2FIdEmlUvjPVtmvdi8RICk6sLcOxMmHdSWUvS6aH1qegi1b5p2aul%2BieilyIVYcw533DSA5e9QxF6CXX%2B9lnwFRn2ijs1FLdrE1cBMxLM3hUqmUpcQfAYoX%2BHIxzsiCJjcFzEZmUAvQPO0RCdy1tbv9sElfG%2F1j35eQ3%2BeNMO%2Fj2bUGOrEBlJLvANktdqzp1uCuSH0eCO%2BAH7r5OuiXzQkXuVQRIgNZh7BOpRQUe17BOSmJ8zLSPk2%2BeXQboT8VMd4A6Mau1gAAEd9KCR4LnvOtxC0LbJzih3f0M%2FMaiygaPuRBP5uQ2BcDX4Sl%2Bq7YMu%2F0cLhr4%2BogeEEXJCHeqQy2aU8sxyJhgoc01%2FK4x7DrFG71HUQ6RfC5YTNFBn%2FDvs67d3KmUiuyhnNYzZTh7aGE8pzGonoe&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240809T210310Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY67WYP4PM%2F20240809%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=49f95b8789decc69cd2507fffa9b06d7f74a8dbcb638af573e4e905e76b5923a&hash=77bdba4d533e0c7c4359d88b85c90a90fb1a7e6fc1ef817544435a0d2e984525&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1352231015301722&tid=spdf-0c1e36cc-69f9-44bf-b0e4-bae0ff444eb9&sid=21a47c3f56535040da8be6c1225a593131dfgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=17155c0557070d5703&rr=8b0ab0fc59411cde&cc=us
ax1[1].annotate('Typical EI value: 0.05', xy=(0.5, 0.05), xytext=(0.5, 0.03), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->')) 
ax1[1].set_title('Emission Index CO', color='C1')
ax1[1].set_ylabel('EI [kg/kg]')
ax1[2].plot(PFR_z, EI_H2O_JetA_PFR, '.-', color='C2')
ax1[2].set_xlabel('Combustor length [m]')
ax1[2].axhline(y=1.34, color='r', linestyle='--')
ax1[2].annotate('Typical EI value: 1.34', xy=(0.5, 1.34), xytext=(0.5, 1.32), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
ax1[2].set_title('Emission Index H2O', color='C2')
ax1[2].set_ylabel('EI [kg/kg]')

#f, ax1 = plt.subplots(3, 1, figsize=(16, 12))
#f.suptitle('JetA Surrogate PFR Emissions')
#subtitle = f'Equivalence ratio: {equiv_ratio}, Temperature: {JetA_PFR.T:.1f} K, Pressure: {JetA_PFR.P/ct.one_atm:.1f} atm, Length: {length} m, Area: {area} m^2'
#plt.figtext(0.5, 0.925, subtitle, ha='center', fontsize=12)
#ax1[0].plot(PFR_z, EI_NO2_JetA_PFR, '.-', color='C0')
#ax1[0].axhline(y=0.01, color='r', linestyle='--')
#ax1[0].annotate('Typical EI value: 0.01', xy=(0.5, 0.01), xytext=(0.5, 0.008), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
#ax1[0].set_title('Emission Index NO2', color='C0')
#ax1[0].set_ylabel('EI [kg/kg]')
#ax1[1].plot(PFR_z, EI_NO_JetA_PFR, '.-', color='C1')
#ax1[1].axhline(y=0.01, color='r', linestyle='--')
#ax1[1].annotate('Typical EI value: 0.01', xy=(0.5, 0.01), xytext=(0.5, 0.008), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
#ax1[1].set_title('Emission Index NO', color='C1')
#ax1[1].set_ylabel('EI [kg/kg]')
#ax1[2].plot(PFR_z, EI_soot_JetA_PFR, '.-', color='C2')
#ax1[2].axhline(y=0.00004, color='r', linestyle='--')
#ax1[2].annotate('Typical EI value: 0.00004', xy=(0.5, 0.00004), xytext=(0.5, -0.00196), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
#ax1[2].set_xlabel('Combustor length [m]')
#ax1[2].set_title('Emission Index C-soot', color='C2')
#ax1[2].set_ylabel('EI [kg/kg]')

plt.figure(figsize=(16, 12))
plt.plot(PFR_z, X_CO2, '.-', label='CO2 Mole Fraction', color='C0')
plt.plot(PFR_z, X_CO, '.-', label='CO Mole Fraction', color='C1')
plt.plot(PFR_z, X_H2O, '.-', label='H2O Mole Fraction', color='C2')
#plt.plot(PFR_z, X_NO2, '.-', label='NO2 Mole Fraction', color='C3')
#plt.plot(PFR_z, X_NO, '.-', label='NO Mole Fraction', color='C4')
#plt.plot(PFR_z, X_soot, '.-', label='Soot Mole Fraction', color='C5')
plt.xlabel('Combustor length [m]')
plt.ylabel('Mole Fraction')
plt.title('Mole Fraction of CO2, CO, H2O, NO2, NO and soot vs. combustor length')
plt.legend()
plt.grid(True)

plt.figure(figsize=(16, 12))
plt.plot(PFR_z, Y_CO2, '.-', label='CO2 Mass Fraction', color='C0')
plt.plot(PFR_z, Y_CO, '.-', label='CO Mass Fraction', color='C1')
plt.plot(PFR_z, Y_H2O, '.-', label='H2O Mass Fraction', color='C2')
#plt.plot(PFR_z, Y_NO2, '.-', label='NO2 Mass Fraction', color='C3')
#plt.plot(PFR_z, Y_NO, '.-', label='NO Mass Fraction', color='C4')
#plt.plot(PFR_z, Y_soot, '.-', label='Soot Mass Fraction', color='C5')
plt.xlabel('Combustor length [m]')
plt.ylabel('Mass Fraction')
plt.title('Mass Fraction of CO2, CO, H2O, NO2, NO and soot vs. combustor length')
plt.legend()
plt.grid(True)

plt.figure(figsize=(16, 12))
plt.plot(PFR_z, M_CO2, '.-', label='CO2 Mass', color='C0')
plt.plot(PFR_z, M_CO, '.-', label='CO Mass', color='C1')
plt.plot(PFR_z, M_H2O, '.-', label='H2O Mass', color='C2')
#plt.plot(PFR_z, M_NO2, '.-', label='NO2 Mass', color='C3')
#plt.plot(PFR_z, M_NO, '.-', label='NO Mass', color='C4')
#plt.plot(PFR_z, M_soot, '.-', label='Soot Mass', color='C5')
plt.xlabel('Combustor length [m]')
plt.ylabel('Mass[kg]')
plt.title('Mass of CO2, CO, H2O, NO2, NO and soot vs. combustor length')
plt.legend()
plt.grid(True)

plt.show()