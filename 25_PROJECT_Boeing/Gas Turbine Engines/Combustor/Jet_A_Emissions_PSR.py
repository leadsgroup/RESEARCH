import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

#######################################################################
# Perfectly Stirred Reactor (PSR)
#######################################################################
#Reference: Yoo, Kwang-Hee & Kim, Jong-Chan & Sung, Hong-Gye & Zhang, Liwei & Yang, Vigor. (2011). Flow Dynamics in Combustors with Multi-Element Swirl Injectors.

JetA_PSR = ct.Solution('JetFuelSurrogate.yaml')
#JetA_PSR = ct.Solution('chem.yaml') 

# Create a Reservoir for the inlet, set to a fuel/air mixture at a specified equivalence ratio
equiv_ratio = 0.5  # lean combustion
JetA_PSR.TP = 2000.0, 25*ct.one_atm
length = 0.6 
area = 0.05 

fuel = 'N-C12H26:0.6, A1CH3:0.2, A1:0.2'
#for new fuel chemical mechanism, n-tetradecane(C14H30)is missing. 
#fuel = 'NC10H22:0.16449, NC12H26:0.34308, NC16H34:0.10335, IC8H18:0.08630, NC7H14:0.07945, C6H5C2H5: 0.07348, C6H5C4H9: 0.05812, C10H7CH3: 0.10972'
oxidizer = 'O2:0.2095, N2:0.7809, AR:0.0093, CO2:0.0003'  

JetA_PSR.set_equivalence_ratio(equiv_ratio, fuel, oxidizer)
inlet_JetA_PSR = ct.Reservoir(JetA_PSR)

# Create the combustor, and fill it initially with a mixture consisting of the equilibrium products of the inlet mixture.
JetA_PSR.equilibrate('HP')
combustor_JetA_PSR = ct.IdealGasReactor(JetA_PSR)
combustor_JetA_PSR.volume = length*area

# Create a reservoir for the exhaust
exhaust_JetA_PSR = ct.Reservoir(JetA_PSR)

# Use a variable mass flow rate to keep the residence time in the reactor constant.
def mdot(t):
    return combustor_JetA_PSR.mass / residence_time_PSR

mass_air = (combustor_JetA_PSR.thermo['O2'].Y + combustor_JetA_PSR.thermo['N2'].Y + combustor_JetA_PSR.thermo['AR'].Y)*combustor_JetA_PSR.mass

inlet_mfc_JetA_PSR = ct.MassFlowController(inlet_JetA_PSR, combustor_JetA_PSR, mdot=mdot)

# A PressureController has a baseline mass flow rate matching the 'primary' MassFlowController.
outlet_mfc_JetA_PSR = ct.PressureController(combustor_JetA_PSR, exhaust_JetA_PSR, primary=inlet_mfc_JetA_PSR, K=0.01)

# The simulation only contains one reactor
sim_JetA_PSR = ct.ReactorNet([combustor_JetA_PSR])

# Initialize total mass emission dictionary
total_emissions_JetA_PSR = {species: 0.0 for species in JetA_PSR.species_names}

# Run a loop over decreasing residence times, until the reactor is extinguished.
states_JetA_PSR = ct.SolutionArray(JetA_PSR, extra=['tres', 'EI_CO2_JetA_PSR', 'EI_CO_JetA_PSR', 'EI_H2O_JetA_PSR', 'X_CO2', 'X_CO', 'X_H2O', 'Y_CO2', 'Y_CO', 'Y_H2O','M_CO2', 'M_CO', 'M_H2O'])
#states_JetA_PSR = ct.SolutionArray(JetA_PSR, extra=['tres', 'EI_CO2_JetA_PSR', 'EI_CO_JetA_PSR', 'EI_H2O_JetA_PSR', 'EI_NO2_JetA_PSR', 'EI_NO_JetA_PSR', 'EI_soot_JetA_PSR'])
#residence_time_PFR = combustor_JetA_PFR.mass / mass_flow_rate_PFR
residence_time_PSR = 0.01  # starting residence time

while combustor_JetA_PSR.T > 2100:
    sim_JetA_PSR.initial_time = 0.0  # reset the integrator
    sim_JetA_PSR.advance_to_steady_state()

    # Ensure combustor.T is a float
    temperature_JetA_PSR = combustor_JetA_PSR.T if np.isscalar(combustor_JetA_PSR.T) else combustor_JetA_PSR.T.item()

    ## Compute mass flow rates
    for species in JetA_PSR.species_names:
        species_index_JetA_PSR = JetA_PSR.species_index(species)
        mdot_species_JetA_PSR = combustor_JetA_PSR.thermo[species].Y * combustor_JetA_PSR.mass / residence_time_PSR
        total_emissions_JetA_PSR[species] += mdot_species_JetA_PSR * residence_time_PSR
    
    # Compute emission index 
    EI_CO2_JetA_PSR = total_emissions_JetA_PSR['CO2'] / (combustor_JetA_PSR.mass - mass_air)
    EI_CO_JetA_PSR = total_emissions_JetA_PSR['CO'] / (combustor_JetA_PSR.mass - mass_air)      
    EI_H2O_JetA_PSR = total_emissions_JetA_PSR['H2O'] / (combustor_JetA_PSR.mass - mass_air)      
    #EI_NO2_JetA_PSR = total_emissions_JetA_PSR['NO2'] / (combustor_JetA_PSR.mass - mass_air)      
    #EI_NO_JetA_PSR = total_emissions_JetA_PSR['NO'] / (combustor_JetA_PSR.mass - mass_air)      
    #EI_soot_JetA_PSR = total_emissions_JetA_PSR['CSOLID'] / (combustor_JetA_PSR.mass - mass_air)
    
    # Extract mole fractions
    X_CO2 = combustor_JetA_PSR.thermo['CO2'].X[0]
    X_CO = combustor_JetA_PSR.thermo['CO'].X[0]
    X_H2O = combustor_JetA_PSR.thermo['H2O'].X[0]  
    
    # Extract mass fractions
    Y_CO2 = combustor_JetA_PSR.thermo['CO2'].Y[0]
    Y_CO = combustor_JetA_PSR.thermo['CO'].Y[0]
    Y_H2O = combustor_JetA_PSR.thermo['H2O'].Y[0]     
    
    M_CO2 = total_emissions_JetA_PSR['CO2'][0]
    M_CO = total_emissions_JetA_PSR['CO'][0]
    M_H2O = total_emissions_JetA_PSR['H2O'][0]
    
    states_JetA_PSR.append(combustor_JetA_PSR.thermo.state, tres=residence_time_PSR, EI_CO2_JetA_PSR=EI_CO2_JetA_PSR, EI_CO_JetA_PSR=EI_CO_JetA_PSR, EI_H2O_JetA_PSR=EI_H2O_JetA_PSR, X_CO2=X_CO2, X_CO=X_CO, X_H2O=X_H2O, Y_CO2=Y_CO2, Y_CO=Y_CO, Y_H2O=Y_H2O, M_CO2=M_CO2, M_CO=M_CO, M_H2O=M_H2O)
    #states_JetA_PSR.append(combustor_JetA_PSR.thermo.state, tres=residence_time_PSR, EI_CO2_JetA_PSR=EI_CO2_JetA_PSR, EI_CO_JetA_PSR=EI_CO_JetA_PSR, EI_H2O_JetA_PSR=EI_H2O_JetA_PSR, EI_NO2_JetA_PSR=EI_NO2_JetA_PSR, EI_NO_JetA_PSR=EI_NO_JetA_PSR, EI_soot_JetA_PSR=EI_soot_JetA_PSR)
    residence_time_PSR *= 0.9 

# Print total emissions for each species
print("\nJetA PSR emissions for each species (in kg):")
for species, total_emission in total_emissions_JetA_PSR.items():
    total_emission_scalar = total_emission if np.isscalar(total_emission) else total_emission.item()
    print(f"{species}: {total_emission_scalar:.6e} kg")

# Plot results
f, ax1 = plt.subplots(3, 1, figsize=(16, 12))
f.suptitle('JetA Surrogate PSR Emissions')
subtitle = f'Equivalence ratio: {equiv_ratio}, Temperature: {JetA_PSR.T:.1f} K, Pressure: {JetA_PSR.P/ct.one_atm:.1f} atm, Length: {length} m, Area: {area} m^2'
plt.figtext(0.5, 0.925, subtitle, ha='center', fontsize=12)
ax1[0].plot(states_JetA_PSR.tres, states_JetA_PSR.EI_CO2_JetA_PSR, '.-', color='C0')
ax1[0].axhline(y=3.16, color='r', linestyle='--')
ax1[0].annotate('Typical EI value: 3.16', xy=(0.5, 3.16), xytext=(0.5, 3.18), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
ax1[0].set_title('Emission Index CO2', color='C0')
ax1[0].set_ylabel('EI [kg/kg]')
ax1[1].plot(states_JetA_PSR.tres, states_JetA_PSR.EI_CO_JetA_PSR, '.-', color='C1')
ax1[1].set_title('Emission Index CO', color='C1')
ax1[1].set_ylabel('EI [kg/kg]')
ax1[2].plot(states_JetA_PSR.tres, states_JetA_PSR.EI_H2O_JetA_PSR, '.-', color='C2')
ax1[2].set_xlabel('residence time [s]')
ax1[2].axhline(y=1.34, color='r', linestyle='--')
ax1[2].annotate('Typical EI value: 1.34', xy=(0.5, 1.34), xytext=(0.5, 1.36), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
ax1[2].set_title('Emission Index H2O', color='C2')
ax1[2].set_ylabel('EI [kg/kg]')



#f, ax1 = plt.subplots(3, 1, figsize=(16, 12))
#f.suptitle('JetA Surrogate PSR Emissions')
#subtitle = f'Equivalence ratio: {equiv_ratio}, Temperature: {JetA_PSR.T:.1f} K, Pressure: {JetA_PSR.P/ct.one_atm:.1f} atm, Length: {length} m, Area: {area} m^2'
#plt.figtext(0.5, 0.925, subtitle, ha='center', fontsize=12)
#ax1[0].plot(states_JetA_PSR.tres, states_JetA_PSR.EI_NO2_JetA_PSR, '.-', color='C0')
##ax1[0].axhline(y=3.16, color='r', linestyle='--')
##ax1[0].annotate('Typical EI value: 3.16', xy=(0.5, 3.16), xytext=(0.5, 3.18), textcoords='data', color='r', arrowprops=dict(facecolor='r', arrowstyle='->'))
#ax1[0].set_title('Emission Index NO2', color='C0')
#ax1[0].set_ylabel('EI [kg/kg]')
#ax1[1].plot(states_JetA_PSR.tres, states_JetA_PSR.EI_NO_JetA_PSR, '.-', color='C1')
#ax1[1].set_title('Emission Index NO', color='C1')
#ax1[1].set_ylabel('EI [kg/kg]')
#ax1[2].plot(states_JetA_PSR.tres, states_JetA_PSR.EI_soot_JetA_PSR, '.-', color='C2')
#ax1[2].set_xlabel('residence time [s]')
#ax1[2].set_title('Emission Index C-soot', color='C2')
#ax1[2].set_ylabel('EI [kg/kg]')

plt.figure(figsize=(16, 12))
plt.plot(states_JetA_PSR.tres, states_JetA_PSR.X_CO2, '.-', label='CO2 Mole Fraction', color='C0')
plt.plot(states_JetA_PSR.tres, states_JetA_PSR.X_CO, '.-', label='CO Mole Fraction', color='C1')
plt.plot(states_JetA_PSR.tres, states_JetA_PSR.X_H2O, '.-', label='H2O Mole Fraction', color='C2')
plt.xlabel('Residence Time [s]')
plt.ylabel('Mole Fraction')
plt.title('Mole Fraction of CO2, CO, and H2O vs. Residence Time')
plt.legend()
plt.grid(True)

plt.figure(figsize=(16, 12))
plt.plot(states_JetA_PSR.tres, states_JetA_PSR.Y_CO2, '.-', label='CO2 Mass Fraction', color='C0')
plt.plot(states_JetA_PSR.tres, states_JetA_PSR.Y_CO, '.-', label='CO Mass Fraction', color='C1')
plt.plot(states_JetA_PSR.tres, states_JetA_PSR.Y_H2O, '.-', label='H2O Mass Fraction', color='C2')
plt.xlabel('Residence Time [s]')
plt.ylabel('Mass Fraction')
plt.title('Mass Fraction of CO2, CO, and H2O vs. Residence Time')
plt.legend()
plt.grid(True)

plt.figure(figsize=(16, 12))
plt.plot(states_JetA_PSR.tres, M_CO2, '.-', label='CO2 Mass', color='C0')
plt.plot(states_JetA_PSR.tres, M_CO, '.-', label='CO Mass', color='C1')
plt.plot(states_JetA_PSR.tres, M_H2O, '.-', label='H2O Mass', color='C2')
plt.xlabel('Residence Time [s]')
plt.ylabel('Mass[kg]')
plt.title('Mass of CO2, CO, and H2O vs. Residence Time')
plt.legend()
plt.grid(True)

plt.show()