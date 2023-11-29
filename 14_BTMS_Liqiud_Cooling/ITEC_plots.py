import matplotlib.pyplot as plt
import numpy as np

# Given data for the first set
mass_flow_rate = [0.1, 0.5, 1, 1.5, 2]
channel_width = [0.001068809, 0.001, 0.001, 0.001, 0.002990134]
angle_of_contact_radians = [0.823273618, 0.829030395, 0.378364017, 0.27953316, 0.772574583]
angle_of_contact_degrees = np.degrees(angle_of_contact_radians)

# Plot Channel Width and Angle vs Mass Flow Rate
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot Channel Width on the left Y-axis
ax1.set_xlabel('Mass Flow Rate of Coolant')
ax1.set_ylabel('Channel Width', color='tab:blue')
ax1.plot(mass_flow_rate, channel_width, color='tab:blue', marker='o')
ax1.tick_params(axis='y', labelcolor='tab:blue')

# Create a secondary axis for Angle of Contact on the right Y-axis
ax2 = ax1.twinx()
ax2.set_ylabel('Angle of Contact (degrees)', color='tab:red')
ax2.plot(mass_flow_rate, angle_of_contact_degrees, color='tab:red', linestyle='--', marker='s')
ax2.tick_params(axis='y', labelcolor='tab:red')

# Display the plot
plt.title('Channel Width and Angle of Contact vs Mass Flow Rate')
plt.show()
# Given data for the second set
mass = [42.64409364, 42.75128108, 20.25481331, 15.32135184, 44.998556]
power = [0.013889811, 0.320564344, 0.568661325, 0.951951201, 1.645415373]

# Plot Power and Mass vs Mass Flow Rate
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot Mass on the left Y-axis
ax1.set_xlabel('Mass Flow Rate of Coolant')
ax1.set_ylabel('Mass (kg)', color='tab:green')
ax1.plot(mass_flow_rate, mass, color='tab:green', marker='^')
ax1.tick_params(axis='y', labelcolor='tab:green')

# Create a secondary axis for Power on the right Y-axis
ax2 = ax1.twinx()
ax2.set_ylabel('Power (W)', color='tab:orange')
ax2.plot(mass_flow_rate, power, color='tab:orange', linestyle='--', marker='*')
ax2.tick_params(axis='y', labelcolor='tab:orange')

# Display the plot
plt.title('Mass and Power vs Mass Flow Rate')
plt.show()

# Given data
mass_flow_rate = [0.11, 0.5, 1, 1.5, 2, 2.5]
height_of_stack = [0.010367271, 0.026737309, 0.040017113, 0.048559426, 0.053602384, 0.053677281]
length_of_hot_side = [1.86420135, 1.647889512, 1.713028869, 1.924355652, 2.376950397, 3.452468908]
length_of_cold_side = [1.378999049, 0.490059006, 0.335372408, 0.298018125, 0.30980726, 0.394728004]

# Create a figure and axis
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot mass flow rate and height_of_stack on the left Y-axis
ax1.set_xlabel('Mass Flow Rate')
ax1.set_ylabel('Height of Stack', color='tab:blue')
ax1.plot(mass_flow_rate, height_of_stack, color='tab:blue', marker='o', label='Height of Stack')
ax1.tick_params(axis='y', labelcolor='tab:blue')

# Create a secondary axis for length_of_hot_side and length_of_cold_side on the right Y-axis
ax2 = ax1.twinx()
ax2.set_ylabel('Length', color='tab:red')
ax2.plot(mass_flow_rate, length_of_hot_side, color='tab:red', linestyle='--', marker='s', label='Length of Hot Side')
ax2.plot(mass_flow_rate, length_of_cold_side, color='tab:green', linestyle='-', marker='^', label='Length of Cold Side')
ax2.tick_params(axis='y', labelcolor='tab:red')

# Display the legend
fig.tight_layout()
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# Display the plot
plt.title('Parameters vs Mass Flow Rate')
plt.show()

mass_flow_rate_air = [0.11, 0.5, 1, 1.5, 2, 2.5]
mass = [61.84426288, 50.10399404, 53.34773164, 64.6219048, 91.59547184, 169.7447745]
power_watt = [47651.40283, 12827.62044, 10297.28773, 11069.56922, 12703.43629, 14682.51306]

# Convert power to kilowatts
power_kw = [p / 1000 for p in power_watt]

# Create a figure and axis
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot mass flow rate of air and mass on the left Y-axis
ax1.set_xlabel('Mass Flow Rate of Air')
ax1.set_ylabel('Mass (kg)', color='tab:blue')
ax1.plot(mass_flow_rate_air, mass, color='tab:blue', marker='o', label='Mass')
ax1.tick_params(axis='y', labelcolor='tab:blue')

# Create a secondary axis for power_kw on the right Y-axis
ax2 = ax1.twinx()
ax2.set_ylabel('Power (kW)', color='tab:red')
ax2.plot(mass_flow_rate_air, power_kw, color='tab:red', linestyle='--', marker='s', label='Power (kW)')
ax2.tick_params(axis='y', labelcolor='tab:red')

# Display the legend
fig.tight_layout()
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# Display the plot
plt.title('Mass and Power vs Mass Flow Rate of Air')
plt.show()