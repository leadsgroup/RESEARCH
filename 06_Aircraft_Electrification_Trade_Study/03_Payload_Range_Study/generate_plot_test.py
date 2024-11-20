import os
import json
import matplotlib.pyplot as plt

# Define the directory path
directory = r"C:\Research\RESEARCH\06_Aircraft_Electrification_Trade_Study\03_Payload_Range_Study\A220_data"

# Check if the directory exists
if not os.path.exists(directory):
    raise FileNotFoundError(f"Directory '{directory}' does not exist. Please check the path.")

# List all JSON files
json_files = [f for f in os.listdir(directory) if f.endswith('.json')]
if not json_files:
    raise FileNotFoundError("No JSON files found in the specified directory.")

# Initialize Set 1 color scheme
cmap = plt.colormaps["Set1"]  # Get the Set1 colormap
colors = [cmap(i / 8) for i in range(9)]  # Generate exactly 9 colors (0 to 8)

# Explicitly map fuel names to specific colors
fuel_colors = {
    "Jet A1": colors[0],   # Red
    "Butanol": colors[1],  # Blue
    "Ethane": colors[2],   # Green
    "Ethanol": colors[3],  # Purple
    "Liquid Natural Gas": colors[4],      # Orange
    "Liquid Petroleum Gas": colors[5],      # Yellow
    "Methane": colors[6],  # Brown
    "Propane": colors[7],  # Pink
    "Propanol": colors[8]  # Gray
}

# Prepare the plot
plt.figure(figsize=(10, 6))  # Set DPI for on-screen display
plt.xlabel("Range (nmi)")
plt.ylabel("Payload (kg)")

# Loop through each JSON file and plot the data
for file in json_files:
    file_path = os.path.join(directory, file)
    
    # Open and parse the JSON file
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    # Assume the JSON contains 'range' and 'payload' keys
    range_data = data.get("range", [])
    payload_data = data.get("payload", [])
    
    if not range_data or not payload_data:
        print(f"Skipping {file}: Missing 'range' or 'payload' data.")
        continue
    
    # Convert range data to nautical miles
    range_data = [r * 0.000539957 for r in range_data]
    
    # Extract label and color
    label = file.replace("_range.json", "")
    color = fuel_colors.get(label, "black")  # Default to black if fuel not explicitly mapped
    
    # Plot the data
    plt.plot(range_data, payload_data, label=label, color=color)

# Add legend and grid
plt.legend()
plt.grid(True)

# Save the plot with high DPI
plt.savefig("A220_payload_range_plot.png", dpi=600)  # Save with 300 DPI

# Show the plot
plt.show()
