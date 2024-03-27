import matplotlib.pyplot as plt

# Data from the table
# Plotting the bar chart with labeled bars
plt.figure(figsize=(10, 6))
plt.barh(["DER", "ICA", "Climb"], [100, 100, 75], color='red', label='Operation of Fan')
plt.barh(["Ascent", "Cruise", "Descent", "Downleg"], [50, 7, 10, 100], color='skyblue', label='Ram Inlet')
plt.barh(["Reserve Climb"], [100], color='red')
plt.barh(["Reserve Loiter", "Reserve Descent"], [100, 100], color='skyblue')
plt.barh(["Baseleg", "Final Approach"], [100, 100], color='red')

plt.xlabel('Percentage of Design Power (%)')
plt.title('Percentage of Design Power Across Flight Segments')

# Create a legend based on the labels
plt.legend()

plt.gca().invert_yaxis()  
plt.grid(axis='x') 

plt.show()