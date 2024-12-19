import serial
import pandas as pd
import time
import numpy as np


ser = serial.Serial('/dev/tty.usbmodem11301', 115200, timeout=1) # Change port as per connection

# Wait for the Arduino to reset
time.sleep(2)

# Initialize an empty list to store the data
data = []

try:
    while True:
        # Read a line from the serial port
        line = ser.readline().decode('utf-8').strip()
        if line:
            print(line)  # Print to console for debugging
            data.append(line.split(','))  # Split the line into a list and append to data
except KeyboardInterrupt:
    print("Program terminated.")
finally:
    ser.close()

# Convert the data list to a pandas DataFrame
df = pd.DataFrame(data, columns=['TC_Temp'])

# Save the DataFrame to a CSV file
df.to_csv('/Users/sai/Documents/Arduino/thermocouple_data.csv', index=False)
print("Data saved to thermocouple_data.csv")