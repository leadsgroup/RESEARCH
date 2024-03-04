# Control System for mass flow rate optimization 

class CoolingSystem:
    
    
    def __init__(self, max_flow_rate, min_flow_rate, max_gradient, timestep):
        self.max_flow_rate = max_flow_rate
        self.min_flow_rate = min_flow_rate
        self.max_gradient = max_gradient
        self.timestep = timestep
        self.temperature_history = []  # Store temperature history over last few time steps

    def update_temperature_history(self, temperature):
        self.temperature_history.append(temperature)
        # If the length of temperature history exceeds the desired number of time steps,
        # remove the oldest temperature data
        if len(self.temperature_history) > self.timestep:
            self.temperature_history.pop(0)

    def calculate_flow_rate(self):
        if len(self.temperature_history) < self.timestep:
            
            return self.min_flow_rate  # Return minimum flow rate if there's not enough data

        # Calculate the gradient of temperature over the last few time steps
        temperature_gradient = (self.temperature_history[-1] - self.temperature_history[0]) / self.timestep

        # Scale the gradient to fit within the range of the maximum gradient
        scaled_gradient = min(max(temperature_gradient, -self.max_gradient), self.max_gradient)

        # Calculate the mass flow rate based on the scaled gradient
        mass_flow_rate = (scaled_gradient / self.max_gradient) * (self.max_flow_rate - self.min_flow_rate) + self.min_flow_rate

        return mass_flow_rate


# Example usage:
max_flow_rate = 1  # Maximum allowable mass flow rate
min_flow_rate = 0.5   # Minimum allowable mass flow rate
max_gradient = 5     # Maximum allowed gradient of temperature change
timestep = 2         # Number of time steps to consider for gradient calculation

cooling_system = CoolingSystem(max_flow_rate, min_flow_rate, max_gradient, timestep)

# Simulated temperature data over time
temperature_data = [30, 32, 35, 38, 80, 42, 41, 39, 37, 35]

for temperature in temperature_data:
    cooling_system.update_temperature_history(temperature)
    mass_flow_rate = cooling_system.calculate_flow_rate()
    print(f"Temperature: {temperature}°C, Mass Flow Rate: {mass_flow_rate} kg/s")
