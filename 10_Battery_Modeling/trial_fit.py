import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define the data
data = np.array([[0.0, 3.844504115757135],
                 [48.51355893779851, 3.8054220116789006],
                 [178.96703485690614, 3.7913979503294954],
                 [368.62800278292434, 3.7861096128139584],
                 [677.2198900901326, 3.711228629722568],
                 [938.3613580044248, 3.6436059191857604],
                 [1294.403664704549, 3.561466663539786],
                 [1615.0809471322596, 3.4471752538636533],
                 [1947.5622053891793, 3.340962923008372],
                 [2333.4662257764444, 3.2196594826575367],
                 [2570.999742032317, 3.135878614479023],
                 [2898.2122057800398, 2.9187753570507264],
                 [3005.8550846993903, 2.754039539410427],
                 [3060.075201488396, 2.6043948312656835],
                 [3084.3867013750464, 2.5018292253934007]])

x_data = data[:, 0]
y_data = data[:, 1]

# Define the sinh model function
def sinh_model(x, a, b, c):
    return a * np.sinh(b * x) + c

# Fit the model to the data
initial_guess = [1, 0.001, 3]  # Initial guess for a, b, and c
params, covariance = curve_fit(sinh_model, x_data, y_data, p0=initial_guess)

# Extract the parameters
a, b, c = params

# Generate model predictions
x_fit = np.linspace(min(x_data), max(x_data), 200)
y_fit = sinh_model(x_fit, a, b, c)

# Plot the data and the fit
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_fit, y_fit, label=f'Fit: y = {a:.3f} * sinh({b:.6f} * x) + {c:.3f}', color='red')
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.legend()
plt.show()

print(f"Fitted parameters: a = {a:.3f}, b = {b:.6f}, c = {c:.3f}")
