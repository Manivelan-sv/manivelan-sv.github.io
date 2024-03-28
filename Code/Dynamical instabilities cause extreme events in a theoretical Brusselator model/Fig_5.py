'''
Dynamical instabilities cause extreme events in a theoretical Brusselator model.
Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

# Read data from the file "1.11.dat"
data = np.loadtxt('data_file.dat')
t = data[:, 0]
x = data[:, 1]

# Plot x against t in the first subplot
plt.figure(figsize=(9, 5))

# Calculate the analytic signal using the Hilbert transform
yh = hilbert(x)

# Calculate the amplitude and phase of the analytic signal
amplitude = np.abs(yh)

# Calculate the phase values and unwrap them
phase = np.angle(yh)
phase_unwrapped = np.unwrap(phase / np.pi)  # Corrected line

# Plot the real and imaginary parts of the analytic signal in the first subplot
plt.subplot(2, 1, 1)
plt.plot(t, np.imag(yh), 'r', label='Imaginary')
plt.plot(t, np.real(yh), 'g', label='Real')
plt.legend()
plt.ylabel('Hilbert', fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=20)  # Adjust tick label size

# Plot the unwrapped phase values in the second subplot
plt.subplot(2, 1, 2)
plt.plot(t, phase_unwrapped)
plt.xlabel('Time', fontsize=25)
plt.ylabel('Δφ', fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=20)  # Adjust tick label size

# Show the plot
plt.show()

