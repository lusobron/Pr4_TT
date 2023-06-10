import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.plotting import plot

# Define the initial and final positions
r1 = np.array([20.0e6, 20.0e6, 0])  # [m]
r2 = np.array([-20.0e6, 10.0e6, 0]) # [m]

# Define the time of flight
tof = 1.0 * 86400  # [s]

# Compute the Lambert solution
ss_i = Orbit.from_vectors(Earth, r1, [0, 0, 0] * u.m / u.s, None)
ss_f = Orbit.from_vectors(Earth, r2, [0, 0, 0] * u.m / u.s, None)
man_lambert = Maneuver.lambert(ss_i, ss_f, tof)

# Extract the initial and final velocities from the maneuver
V1, V2 = man_lambert.get_total_cost()[:2]

# Print the results
print(f"Initial velocity: {V1}")
print(f"Final velocity: {V2}")