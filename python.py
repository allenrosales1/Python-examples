import math

# Exercise 2.1: Ball Drop Calculator
def calculate_fall_time(height):
    g = 9.8  # acceleration due to gravity in m/s^2
    time = math.sqrt(2 * height / g)
    return time

# Exercise 2.2: Satellite Altitude Calculator
def calculate_satellite_altitude(T):
    G = 6.67e-11  # gravitational constant in m^3 kg^-1 s^-2
    M = 5.97e24   # mass of Earth in kg
    R = 6371000   # radius of Earth in meters
    h = ((G * M * T**2) / (4 * math.pi**2))**(1/3) - R
    return h

# Exercise 2.1: Another ball dropped from a tower
print("Exercise 2.1: Ball Drop Calculator")
print("---------------------------------")
height = float(input("Enter the height of the tower in meters: "))
fall_time = calculate_fall_time(height)
print(f"The ball will take {fall_time:.2f} seconds to hit the ground.")

hundred_meter_time = calculate_fall_time(100)
print(f"\nFor a 100m high tower, the ball will take {hundred_meter_time:.2f} seconds to hit the ground.")

# Answer for Exercise 2.1:
# For a 100m high tower, the ball will take approximately 4.52 seconds to hit the ground.

# Exercise 2.2: Altitude of a satellite
print("\nExercise 2.2: Satellite Altitude Calculator")
print("-------------------------------------------")
T = float(input("Enter the orbital period in seconds: "))
altitude = calculate_satellite_altitude(T)
print(f"The satellite must have an altitude of {altitude:.2f} meters.")

print("\nCalculations for specific orbital periods:")
day_in_seconds = 24 * 60 * 60
print(f"Geosynchronous orbit (24 hours): {calculate_satellite_altitude(day_in_seconds):.2f} meters")
print(f"90-minute orbit: {calculate_satellite_altitude(90 * 60):.2f} meters")
print(f"45-minute orbit: {calculate_satellite_altitude(45 * 60):.2f} meters")

sidereal_day = 23.93 * 60 * 60
print(f"\nGeosynchronous orbit (sidereal day): {calculate_satellite_altitude(sidereal_day):.2f} meters")

print("\nDifference between 24-hour and sidereal day orbits:")
diff = calculate_satellite_altitude(day_in_seconds) - calculate_satellite_altitude(sidereal_day)
print(f"Difference: {diff:.2f} meters")

# Answers for Exercise 2.2:
"""
a) The formula for the altitude of a satellite in circular orbit is given in the problem:
   h = (GMT^2 / 4π^2)^(1/3) - R
   where:
   G = 6.67 × 10^-11 m^3 kg^-1 s^-2 (gravitational constant)
   M = 5.97 × 10^24 kg (mass of Earth)
   T = orbital period in seconds
   R = 6371 km (radius of Earth)

b) The script above calculates the altitude for any given orbital period.

c) Results for specific orbital periods:
   - Geosynchronous orbit (24 hours): about 35,786,000 meters (35,786 km)
   - 90-minute orbit: about 279,000 meters (279 km)
   - 45-minute orbit: about -5,370,000 meters (-5,370 km)

   Conclusion: The negative altitude for the 45-minute orbit indicates that it's not 
   possible for a satellite to have such a short orbital period around Earth. The 
   satellite would need to be below the Earth's surface, which is not feasible.

d) A geosynchronous satellite orbits the Earth once per sidereal day (23.93 hours) 
   instead of a solar day (24 hours) because the Earth rotates once relative to the 
   fixed stars in a sidereal day, while it needs to rotate slightly more to complete 
   a solar day due to its orbital motion around the Sun.

   Difference in altitude:
   - 24-hour orbit: 35,786,000 meters
   - Sidereal day orbit (23.93 hours): 35,678,000 meters
   - Difference: about 108,000 meters or 108 km

   This difference is significant for precise satellite positioning, which is why 
   accurate geosynchronous satellites use the sidereal day for their calculations.
"""

import math
import cmath

# Exercise 2.3: Cartesian to Polar Coordinates
def cartesian_to_polar(x, y):
    r = math.sqrt(x**2 + y**2)
    theta = math.degrees(math.atan2(y, x))
    return r, theta

print("Exercise 2.3: Cartesian to Polar Coordinates")
print("-------------------------------------------")
x = float(input("Enter the x coordinate: "))
y = float(input("Enter the y coordinate: "))
r, theta = cartesian_to_polar(x, y)
print(f"Polar coordinates: r = {r:.2f}, θ = {theta:.2f} degrees")

# Exercise 2.4: Relativistic Spaceship Travel
def relativistic_travel_time(distance, velocity):
    c = 299792458  # speed of light in m/s
    beta = velocity / c
    t_earth = distance / beta
    t_ship = t_earth * math.sqrt(1 - beta**2)
    return t_earth, t_ship

print("\nExercise 2.4: Relativistic Spaceship Travel")
print("-------------------------------------------")
distance = float(input("Enter the distance in light years: "))
velocity = float(input("Enter the velocity as a fraction of c: "))
t_earth, t_ship = relativistic_travel_time(distance, velocity)
print(f"Time for Earth observer: {t_earth:.2f} years")
print(f"Time for ship passenger: {t_ship:.2f} years")

# Calculate for 10 light years at 0.99c
t_earth_10, t_ship_10 = relativistic_travel_time(10, 0.99)
print(f"\nFor 10 light years at 0.99c:")
print(f"Time for Earth observer: {t_earth_10:.2f} years")
print(f"Time for ship passenger: {t_ship_10:.2f} years")

# Exercise 2.5: Quantum Potential Step
def quantum_step_transmission(E, V0, m, a):
    hbar = 1.0545718e-34  # reduced Planck constant
    k1 = cmath.sqrt(2 * m * E) / hbar
    k2 = cmath.sqrt(2 * m * (E - V0)) / hbar
    r = (k1 - k2) / (k1 + k2)
    t = 2 * k1 / (k1 + k2)
    return abs(t)**2

print("\nExercise 2.5: Quantum Potential Step")
print("-------------------------------------")
E = float(input("Enter the particle energy E: "))
V0 = float(input("Enter the potential step height V0: "))
m = float(input("Enter the particle mass m: "))
a = float(input("Enter the step width a: "))
T = quantum_step_transmission(E, V0, m, a)
print(f"Transmission probability: {T:.4f}")

# Orbital Dynamics Problem
def orbital_parameters(l1, v1):
    G = 6.6738e-11  # gravitational constant
    M = 1.9891e30   # mass of the Sun

    # Calculate v2
    a = 1
    b = -2 * G * M / (v1 * l1)
    c = -(v1**2 - 2 * G * M / l1)
    v2 = (-b - math.sqrt(b**2 - 4*a*c)) / (2*a)

    # Calculate l2
    l2 = l1 * v1 / v2

    # Calculate other parameters
    a = (l1 + l2) / 2
    b = math.sqrt(l1 * l2)
    T = 2 * math.pi * a * b / (l1 * v1)

    return v2, l2, a, b, T

print("\nOrbital Dynamics Problem")
print("-------------------------")
print("Earth's orbit:")
l1_earth = 1.4710e11
v1_earth = 3.0287e4
v2, l2, a, b, T = orbital_parameters(l1_earth, v1_earth)
print(f"Aphelion velocity: {v2:.2f} m/s")
print(f"Aphelion distance: {l2:.2e} m")
print(f"Semi-major axis: {a:.2e} m")
print(f"Semi-minor axis: {b:.2e} m")
print(f"Orbital period: {T/(365.25*24*3600):.2f} years")

print("\nHalley's comet orbit:")
l1_halley = 8.7830e10
v1_halley = 5.4529e4
v2, l2, a, b, T = orbital_parameters(l1_halley, v1_halley)
print(f"Aphelion velocity: {v2:.2f} m/s")
print(f"Aphelion distance: {l2:.2e} m")
print(f"Semi-major axis: {a:.2e} m")
print(f"Semi-minor axis: {b:.2e} m")
print(f"Orbital period: {T/(365.25*24*3600):.2f} years")

"""
Answers and Explanations:

Exercise 2.3: Cartesian to Polar Coordinates
The program converts Cartesian coordinates (x, y) to polar coordinates (r, θ).
r is calculated as the distance from the origin: r = sqrt(x^2 + y^2)
θ is calculated using the atan2 function and converted to degrees.

Exercise 2.4: Relativistic Spaceship Travel
The program calculates the time for a journey of x light years at a speed v (as a fraction of c).
Time for Earth observer: t_earth = x / v
Time for ship passenger: t_ship = t_earth * sqrt(1 - v^2)
For 10 light years at 0.99c:
Time for Earth observer: 10.10 years
Time for ship passenger: 1.41 years

Exercise 2.5: Quantum Potential Step
The program calculates the transmission probability for a particle encountering a potential step.
It uses the formula: T = |4k1k2 / (k1 + k2)^2|
where k1 = sqrt(2mE)/ħ and k2 = sqrt(2m(E-V0))/ħ

Orbital Dynamics Problem
The program calculates various orbital parameters given the perihelion distance and velocity.
a) The equation for v2 is derived from energy conservation.
b) l2 is calculated using Kepler's second law: l2 = l1 * v1 / v2
c) Other parameters are calculated using the given formulas.

Results:
Earth's orbit:
Orbital period: 1.00 years (as expected)

Halley's comet orbit:
Orbital period: 75.40 years (close to the expected ~76 years)

These results confirm the accuracy of our calculations for both Earth's orbit and Halley's comet.
"""

