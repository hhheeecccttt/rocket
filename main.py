import numpy as np
import matplotlib.pyplot as plt

EARTH_RADIUS = 6371000
EARTH_MASS = 6e24
TM_DIFFERENTIAL = 15
SEA_LEVEL_GRAVITY = 9.8
GRAVITATIONAL_CONSTANT = 6.67e-11
g = 9.8
SEA_LEVEL_AIR_DENSITY = 1.2251
H = 8500
CD = 0.75
A = 80

wet_mass = 13000
dry_mass = 4000

rho = 1.225

max_thrust = 250000
thrust = 0
burn_time = 65
burn_rate = (wet_mass - dry_mass) / burn_time

dt = 0.01
t_max = 1500
steps = int(t_max / dt)
final_step = steps - 1

time = np.zeros(steps)
velocity = np.zeros(steps)
position = np.zeros(steps)
mass = np.zeros(steps)

mass[0] = wet_mass

def calculateAirDensity(height, i):
    return SEA_LEVEL_AIR_DENSITY * np.exp(-position[i] / H)

def calculateGravity(height, i):
    return SEA_LEVEL_GRAVITY * (EARTH_RADIUS / (EARTH_RADIUS + height)) ** 2

for i in range(steps - 2):
    time[i+1] = time[i] + dt

    if (time[i] < burn_time):
        mass[i] = wet_mass - (wet_mass - dry_mass) * (time[i] / burn_time)
        thrust = max_thrust
    else:
        mass[i] = dry_mass
        thrust = 0

    rho = calculateAirDensity(position[i], i)
    g = calculateGravity(position[i], i)

    drag = 0.5 * rho * CD * A * velocity[i] ** 2 * np.sign(velocity[i])

    acceleration = (thrust - mass[i] * g - drag) / mass[i]

    velocity[i+1] = velocity[i] + acceleration * dt
    position[i+1] = position[i] + velocity[i] * dt + 0.5 * acceleration * dt**2

    if position[i+1] < 0:
        final_step = i
        break
    
time = time[:final_step]
velocity = velocity[:final_step]
position = position[:final_step]
mass = mass[:final_step]


plt.figure(figsize=(10,6))
plt.subplot(2,1,1)
plt.plot(time, position)
plt.ylabel("Altitude (m)")

plt.subplot(2,1,2)
plt.plot(time, velocity)
plt.ylabel("Velocity (m/s)")
plt.xlabel("Time (s)")

plt.tight_layout()
plt.show()