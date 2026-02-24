import numpy as np
import matplotlib.pyplot as plt

EARTH_RADIUS = 6371000
EARTH_MASS = 6e24
SEA_LEVEL_GRAVITY = 9.8
GRAVITATIONAL_CONSTANT = 6.67e-11
SEA_LEVEL_AIR_DENSITY = 1.2251
H = 8500

CD = 0.125
A = 2.1
EXHAUST_VELOCITY = 2000

wet_mass = 13000
dry_mass = 4000
fuel_mass = wet_mass - dry_mass

thrust = 0
burn_time = 65
mass_flow_rate = fuel_mass / burn_time

rho = 1.225
g = 9.8

dt = 0.01
t_max = 3000
steps = int(t_max / dt)
final_step = steps - 1

time = np.zeros(steps)
velocity = np.zeros(steps)
position = np.zeros(steps)
acceleration = np.zeros(steps)
mass = np.zeros(steps)

mass[0] = wet_mass

def calculateAirDensity(height):
    return SEA_LEVEL_AIR_DENSITY * np.exp(-height / H)

def calculateGravity(height):
    return SEA_LEVEL_GRAVITY * (EARTH_RADIUS / (EARTH_RADIUS + height)) ** 2

def calculateDerivatives(d, v, m, T):
    rho = calculateAirDensity(d)
    g = calculateGravity(d)

    drag = np.sign(v) * 0.5 * rho * CD * A * v ** 2 
    acc = (T - m * g - drag) / m

    return v, acc

def RK4(d, v, m, T):
    k1v, k1a = calculateDerivatives(d, v, m, T)
    k2v, k2a = calculateDerivatives(d + 0.5 * dt * k1v, v + 0.5 * dt * k1a, m, T)
    k3v, k3a = calculateDerivatives(d + 0.5 * dt * k2v, v + 0.5 * dt * k2a, m, T)
    k4v, k4a = calculateDerivatives(d + dt * k3v, v + dt * k3a, m, T)

    pos = d + (dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)
    vel = v + (dt / 6) * (k1a + 2 * k2a + 2 * k3a + k4a)
    
    return pos, vel

for i in range(steps - 1):
    time[i+1] = time[i] + dt

    if (time[i] < burn_time):
        thrust = mass_flow_rate * EXHAUST_VELOCITY
        mass[i + 1] = mass[i] - mass_flow_rate * dt
    else:
        mass[i + 1] = dry_mass
        thrust = 0

    _, acceleration[i+1] = calculateDerivatives(position[i], velocity[i], mass[i],thrust)
    position[i+1], velocity[i+1] = RK4(position[i], velocity[i], mass[i], thrust)

    if position[i+1] < 0:
        final_step = i
        break

time = time[1:final_step]
velocity = velocity[1:final_step]
position = position[1:final_step]
acceleration = acceleration[1:final_step]
mass = mass[1:final_step]

plt.figure(figsize=(10,6))
plt.subplot(4,1,1)
plt.plot(time, position)
plt.ylabel("Altitude (m)")

plt.subplot(4,1,2)
plt.plot(time, velocity)
plt.ylabel("Velocity (m/s)")
plt.xlabel("Time (s)")

plt.subplot(4,1,3)
plt.plot(time, acceleration)
plt.ylabel("Acceleration (m/s2)")
plt.xlabel("Time (s)")

plt.subplot(4,1,4)
plt.plot(time, mass)
plt.ylabel("Mass (kg)")
plt.xlabel("Time (s)")

plt.tight_layout()
plt.show()