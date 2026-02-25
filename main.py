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
launch_angle = 80
theta = np.radians(launch_angle)

dt = 0.1
t_max = 3000
steps = int(t_max / dt)
final_step = steps - 1

time = np.zeros(steps)
position = np.zeros((steps, 2))
velocity = np.zeros((steps, 2))
acceleration = np.zeros((steps, 2))
mass = np.zeros(steps)

mass[0] = wet_mass
velocity[0] = [0, 0]
position[0] = [0, 0]

def calculateAirDensity(height):
    return SEA_LEVEL_AIR_DENSITY * np.exp(-height / H)

def calculateGravityVector(height):
    h = height[1]
    g = SEA_LEVEL_GRAVITY * (EARTH_RADIUS / (EARTH_RADIUS + h))**2
    return np.array([0, -g])

def calculateDragVector(v, d):
    rho = calculateAirDensity(d[1])
    speed = np.sqrt(v[0]*v[0] + v[1]*v[1])
    if speed == 0:
        return np.array([0.0, 0.0])
    drag = 0.5 * rho * CD * A * speed**2
    return -drag * (v / speed)

def calculateAcceleration(d, v, m, T):
    gravity = calculateGravityVector(d)
    drag = calculateDragVector(v, d)
    total_force = T + drag + gravity * m
    return total_force / m

def RK4(d, v, m, T):
    k1v = v
    k1a = calculateAcceleration(d, v, m, T)

    k2v = v + 0.5*dt*k1a
    k2a = calculateAcceleration(d + 0.5 * dt * k1v, k2v, m, T)

    k3v = v + 0.5*dt*k2a
    k3a = calculateAcceleration(d + 0.5 * dt * k2v, k3v, m, T)

    k4v = v + dt * k3a
    k4a = calculateAcceleration(d + dt * k3v, k4v, m, T)

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

    thrust_vector = thrust * np.array([np.cos(theta), np.sin(theta)])

    acceleration[i+1] = calculateAcceleration(position[i], velocity[i], mass[i],thrust_vector)
    position[i+1], velocity[i+1] = RK4(position[i], velocity[i], mass[i], thrust_vector)

    if position[i+1,1] < 0:
        position[i+1,1] = 0
        final_step = i+1
        break


time = time[1:final_step]
velocity = velocity[1:final_step]
position = position[1:final_step]
acceleration = acceleration[1:final_step]
mass = mass[1:final_step]

plt.figure(figsize=(10,6))
plt.subplot(2,1,1)
plt.plot(position[:,0], position[:,1])
plt.xlabel("Horizontal distance (m)")
plt.ylabel("Altitude (m)")
plt.title(f"Rocket Trajectory at {launch_angle}° Launch")

plt.subplot(2,1,2)
plt.plot(time, velocity[:,1])
plt.xlabel("Time (s)")
plt.ylabel("Vertical velocity (m/s)")
plt.title("Vertical Velocity vs Time")

plt.tight_layout()
plt.show()