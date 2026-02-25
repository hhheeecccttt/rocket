import numpy as np
import json
import math
import matplotlib.pyplot as plt

EARTH_RADIUS = 6371000
EARTH_MASS = 6e24
SEA_LEVEL_GRAVITY = 9.8
GRAVITATIONAL_CONSTANT = 6.67e-11
SEA_LEVEL_AIR_DENSITY = 1.2251
H = 8500
SPEED_OF_SOUND = 340

T0 = 288.15 
LAPSE = -0.0065  
R = 287             
GAMMA = 1.4

index = 0

with open('stats.json', 'r') as file:
    data = json.load(file)
    rocket = data[index]
    
    name = rocket["name"]
    wetMass = rocket["wetMass"]
    dryMass = rocket["dryMass"]
    burnTime = rocket["burnTime"]
    exhaustVelocity = rocket["exhaustVelocity"]
    dragCoefficient = rocket["dragCoefficient"]
    frontalArea = rocket["frontalArea"]

print(name)
fuelMass = wetMass - dryMass
thrust = 0
massFlowRate = fuelMass / burnTime

rho = 1.225
g = 9.8
launchAngle = 90
theta = np.radians(launchAngle)
mach = 340

dt = 1
t_max = 1000
steps = int(t_max / dt)
final_step = steps - 1

time = np.zeros(steps)
position = np.zeros((steps, 2))
velocity = np.zeros((steps, 2))
acceleration = np.zeros((steps, 2))
mass = np.zeros(steps)
angle = np.zeros(steps)

mass[0] = wetMass
angle[0] = launchAngle
velocity[0] = [0, 0]
position[0] = [0, 0]

def calculateAirDensity(height):
    return SEA_LEVEL_AIR_DENSITY * np.exp(-height / H)

def calculateGravityVector(height):
    h = height[1]
    g = SEA_LEVEL_GRAVITY * (EARTH_RADIUS / (EARTH_RADIUS + h))**2
    return np.array([0, -g])

def calculateDragVector(v, d):
    height = d[1]
    rho = calculateAirDensity(height)
    speed = np.linalg.norm(v)

    if speed == 0:
        return np.array([0.0, 0.0])

    mach = calculateMach(v, height)
    Cd = calculateDragCoefficient(mach)

    drag = 0.5 * rho * Cd * frontalArea * speed**2
    return -drag * (v / speed)

def calculateAcceleration(d, v, m, T):
    gravity = calculateGravityVector(d)
    drag = calculateDragVector(v, d)
    total_force = T + drag + gravity * m
    return total_force / m

def speedOfSound(height):
    h = max(height, 0)

    if h < 11000:
        T = T0 + LAPSE * h
    else:
        T = 216.65

    return np.sqrt(GAMMA * R * T)
    
def calculateMach(v, height):
    speed = np.linalg.norm(v)
    a = speedOfSound(height)
    return speed / a

def calculateDragCoefficient(M):
    Cd0 = dragCoefficient
    
    rise = 0.25 / (1 + np.exp(-20*(M - 0.9)))
    decay = np.exp(-(M - 1.05)/0.7)
    
    if M < 0.9:
        return Cd0
    else:
        return Cd0 + rise * decay

def RK4(d, v, m, T):
    k1v = v
    k1a = calculateAcceleration(d, v, m, T)

    k2v = v + 0.5 * dt * k1a
    k2a = calculateAcceleration(d + 0.5 * dt * k1v, k2v, m, T)

    k3v = v + 0.5 * dt * k2a
    k3a = calculateAcceleration(d + 0.5 * dt * k2v, k3v, m, T)

    k4v = v + dt * k3a
    k4a = calculateAcceleration(d + dt * k3v, k4v, m, T)

    pos = d + (dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)
    vel = v + (dt / 6) * (k1a + 2 * k2a + 2 * k3a + k4a)
    
    return pos, vel

for i in range(steps - 1):
    time[i+1] = time[i] + dt

    if (time[i] < burnTime):
        thrust = massFlowRate * exhaustVelocity
        mass[i + 1] = mass[i] - massFlowRate * dt
    else:
        mass[i + 1] = dryMass
        thrust = 0

    theta -= 0.1 * dt
    angle[i + 1] = theta
    thrust_vector = thrust * np.array([np.cos(theta), np.sin(theta)])

    acceleration[i + 1] = calculateAcceleration(position[i], velocity[i], mass[i],thrust_vector)
    position[i + 1], velocity[i + 1] = RK4(position[i], velocity[i], mass[i], thrust_vector)

    if position[i+1,1] < 0:
        position[i+1,1] = 0
        final_step = i+1
        break


time = time[1:final_step]
velocity = velocity[1:final_step]
position = position[1:final_step]
acceleration = acceleration[1:final_step]
mass = mass[1:final_step]
angle = angle[1:final_step]


mach_values = np.linspace(0, 3, 500)
cd_values = [calculateDragCoefficient(M) for M in mach_values]

plt.figure(figsize=(10,6))
plt.subplot(3,1,1)
plt.plot(position[:,0], position[:,1])
plt.xlabel("Horizontal distance (m)")
plt.ylabel("Altitude (m)")
plt.title(f"Rocket Trajectory at {launchAngle}° Launch")

plt.subplot(3,1,2)
plt.plot(time, velocity[:,1])
plt.xlabel("Time (s)")
plt.ylabel("Vertical velocity (m/s)")
plt.title("Vertical Velocity vs Time")

plt.subplot(3,1,3)
plt.plot(time, -angle * 180 / math.pi + 90)
plt.xlabel("Time (s)")
plt.ylabel("Angle (deg)")
plt.title("Angle vs Time")

'''
theta = np.linspace(0, 2*np.pi, 1000)
earth_x = EARTH_RADIUS * np.cos(theta)
earth_y = EARTH_RADIUS * np.sin(theta)

plt.plot(earth_x, earth_y)
plt.plot(position[:,0], position[:,1])
plt.gca().set_aspect('equal')
'''

plt.tight_layout()
plt.show()