import json
import math
import atmosphere

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

EARTH_RADIUS = 6371000
EARTH_MASS = 6e24
EARTH_ANGULAR_VELOCITY = 7.2921159e-5
SEA_LEVEL_GRAVITY = 9.8
GRAVITATIONAL_CONSTANT = 6.67e-11
SEA_LEVEL_AIR_DENSITY = 1.2251
H = 8500
SPEED_OF_SOUND = 340

T0 = 288.15
LAPSE = -0.0065
R = 287
GAMMA = 1.4

index = 1

with open("stats.json", "r") as file:
    data = json.load(file)
    rocket = data[index]

    stages = rocket["stages"]
    initialStage = stages[0]

    name = rocket["name"]
    wetMass = sum(stage["wetMass"] for stage in stages)

    dragCoefficient = rocket["dragCoefficient"]
    frontalArea = rocket["frontalArea"]

    VERTICAL_CLIMB_TIME = rocket["VERTICAL_CLIMB_TIME"]
    PITCH_KICK_DEGREES = rocket["PITCH_KICK_DEGREES"]
    PITCH_KICK_DURATION = rocket["PITCH_KICK_DURATION"]

currentStageIndex = 0
currentStage = stages[currentStageIndex]
stageTime = 0
time = 0

fuelMass = currentStage["wetMass"] - currentStage["dryMass"]
massFlowRate = fuelMass / currentStage["burnTime"]
thrust = 0

print(name)

rho = 1.225
g = 9.8
launchAngle = 90
theta = np.radians(launchAngle)
mach = 340

dt = 1
t_max = 10000
steps = int(t_max / dt)
final_step = steps - 1

time = np.zeros(steps)
mass = np.zeros(steps)
angle = np.zeros(steps)
velocity = np.zeros((steps, 2))
position = np.zeros((steps, 2))
altitude = np.zeros(steps)

mass[0] = wetMass
angle[0] = launchAngle
velocity[0] = [EARTH_ANGULAR_VELOCITY * EARTH_RADIUS, 0]
position[0] = [0, EARTH_RADIUS]

def calculateThrust(seaISP, vacuumISP, height):
    pressure_ratio = atmosphere.density(height) / 1.225
    isp = seaISP + (vacuumISP - seaISP) * (1 - pressure_ratio)
    return massFlowRate * isp * SEA_LEVEL_GRAVITY

def calculateAltitude(d):
    return np.linalg.norm(d) - EARTH_RADIUS

def calculateAirDensity(height):
    return SEA_LEVEL_AIR_DENSITY * np.exp(-height / H)

def calculateGravityVector(d):
    r = np.linalg.norm(d)
    g_mag = GRAVITATIONAL_CONSTANT * EARTH_MASS / r**2
    return -g_mag * (d / r)

def calculateDragVector(v, d):
    height = calculateAltitude(d)

    if height > 80000:
        return np.array([0.0, 0.0])

    relativeSpeed = calculateRelativeVelocity(position[i], velocity[i])
    speed = np.linalg.norm(relativeSpeed)
    if speed == 0:
        return np.array([0.0, 0.0])

    rho = atmosphere.density(height)
    mach = calculateMach(relativeSpeed, height)
    Cd = calculateDragCoefficient(mach)

    drag = 0.5 * rho * Cd * frontalArea * speed**2
    return -drag * (relativeSpeed / speed)

def calculateAcceleration(d, v, m, T):
    gravity = calculateGravityVector(d)
    drag = calculateDragVector(v, d)
    return (T + drag) / m + gravity

def calculateSpeedOfSound(height):
    h = max(height, 0)

    if h < 11000:
        T = T0 + LAPSE * h
    else:
        T = 216.65

    return np.sqrt(GAMMA * R * T)

def calculateMach(v, height):
    speed = np.linalg.norm(v)
    a = calculateSpeedOfSound(height)
    return speed / a

def calculateDragCoefficient(M):
    Cd0 = dragCoefficient

    if M < 0.9:
        return Cd0
    else:
        rise = 0.25 / (1 + np.exp(-20 * (M - 0.9)))
        decay = np.exp(-(M - 1.05) / 0.5)
        return Cd0 + rise * decay * 2

def calculateRelativeVelocity(d, v):
    r = np.linalg.norm(d)
    r_hat = d / r
    t_hat = np.array([r_hat[1], -r_hat[0]])
    air_vel = EARTH_ANGULAR_VELOCITY * EARTH_RADIUS * t_hat
    return v - air_vel

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

def orbitalParameters():
    positionVector = position[-1]
    velocityVector = velocity[-1]
    orbitalRadius = np.linalg.norm(positionVector)
    orbitalVelocity = np.linalg.norm(velocityVector)

    mu = GRAVITATIONAL_CONSTANT * EARTH_MASS
    specificEnergy = 0.5 * orbitalVelocity ** 2 - mu / orbitalRadius
    angularMomentum = np.cross(positionVector, velocityVector)
    eccentricity = math.sqrt(1 + (2 * specificEnergy * (angularMomentum ** 2) / (mu ** 2)))
    semiMajorAxis = -mu / (2 * specificEnergy)
    semiMinorAxis = semiMajorAxis * np.sqrt(1 - eccentricity**2)
    periapsis = semiMajorAxis * (1 - eccentricity) - EARTH_RADIUS
    apoapsis  = semiMajorAxis * (1 + eccentricity) - EARTH_RADIUS

    print(f"Semi-major axis: {semiMajorAxis / 1000:.1f} km")
    print(f"Eccentricity: {eccentricity:.4f}")
    print(f"Periapsis altitude: {periapsis/1000:.1f} km")
    print(f"Apoapsis altitude: {apoapsis/1000:.1f} km")

    theta = np.linspace(0, 2 * np.pi, 1000)
    earth_x = EARTH_RADIUS * np.cos(theta)
    earth_y = EARTH_RADIUS * np.sin(theta)

    x = semiMajorAxis * np.cos(theta) - semiMajorAxis * eccentricity
    y = semiMinorAxis * np.sin(theta)
    e_vec = np.array([velocityVector[1] * angularMomentum, -velocityVector[0] * angularMomentum]) / mu - positionVector / orbitalRadius
    orbitAngle = math.atan2(e_vec[1], e_vec[0])
    x_rot = x * np.cos(orbitAngle) - y * np.sin(orbitAngle)
    y_rot = x * np.sin(orbitAngle) + y * np.cos(orbitAngle)
    
    fig, ax = plt.subplots()
    ax.plot(earth_x, earth_y)
    ax.plot(x_rot, y_rot, 'r--', label='Final orbit')
    ax.set_aspect("equal")

    trail, = ax.plot([], [], 'b-')
    dot, = ax.plot([], [], 'bo', markersize=5)

    def update(frame):
        trail.set_data(position[:frame, 0], position[:frame, 1])
        dot.set_data([position[frame, 0]], [position[frame, 1]])
        return trail, dot

    ani = FuncAnimation(fig, update, frames=range(0, len(position), 10), interval=1, blit=True)
    plt.tight_layout()
    plt.show()

for i in range(steps - 1):
    stageTime += dt

    if stageTime < currentStage["burnTime"]:
        mass[i + 1] = mass[i] - massFlowRate * dt
        thrust = calculateThrust(currentStage["seaISP"], currentStage["vacuumISP"], calculateAltitude(position[i]))
    elif currentStageIndex + 1 < len(stages):
        print(f"Stage {currentStageIndex + 1} sep at t={time[i]:.1f}s")
        mass[i] -= currentStage["dryMass"]
        currentStageIndex += 1
        stageTime = 0
        currentStage = stages[currentStageIndex]

        fuelMass = currentStage["wetMass"] - currentStage["dryMass"]
        massFlowRate = fuelMass / currentStage["burnTime"]
        thrust = calculateThrust(currentStage["seaISP"], currentStage["vacuumISP"], calculateAltitude(position[i]))
        mass[i + 1] = mass[i]
    else:
        thrust = 0
        mass[i + 1] = mass[i]

    r = np.linalg.norm(position[i])
    v = np.linalg.norm(velocity[i])

    if v >= np.sqrt(GRAVITATIONAL_CONSTANT * EARTH_MASS / r):
        thrust = 0

    if time[i] < VERTICAL_CLIMB_TIME:
        theta = np.radians(launchAngle)
    elif time[i] < VERTICAL_CLIMB_TIME + PITCH_KICK_DURATION:
        theta -= np.radians(PITCH_KICK_DEGREES / PITCH_KICK_DURATION) * dt
    else:
        relative_vel = calculateRelativeVelocity(position[i], velocity[i])
        theta = math.atan2(relative_vel[1], relative_vel[0])

    time[i + 1] = time[i] + dt
    altitude[i + 1] = calculateAltitude(position[i])
    angle[i + 1] = theta
    thrust_vector = thrust * np.array([np.cos(theta), np.sin(theta)])
    position[i + 1], velocity[i + 1] = RK4(
        position[i], velocity[i], mass[i], thrust_vector
    )

    if calculateAltitude(position[i]) < 0:
        final_step = i - 1
        break

orbitalParameters()
