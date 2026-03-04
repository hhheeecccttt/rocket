import json
import math
import atmosphere
from constants import *
import output
from rocket import Rocket, Stage

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

index = 0

with open("stats.json", "r") as file:
    data = json.load(file)
    rocket_data = data[index]

rocket = Rocket.from_dict(rocket_data)

currentStageIndex = 0
stageTime = 0
time = 0

fuelMass = rocket.currentStage.wetMass - rocket.currentStage.dryMass
massFlowRate = fuelMass / rocket.currentStage.burnTime
thrust = 0

print(rocket.name)

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

mass[0] = rocket.wetMass
angle[0] = launchAngle
velocity[0] = [EARTH_ANGULAR_VELOCITY * EARTH_RADIUS, 0]
position[0] = [0, EARTH_RADIUS]

def calculateThrust(seaISP, vacuumISP, height):
    pressure_ratio = atmosphere.density(height) / 1.225
    isp = seaISP + (vacuumISP - seaISP) * (1 - pressure_ratio)
    return massFlowRate * isp * SEA_LEVEL_GRAVITY

def calculateAltitude(d):
    return np.linalg.norm(d) - EARTH_RADIUS

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
    mach = atmosphere.calculateMach(relativeSpeed, height)
    Cd = calculateDragCoefficient(mach)

    drag = 0.5 * rho * Cd * rocket.frontalArea * speed**2
    return -drag * (relativeSpeed / speed)

def calculateAcceleration(d, v, m, T):
    gravity = calculateGravityVector(d)
    drag = calculateDragVector(v, d)
    return (T + drag) / m + gravity

def calculateDragCoefficient(M):
    Cd0 = rocket.dragCoefficient

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

for i in range(steps - 1):
    rocket.stageTime += dt
    currentStage = rocket.currentStage
    massFlowRate = currentStage.massFlowRate

    if rocket.stageTime < currentStage.burnTime:
        mass[i + 1] = mass[i] - massFlowRate * dt
        thrust = calculateThrust(
            currentStage.seaISP,
            currentStage.vacuumISP,
            calculateAltitude(position[i])
        )

    elif rocket.next_stage():
        mass[i] -= currentStage.dryMass
        currentStage = rocket.currentStage
        fuelMass = currentStage.fuelMass
        massFlowRate = currentStage.massFlowRate

        thrust = calculateThrust(
            currentStage.seaISP,
            currentStage.vacuumISP,
            calculateAltitude(position[i])
        )

        mass[i + 1] = mass[i]

    else:
        thrust = 0
        mass[i + 1] = mass[i]

    r = np.linalg.norm(position[i])
    v = np.linalg.norm(velocity[i])

    if v >= np.sqrt(GRAVITATIONAL_CONSTANT * EARTH_MASS / r):
        thrust = 0

    if time[i] < rocket.VERTICAL_CLIMB_TIME:
        theta = np.radians(launchAngle)
    elif time[i] < rocket.VERTICAL_CLIMB_TIME + rocket.PITCH_KICK_DURATION:
        theta -= np.radians(rocket.PITCH_KICK_DEGREES / rocket.PITCH_KICK_DURATION) * dt
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

output.output(position, velocity)
