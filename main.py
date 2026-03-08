import json
import math
import utilities
import atmosphere
from constants import *
import output
import export
from rocket import Rocket, Stage
from planets import Earth
import numpy as np

index = 0
earth = Earth()

with open("stats.json", "r") as file:
    data = json.load(file)
    rocket_data = data[index]

rocket = Rocket.from_dict(rocket_data)

print(rocket.name)

launchAngle = 90
theta = np.radians(launchAngle)

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

mass[0] = rocket.mass
angle[0] = launchAngle
velocity[0] = [earth.ANGULAR_VELOCITY * earth.RADIUS, 0]
position[0] = [0, earth.RADIUS]

def calculateGravityVector(pos, planet):
    r = np.linalg.norm(pos)
    g_mag = GRAVITATIONAL_CONSTANT * planet.MASS / r**2
    return -g_mag * (pos / r)

def calculateDragVector(d, v):
    height = utilities.calculateAltitude(d, earth)

    if height > 80000:
        return np.array([0.0, 0.0])

    relativeSpeed = utilities.calculateRelativeVelocity(position[i], velocity[i], earth)
    speed = np.linalg.norm(relativeSpeed)
    if speed == 0:
        return np.array([0.0, 0.0])

    rho = atmosphere.density(height)
    mach = atmosphere.calculateMach(relativeSpeed, height)
    Cd = calculateDragCoefficient(mach)

    drag = 0.5 * rho * Cd * rocket.frontalArea * speed**2
    return -drag * (relativeSpeed / speed)

def calculateAcceleration(d, v, m, T, planet):
    gravity = calculateGravityVector(d, planet)
    drag = calculateDragVector(d, v)
    return (T + drag) / m + gravity

def calculateDragCoefficient(M):
    Cd0 = rocket.dragCoefficient

    if M < 0.9:
        return Cd0
    else:
        rise = 0.25 / (1 + np.exp(-20 * (M - 0.9)))
        decay = np.exp(-(M - 1.05) / 0.5)
        return Cd0 + rise * decay * 2

def RK4(pos, vel, m, T, planet):
    k1v = vel
    k1a = calculateAcceleration(pos, vel, m, T, planet)

    k2v = vel + 0.5 * dt * k1a
    k2a = calculateAcceleration(pos + 0.5 * dt * k1v, k2v, m, T, planet)

    k3v = vel + 0.5 * dt * k2a
    k3a = calculateAcceleration(pos + 0.5 * dt * k2v, k3v, m, T, planet)

    k4v = vel + dt * k3a
    k4a = calculateAcceleration(pos + dt * k3v, k4v, m, T, planet)

    newPos = pos + (dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)
    newVel = vel + (dt / 6) * (k1a + 2 * k2a + 2 * k3a + k4a)

    return newPos, newVel

for i in range(steps - 1):
    height = utilities.calculateAltitude(position[i], earth)
    thrust = rocket.update(dt, height)
    mass[i + 1] = rocket.mass

    r = np.linalg.norm(position[i])
    v = np.linalg.norm(velocity[i])

    if v >= np.sqrt(GRAVITATIONAL_CONSTANT * earth.MASS / r):
        thrust = 0

    if time[i] < rocket.VERTICAL_CLIMB_TIME:
        theta = np.radians(launchAngle)
    elif time[i] < rocket.VERTICAL_CLIMB_TIME + rocket.PITCH_KICK_DURATION:
        theta -= np.radians(rocket.PITCH_KICK_DEGREES / rocket.PITCH_KICK_DURATION) * dt
    else:
        relative_vel = utilities.calculateRelativeVelocity(position[i], velocity[i], earth)
        theta = math.atan2(relative_vel[1], relative_vel[0])

    time[i + 1] = time[i] + dt
    altitude[i + 1] = height
    angle[i + 1] = theta
    thrust_vector = thrust * np.array([np.cos(theta), np.sin(theta)])
    position[i + 1], velocity[i + 1] = RK4(position[i], velocity[i], mass[i], thrust_vector, earth)

    if height < 0:
        final_step = i - 1
        break

output.output(position, velocity, earth)
export.exportTrajectory(position, velocity, earth)