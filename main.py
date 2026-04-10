import json
import math
import utilities
import atmosphere
from constants import *
import output
import export
import guidance
from rocket import Rocket, Stage
from planets import Earth, Moon
import numpy as np

index = 0
earth = Earth()
moon = Moon()

bodies = [earth, moon]

with open("stats.json", "r") as file:
    data = json.load(file)
    rocket_data = data[index]

rocket = Rocket.from_dict(rocket_data)

print(rocket.name)

dt = 1
t_max = 1500000
steps = int(t_max / dt)
final_step = steps

time     = [0.0]
mass     = [rocket.mass]
position = [np.array([earth.RADIUS, 0.0, 0.0])]
velocity = [np.array([0.0, earth.ANGULAR_VELOCITY * earth.RADIUS, 0.0])]

velocity[0] = np.array([0.0, earth.ANGULAR_VELOCITY * earth.RADIUS, 0.0])
position[0] = np.array([earth.RADIUS, 0.0, 0.0])

def calculateDragVector(d, v, rocket):
    height = utilities.calculateAltitude(d, earth)

    if height > 80000:
        return np.array([0.0, 0.0, 0.0])

    relativeSpeed = utilities.calculateRelativeVelocity(d, v, earth)
    speed = np.linalg.norm(relativeSpeed)
    if speed == 0:
        return np.array([0.0, 0.0, 0.0])

    rho = atmosphere.density(height)
    mach = atmosphere.calculateMach(relativeSpeed, height)
    Cd = calculateDragCoefficient(mach, rocket)

    drag = 0.5 * rho * Cd * rocket.frontalArea * speed**2
    return -drag * (relativeSpeed / speed)

def calculateGravityVector(pos, bodies):
    total = np.zeros(3)
    for body in bodies:
        positionVector = pos - body.position
        distance = np.linalg.norm(positionVector)
        g_mag = GRAVITATIONAL_CONSTANT * body.MASS / distance ** 2
        total += -g_mag * (positionVector / distance)
    return total

def calculateAcceleration(d, v, m, T, bodies, rocket):
    gravity = calculateGravityVector(d, bodies)
    drag = calculateDragVector(d, v, rocket)
    return (T + drag) / m + gravity

def calculateDragCoefficient(M, rocket):
    Cd0 = rocket.dragCoefficient

    if M < 0.9:
        return Cd0
    else:
        rise = 0.25 / (1 + np.exp(-20 * (M - 0.9)))
        decay = np.exp(-(M - 1.05) / 0.5)
        return Cd0 + rise * decay * 2

def RK4(pos, vel, m, T, bodies, rocket, dt):
    k1v = vel
    k1a = calculateAcceleration(pos, vel, m, T, bodies, rocket)

    k2v = vel + 0.5 * dt * k1a
    k2a = calculateAcceleration(pos + 0.5 * dt * k1v, k2v, m, T, bodies, rocket)

    k3v = vel + 0.5 * dt * k2a
    k3a = calculateAcceleration(pos + 0.5 * dt * k2v, k3v, m, T, bodies, rocket)

    k4v = vel + dt * k3a
    k4a = calculateAcceleration(pos + dt * k3v, k4v, m, T, bodies, rocket)

    newPos = pos + (dt / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)
    newVel = vel + (dt / 6) * (k1a + 2 * k2a + 2 * k3a + k4a)

    return newPos, newVel

for i in range(steps - 1):
    height = utilities.calculateAltitude(position[i], earth)
    thrust = rocket.update(dt, height)
    if (i > 650):
        thrust = 0

    r = np.linalg.norm(position[i])
    v = np.linalg.norm(velocity[i])

    #if v >= np.sqrt(GRAVITATIONAL_CONSTANT * earth.MASS / r):
    #    thrust = 0

    if time[i] < rocket.VERTICAL_CLIMB_TIME:
        direction = position[i] / r
    elif time[i] < rocket.VERTICAL_CLIMB_TIME + rocket.PITCH_KICK_DURATION:
        t_norm = (time[i] - rocket.VERTICAL_CLIMB_TIME) / rocket.PITCH_KICK_DURATION
        pitch_angle = np.radians(90 - rocket.PITCH_KICK_DEGREES * t_norm)
        radial = position[i] / r
        az = np.radians(rocket.launchAzimuth)
        north = np.array([0, 0, 1])
        north = north - np.dot(north, radial) * radial
        north /= np.linalg.norm(north)
        east_vec = np.array([0, 1, 0])
        east_vec = east_vec - np.dot(east_vec, radial) * radial
        east_vec /= np.linalg.norm(east_vec)
        east = np.sin(az) * east_vec + np.cos(az) * north
        east = east - np.dot(east, radial) * radial
        east /= np.linalg.norm(east)
        direction = np.sin(pitch_angle) * radial + np.cos(pitch_angle) * east
    else:
        v_rel = velocity[i] - velocity[0]
        if np.linalg.norm(v_rel) < 1.0:
            direction = position[i] / r
        else:
            direction = v_rel / np.linalg.norm(v_rel)

    thrust_vector = thrust * direction
    moon.update(time[i])

    newPos, newVel = RK4(position[i], velocity[i], mass[i], thrust_vector, bodies, rocket, dt)

    position.append(newPos)
    velocity.append(newVel)
    time.append(time[-1] + dt)
    mass.append(rocket.mass)


    if height < 0:
        final_step = i - 1
        break

position = np.array(position)
velocity = np.array(velocity)

semiMajorAxis, eccentricity = output.output(position, velocity, earth, bodies)

periapsisRadius = semiMajorAxis * (1 -  eccentricity)
dv_tli = guidance.calculateTLI(periapsisRadius, moon.ORBITAL_RADIUS,  semiMajorAxis, earth)
print(f"TLI delta-v required: {dv_tli:.1f} m/s")

export.exportTrajectory(position, velocity, earth)