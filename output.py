import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from constants import *
import math

def output(position, velocity):
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
