import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from constants import *
import math

def output(position, velocity, planet):
    positionVector = position[-1]
    velocityVector = velocity[-1]
    orbitalRadius = np.linalg.norm(positionVector)
    orbitalVelocity = np.linalg.norm(velocityVector)

    mu = GRAVITATIONAL_CONSTANT * planet.MASS
    specificEnergy = 0.5 * orbitalVelocity ** 2 - mu / orbitalRadius
    angularMomentum = np.cross(positionVector, velocityVector)
<<<<<<< HEAD
    eccentricity = math.sqrt(1 + (2 * specificEnergy * (angularMomentum ** 2) / (mu ** 2)))
    semiMajorAxis = -mu / (2 * specificEnergy)
    semiMinorAxis = semiMajorAxis * np.sqrt(1 - eccentricity**2)
    periapsis = semiMajorAxis * (1 - eccentricity) - planet.RADIUS
    apoapsis  = semiMajorAxis * (1 + eccentricity) - planet.RADIUS
=======
    nodeVector = np.cross([0, 0, 1], angularMomentum)
>>>>>>> 3D-Orbitals

    

    eccentricity = max(0.0, math.sqrt(abs(1 + (2 * specificEnergy * np.linalg.norm(angularMomentum) ** 2) / mu**2)))
    eccentricityVector = np.cross(velocityVector, angularMomentum) / mu - positionVector / orbitalRadius
    semiMajorAxis = -mu / (2 * specificEnergy)
    inclination = np.degrees(np.arccos(np.clip(angularMomentum[2] / np.linalg.norm(angularMomentum), -1, 1)))
    if np.linalg.norm(nodeVector) < 1e-6:
        RAAN = 0.0
        AOP = 0.0
    else:
        RAAN = np.degrees(np.arccos(np.clip(nodeVector[0] / np.linalg.norm(nodeVector), -1, 1)))
        AOP = np.degrees(np.arccos(np.clip(np.dot(nodeVector, eccentricityVector) / (np.linalg.norm(nodeVector) * eccentricity), -1, 1)))
    trueAnomaly = np.degrees(np.arccos(np.clip(np.dot(eccentricityVector, positionVector) / (eccentricity * orbitalRadius), -1, 1)))
    
    

    print(f"Semi-major axis (a): {semiMajorAxis / 1000:.1f} km")
    print(f"Eccentricity (e): {eccentricity:.4f}")
    print(f"Inclination (i): {inclination:.2f}°")
    print(f"Right Ascension of the Ascending Node (Ω): {RAAN:.2f}°")
    print(f"Argument of Periapsis (ω): {AOP}°")
    print(f"True Anomaly (v / θ): {trueAnomaly:.2f}°")

    if eccentricity < 1.0:
        periapsis = semiMajorAxis * (1 - eccentricity) - planet.RADIUS
        apoapsis  = semiMajorAxis * (1 + eccentricity) - planet.RADIUS
        print(f"Periapsis altitude: {periapsis/1000:.1f} km")
        print(f"Apoapsis altitude: {apoapsis/1000:.1f} km")

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d') 

    R = planet.RADIUS * 2
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_zlim(-R, R)

    u, v = np.mgrid[0:2*np.pi:18j, 0:np.pi:9j]
    ex = planet.RADIUS * np.cos(u) * np.sin(v)
    ey = planet.RADIUS * np.sin(u) * np.sin(v)
    ez = planet.RADIUS * np.cos(v)
    ax.plot_wireframe(ex, ey, ez, color='blue', alpha=0.2, linewidth=0.5)

    ax.plot(position[:, 0], position[:, 1], position[:, 2],
            'r-', linewidth=1.5, label='Trajectory')
    ax.plot(*position[0], 'go', markersize=6, label='Launch')
    ax.plot(*position[-1], 'rs', markersize=6, label='Final position')

    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(planet.RADIUS * np.cos(theta),
            planet.RADIUS * np.sin(theta),
            np.zeros(200), 'b--', alpha=0.4, label='Equator')

    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('Rocket Trajectory (ECI Frame)')
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    plt.tight_layout()
    plt.show()