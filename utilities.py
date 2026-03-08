import numpy as np
from constants import *

def calculateAltitude(d, planet):
    return np.linalg.norm(d) - planet.RADIUS

def calculateRelativeVelocity(pos, vel, planet):
    omega = np.array([0, 0, planet.ANGULAR_VELOCITY])
    air_vel = np.cross(omega, pos)
    return vel - air_vel