import numpy as np
from constants import *

def calculateAltitude(d, planet):
    return np.linalg.norm(d) - planet.RADIUS

def calculateRelativeVelocity(pos, vel, planet):
    r = np.linalg.norm(pos)
    r_hat = pos / r
    t_hat = np.array([r_hat[1], -r_hat[0]])
    air_vel = planet.ANGULAR_VELOCITY * planet.RADIUS * t_hat
    return vel - air_vel