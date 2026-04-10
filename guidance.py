from constants import GRAVITATIONAL_CONSTANT
import numpy as np

def calculateTLI(periapsisRadius, moonOrbitalRadius, currentSMA, planet):
    mu = GRAVITATIONAL_CONSTANT * planet.MASS
    a_transfer = (periapsisRadius + moonOrbitalRadius) / 2
    v_transfer = np.sqrt(mu * (2 / periapsisRadius - 1 / a_transfer))
    v_current  = np.sqrt(mu * (2 / periapsisRadius - 1 / currentSMA))
    return v_transfer - v_current