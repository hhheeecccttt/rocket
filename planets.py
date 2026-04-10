import numpy as np


class Earth():
    NAME = "Earth"
    RADIUS = 6371000
    MASS = 6e24
    ANGULAR_VELOCITY = 7.2921159e-5
    SEA_LEVEL_GRAVITY = 9.8
    SEA_LEVEL_AIR_DENSITY = 1.2251
    position = np.array([0.0, 0.0, 0.0])

class Moon():
    NAME = "Moon"
    RADIUS = 1737400
    ORBITAL_RADIUS = 384400000
    ANGULAR_VELOCITY = 2 * 3.14159265 / (27.3 * 24 * 3600)
    MASS = -7.342e22
    SEA_LEVEL_GRAVITY = 1.625
    SEA_LEVEL_AIR_DENSITY = 0

    def __init__(self, initial_angle=89.5):
        self.initial_angle = initial_angle
        self.position = np.zeros(3)

    def update(self, t):
        angle = self.initial_angle + self.ANGULAR_VELOCITY * t
        self.position = np.array([
            self.ORBITAL_RADIUS * np.cos(angle),
            self.ORBITAL_RADIUS * np.sin(angle),
            0.0
        ])