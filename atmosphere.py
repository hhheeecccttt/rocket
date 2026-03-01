import math

SEA_LEVEL_GRAVITY = 9.80665
SPECIFIC_GAS_CONSTANT = 287.05
SOLAR_FLUX_INDEX = 140
GEOMAGNETIC_INDEX = 5

#Base Height, Base Temp, Base Density, Lapse Rate 
layers = [
    [0,       288.15, 1.225,   -0.0065],   # Troposphere
    [11000,   216.65, 0.3639,   0.0],      # Lower Stratosphere
    [20000,   216.65, 0.0880,   0.001],    # Upper Stratosphere
    [32000,   228.65, 0.0132,   0.0028],   # Stratopause
    [47000,   270.65, 0.00143,  0.0],      # Lower Mesosphere
    [51000,   270.65, 0.00086, -0.0028],   # Upper Mesosphere
    [71000,   214.65, 6.4e-5,  -0.002],    # Mesopause
    [86000,   186.65, 2.0e-7,   0.0],      # Thermosphere
]

def density(height):
    for i in range(len(layers) - 1):
        baseHeight, baseTemperature, baseDensity, lapseRate = layers[i]
        h1 = layers[i + 1][0]
        if height < h1:
            break
    else:
        if (height < 600000):
            scaleHeight = (900 + 2.5 * (SOLAR_FLUX_INDEX - 70) + 1.5 * GEOMAGNETIC_INDEX) / (27 - 0.012 * (height / 1000 - 200)) * 1000
            return 2.0e-7 * math.exp(-(height - 86000) / scaleHeight)
        else:
            return 1e-10 * math.exp(-(height - 600000) / 50000)

    if lapseRate == 0:
        return baseDensity * math.exp(-SEA_LEVEL_GRAVITY * (height - baseHeight) / (SPECIFIC_GAS_CONSTANT * baseTemperature))
    else:
        T = baseTemperature + lapseRate * (height - baseHeight)
        exponent = - (SEA_LEVEL_GRAVITY / (SPECIFIC_GAS_CONSTANT * lapseRate) + 1)
        return baseDensity * math.pow(T / baseTemperature, exponent)