import numpy as np
import matplotlib.pyplot as plt

EARTH_RADIUS = 6371000
TM_DIFFERENTIAL = 15
g = 9.8
RHO = 1.2251
CD = 0.75
A = 80

wet_mass = 2214000
dry_mass = 137000

max_thrust = 34500000
thrust = 0
burn_time = 6
burn_rate = (wet_mass - dry_mass) / burn_time

dt = 0.01
t_max = 2000
steps = int(t_max / dt)
final_step = 0

time = np.zeros(steps)
velocity = np.zeros(steps)
position = np.zeros(steps)
mass = np.zeros(steps)

mass[0] = wet_mass

for i in range(steps - 2):
    time[i+1] = time[i] + dt

    if (time[i] < burn_time):
        mass[i] = wet_mass - (wet_mass - dry_mass) * (time[i] / burn_time)
        thrust = max_thrust
    else:
        mass[i] = dry_mass
        thrust = 0

    drag = 0.5 * RHO * CD * A * velocity[i] ** 2 * np.sign(velocity[i])

    acceleration = (thrust - mass[i] * g - drag) / mass[i]

    velocity[i+1] = velocity[i] + acceleration * dt
    position[i+1] = position[i] + velocity[i] * dt

    if position[i+1] < 0:
        final_step = i
        break
    
time = time[:final_step]
velocity = velocity[:final_step]
position = position[:final_step]
mass = mass[:final_step]


plt.figure(figsize=(10,6))
plt.subplot(2,1,1)
plt.plot(time, position)
plt.ylabel("Altitude (m)")

plt.subplot(2,1,2)
plt.plot(time, velocity)
plt.ylabel("Velocity (m/s)")
plt.xlabel("Time (s)")

plt.tight_layout()
plt.show()