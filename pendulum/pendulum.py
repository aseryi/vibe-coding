# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Pendulum parameters
L = 1.0  # length (m)
g = 9.81  # gravity (m/s^2)
theta0 = np.pi / 4  # initial angle (45 degrees)
omega0 = 0  # initial angular velocity

# Time settings
dt = 0.02
t_max = 20

# State variables
theta = theta0
omega = omega0

# Set up the figure
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 0.5)
ax.set_aspect('equal')
ax.set_title('Simple Pendulum')

line, = ax.plot([], [], 'o-', lw=2, markersize=20)

def update(frame):
    global theta, omega

    # Simple pendulum equation: d²θ/dt² = -(g/L) * sin(θ)
    alpha = -(g / L) * np.sin(theta)
    omega += alpha * dt
    theta += omega * dt

    # Calculate pendulum position
    x = L * np.sin(theta)
    y = -L * np.cos(theta)

    line.set_data([0, x], [0, y])
    return line,

ani = FuncAnimation(fig, update, frames=int(t_max / dt),
                    interval=dt * 1000, blit=True)
plt.show()
