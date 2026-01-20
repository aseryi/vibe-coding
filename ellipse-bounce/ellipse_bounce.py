# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Ellipse

# Ellipse parameters
a = 2.0  # semi-major axis (horizontal)
b = 1.5  # semi-minor axis (vertical)

# Particle initial state
x, y = 0.5, 0.3
vx, vy = 0.03, 0.02  # velocity

# Set up the figure
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(-a - 0.3, a + 0.3)
ax.set_ylim(-b - 0.3, b + 0.3)
ax.set_aspect('equal')
ax.set_title('Particle Bouncing in an Ellipse')

# Draw the ellipse boundary
ellipse = Ellipse((0, 0), 2*a, 2*b, fill=False, color='blue', linewidth=2)
ax.add_patch(ellipse)

# Particle dot
particle, = ax.plot([], [], 'ro', markersize=10)

# Trail (optional - shows path)
trail_x, trail_y = [], []
trail, = ax.plot([], [], 'r-', alpha=0.3, linewidth=1)

def is_outside_ellipse(px, py):
    """Check if point is outside or on the ellipse boundary."""
    return (px**2 / a**2 + py**2 / b**2) >= 1.0

def reflect_velocity(px, py, vx, vy):
    """Reflect velocity vector off the ellipse at point (px, py)."""
    # Normal vector to ellipse at (px, py): gradient of x^2/a^2 + y^2/b^2
    nx = 2 * px / a**2
    ny = 2 * py / b**2
    # Normalize
    n_len = np.sqrt(nx**2 + ny**2)
    nx, ny = nx / n_len, ny / n_len

    # Reflect: v' = v - 2(v.n)n
    dot = vx * nx + vy * ny
    vx_new = vx - 2 * dot * nx
    vy_new = vy - 2 * dot * ny
    return vx_new, vy_new

def update(frame):
    global x, y, vx, vy

    # Move particle
    x += vx
    y += vy

    # Check for collision with ellipse boundary
    if is_outside_ellipse(x, y):
        # Reflect velocity
        vx, vy = reflect_velocity(x, y, vx, vy)

        # Push particle back inside (find point on ellipse)
        scale = np.sqrt(x**2 / a**2 + y**2 / b**2)
        x = x / scale * 0.99
        y = y / scale * 0.99

    # Update trail
    trail_x.append(x)
    trail_y.append(y)
    # Keep trail limited
    if len(trail_x) > 500:
        trail_x.pop(0)
        trail_y.pop(0)

    particle.set_data([x], [y])
    trail.set_data(trail_x, trail_y)
    return particle, trail

ani = FuncAnimation(fig, update, frames=2000, interval=16, blit=True)
plt.show()
