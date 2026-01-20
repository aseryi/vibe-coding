# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Pendulum parameters
L = 1.0  # length (m)
g = 9.81  # gravity (m/s^2)
theta0 = np.pi / 4  # initial angle (45 degrees)
omega0 = 0  # initial angular velocity

# Time settings
dt = 0.02
t_max = 30  # 30 seconds of animation

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

# Add time display
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=10)

frame_count = 0

def update(frame):
    global theta, omega, frame_count

    # Simple pendulum equation: d²θ/dt² = -(g/L) * sin(θ)
    alpha = -(g / L) * np.sin(theta)
    omega += alpha * dt
    theta += omega * dt

    # Calculate pendulum position
    x = L * np.sin(theta)
    y = -L * np.cos(theta)

    line.set_data([0, x], [0, y])

    # Update time display
    current_time = frame * dt
    time_text.set_text(f't = {current_time:.1f}s')

    frame_count += 1
    if frame_count % 100 == 0:
        print(f"Rendering frame {frame_count}/{total_frames} ({100*frame_count/total_frames:.1f}%)")

    return line, time_text

total_frames = int(t_max / dt)
print(f"Recording pendulum animation...")
print(f"Duration: {t_max} seconds")
print(f"Total frames: {total_frames}")
print(f"FPS: {1/dt:.0f}")

ani = FuncAnimation(fig, update, frames=total_frames, blit=True, repeat=False)

# Save as GIF
output_path = "/Users/seryi/Library/CloudStorage/GoogleDrive-andrei.seryi@gmail.com/My Drive/Claude/pendulum.gif"
print(f"\nSaving to: {output_path}")
writer = PillowWriter(fps=int(1/dt))
ani.save(output_path, writer=writer)

print("\nDone!")
plt.close()
