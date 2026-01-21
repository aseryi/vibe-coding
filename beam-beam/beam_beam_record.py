# -*- coding: utf-8 -*-
"""
Beam-Beam Collision Simulation - GIF Recording
Records the simulation as an animated GIF
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Simulation parameters
N_SLICES = 50          # Number of slices per bunch
N_PARTICLES = 21       # Number of particles per slice
DZ = 1.0               # Longitudinal spacing between slices
DY = 0.5               # Initial vertical spacing between particles in a slice
DY_KICK = 0.0025       # Velocity kick coefficient (tune for ~1 oscillation)

# Initial vertical offset between bunches (set to 0 for centered collision)
Y_OFFSET = 2.0  # electrons at -1, positrons at +1

# Derived parameters
BUNCH_LENGTH = (N_SLICES - 1) * DZ
Y_MAX = (N_PARTICLES - 1) * DY / 2

class Bunch:
    def __init__(self, z_start, direction, y_offset=0.0):
        self.direction = direction
        self.n_slices = N_SLICES
        self.n_particles = N_PARTICLES
        self.z = np.zeros((N_SLICES, N_PARTICLES))
        self.y = np.zeros((N_SLICES, N_PARTICLES))
        self.vy = np.zeros((N_SLICES, N_PARTICLES))

        for i in range(N_SLICES):
            slice_z = z_start + i * DZ
            for j in range(N_PARTICLES):
                self.z[i, j] = slice_z
                self.y[i, j] = -Y_MAX + j * DY + y_offset
                self.vy[i, j] = 0.0

    def move(self):
        self.z += self.direction * DZ

    def update_y(self):
        self.y += self.vy

    def get_slice_z(self, slice_idx):
        return self.z[slice_idx, 0]

    def get_slice_y_positions(self, slice_idx):
        return self.y[slice_idx, :]

    def apply_kick(self, slice_idx, kicks):
        self.vy[slice_idx, :] += kicks


def calculate_field(y_position, opposite_slice_y):
    n_above = np.sum(opposite_slice_y > y_position)
    n_below = np.sum(opposite_slice_y < y_position)
    return n_above - n_below


def simulate_interaction(bunch1, bunch2):
    for i in range(bunch1.n_slices):
        z1 = bunch1.get_slice_z(i)
        for j in range(bunch2.n_slices):
            z2 = bunch2.get_slice_z(j)
            if abs(z1 - z2) < 0.01:
                y1_particles = bunch1.get_slice_y_positions(i)
                y2_particles = bunch2.get_slice_y_positions(j)

                kicks1 = np.zeros(N_PARTICLES)
                for k, y1 in enumerate(y1_particles):
                    field = calculate_field(y1, y2_particles)
                    kicks1[k] = DY_KICK * field

                kicks2 = np.zeros(N_PARTICLES)
                for k, y2 in enumerate(y2_particles):
                    field = calculate_field(y2, y1_particles)
                    kicks2[k] = DY_KICK * field

                bunch1.apply_kick(i, kicks1)
                bunch2.apply_kick(j, kicks2)


# Initialize bunches
electron_z_start = -N_SLICES
positron_z_start = 1

STEPS_PER_COLLISION = int((2 + 2 * BUNCH_LENGTH) / 2) + 20
NUM_COLLISIONS = 3  # Record 3 collisions
TOTAL_FRAMES = STEPS_PER_COLLISION * NUM_COLLISIONS

electrons = None
positrons = None

def reset_bunches():
    global electrons, positrons
    electrons = Bunch(electron_z_start, direction=+1, y_offset=-Y_OFFSET/2)
    positrons = Bunch(positron_z_start, direction=-1, y_offset=+Y_OFFSET/2)

reset_bunches()

print("Beam-Beam Simulation - GIF Recording")
print(f"Particles per bunch: {N_SLICES * N_PARTICLES}")
print(f"Steps per collision: {STEPS_PER_COLLISION}")
print(f"Number of collisions: {NUM_COLLISIONS}")
print(f"Total frames: {TOTAL_FRAMES}")
print(f"DY_KICK: {DY_KICK}, Y_OFFSET: {Y_OFFSET}")

# Set up visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax_main = axes[0, 0]
ax_main.set_xlim(-60, 60)
ax_main.set_ylim(-10, 10)
ax_main.set_xlabel('Z (longitudinal)')
ax_main.set_ylabel('Y (vertical)')
ax_main.set_title('Beam-Beam Collision')
ax_main.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax_main.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

electron_scatter = ax_main.scatter([], [], c='blue', s=2, alpha=0.7, label='Electrons')
positron_scatter = ax_main.scatter([], [], c='red', s=2, alpha=0.7, label='Positrons')
ax_main.legend(loc='upper right')

ax_hist = axes[0, 1]
ax_phase = axes[1, 0]
ax_centroid = axes[1, 1]

time_history = []
electron_centroid_history = []
positron_centroid_history = []

time_text = ax_main.text(0.02, 0.95, '', transform=ax_main.transAxes, fontsize=10)
param_text = ax_main.text(0.98, 0.95, f'DY_KICK = {DY_KICK}\nY_OFFSET = {Y_OFFSET}',
                          transform=ax_main.transAxes, fontsize=10, ha='right', va='top')

frame_count = 0

def update(frame):
    global electrons, positrons, frame_count

    frame_in_cycle = frame % STEPS_PER_COLLISION

    if frame_in_cycle == 0:
        reset_bunches()
        time_history.clear()
        electron_centroid_history.clear()
        positron_centroid_history.clear()

    if frame_in_cycle > 0:
        electrons.move()
        positrons.move()
        simulate_interaction(electrons, positrons)
        electrons.update_y()
        positrons.update_y()

    e_positions = np.column_stack((electrons.z.flatten(), electrons.y.flatten()))
    p_positions = np.column_stack((positrons.z.flatten(), positrons.y.flatten()))
    electron_scatter.set_offsets(e_positions)
    positron_scatter.set_offsets(p_positions)

    ax_hist.clear()
    ax_hist.hist(electrons.y.flatten(), bins=40, range=(-10, 10), alpha=0.5, color='blue', label='e-')
    ax_hist.hist(positrons.y.flatten(), bins=40, range=(-10, 10), alpha=0.5, color='red', label='e+')
    ax_hist.set_xlim(-10, 10)
    ax_hist.set_ylim(0, 200)
    ax_hist.set_xlabel('Y position')
    ax_hist.set_ylabel('Count')
    ax_hist.set_title('Vertical Distribution')
    ax_hist.legend()

    ax_phase.clear()
    central_slice = N_SLICES // 2
    ax_phase.scatter(electrons.y[central_slice, :], electrons.vy[central_slice, :],
                     c='blue', s=20, alpha=0.7, label='e-')
    ax_phase.scatter(positrons.y[central_slice, :], positrons.vy[central_slice, :],
                     c='red', s=20, alpha=0.7, label='e+')
    ax_phase.set_xlim(-10, 10)
    ax_phase.set_ylim(-0.5, 0.5)
    ax_phase.set_xlabel('Y position')
    ax_phase.set_ylabel('VY')
    ax_phase.set_title('Phase Space (Central Slice)')
    ax_phase.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax_phase.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax_phase.legend()

    time_history.append(frame_in_cycle)
    electron_centroid_history.append(np.mean(electrons.y))
    positron_centroid_history.append(np.mean(positrons.y))

    ax_centroid.clear()
    ax_centroid.plot(time_history, electron_centroid_history, 'b-', label='e- centroid')
    ax_centroid.plot(time_history, positron_centroid_history, 'r-', label='e+ centroid')
    ax_centroid.set_xlim(0, STEPS_PER_COLLISION)
    ax_centroid.set_ylim(-2, 2)
    ax_centroid.set_xlabel('Time step')
    ax_centroid.set_ylabel('Y centroid')
    ax_centroid.set_title('Bunch Centroid Evolution')
    ax_centroid.legend()
    ax_centroid.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

    collision_num = frame // STEPS_PER_COLLISION + 1
    time_text.set_text(f'Collision #{collision_num}, Step: {frame_in_cycle}/{STEPS_PER_COLLISION}')

    frame_count += 1
    if frame_count % 50 == 0:
        print(f"Rendering frame {frame_count}/{TOTAL_FRAMES} ({100*frame_count/TOTAL_FRAMES:.1f}%)")

    return electron_scatter, positron_scatter

print("\nRendering GIF...")
ani = FuncAnimation(fig, update, frames=TOTAL_FRAMES, blit=False, repeat=False)

output_path = "/Users/seryi/Library/CloudStorage/GoogleDrive-andrei.seryi@gmail.com/My Drive/Claude/beam_beam.gif"
writer = PillowWriter(fps=20)
ani.save(output_path, writer=writer)

print(f"\nGIF saved to: {output_path}")
plt.close()
