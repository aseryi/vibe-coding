# -*- coding: utf-8 -*-
"""
Beam-Beam Collision Simulation
2D simulation of electron-positron bunch collision with beam-beam forces
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle

# Simulation parameters
N_SLICES = 50          # Number of slices per bunch
N_PARTICLES = 21       # Number of particles per slice
DZ = 1.0               # Longitudinal spacing between slices
DY = 0.5               # Initial vertical spacing between particles in a slice
DY_KICK = 0.0025       # Velocity kick coefficient (tune for ~1 oscillation)

# Initial vertical offset between bunches (set to 0 for centered collision)
# Total separation is Y_OFFSET: electrons get -Y_OFFSET/2, positrons get +Y_OFFSET/2
Y_OFFSET = 2.0  # electrons at -1, positrons at +1

# Derived parameters
BUNCH_LENGTH = (N_SLICES - 1) * DZ  # Total length of each bunch
Y_MAX = (N_PARTICLES - 1) * DY / 2  # Max Y extent: ±5 for 21 particles with dY=0.5

class Bunch:
    """Represents a bunch of particles"""
    def __init__(self, z_start, direction, y_offset=0.0):
        """
        Initialize bunch
        z_start: Z position of first slice
        direction: +1 for moving right, -1 for moving left
        y_offset: vertical offset of bunch center
        """
        self.direction = direction
        self.n_slices = N_SLICES
        self.n_particles = N_PARTICLES

        # Initialize particle positions
        # Shape: (N_SLICES, N_PARTICLES) for both z and y
        self.z = np.zeros((N_SLICES, N_PARTICLES))
        self.y = np.zeros((N_SLICES, N_PARTICLES))
        self.vy = np.zeros((N_SLICES, N_PARTICLES))  # Vertical velocity

        # Set initial positions
        for i in range(N_SLICES):
            slice_z = z_start + i * DZ
            for j in range(N_PARTICLES):
                self.z[i, j] = slice_z
                self.y[i, j] = -Y_MAX + j * DY + y_offset
                self.vy[i, j] = 0.0

    def move(self):
        """Move bunch by one step in Z direction"""
        self.z += self.direction * DZ

    def update_y(self):
        """Update Y positions based on velocities"""
        self.y += self.vy

    def get_slice_z(self, slice_idx):
        """Get Z position of a slice (all particles in slice have same Z)"""
        return self.z[slice_idx, 0]

    def get_slice_y_positions(self, slice_idx):
        """Get Y positions of all particles in a slice"""
        return self.y[slice_idx, :]

    def apply_kick(self, slice_idx, kicks):
        """Apply velocity kicks to particles in a slice"""
        self.vy[slice_idx, :] += kicks


def calculate_field(y_position, opposite_slice_y):
    """
    Calculate electric field at y_position due to opposite slice
    E = (number of particles with Y > y_position) - (number with Y < y_position)
    """
    n_above = np.sum(opposite_slice_y > y_position)
    n_below = np.sum(opposite_slice_y < y_position)
    return n_above - n_below


def simulate_interaction(bunch1, bunch2):
    """
    Check for overlapping slices and apply beam-beam kicks
    """
    for i in range(bunch1.n_slices):
        z1 = bunch1.get_slice_z(i)

        for j in range(bunch2.n_slices):
            z2 = bunch2.get_slice_z(j)

            # Check if slices overlap (within small tolerance)
            if abs(z1 - z2) < 0.01:
                # Get Y positions of particles in each slice
                y1_particles = bunch1.get_slice_y_positions(i)
                y2_particles = bunch2.get_slice_y_positions(j)

                # Calculate kicks for bunch1 particles due to bunch2 slice
                kicks1 = np.zeros(N_PARTICLES)
                for k, y1 in enumerate(y1_particles):
                    field = calculate_field(y1, y2_particles)
                    kicks1[k] = DY_KICK * field  # Positive for focusing

                # Calculate kicks for bunch2 particles due to bunch1 slice
                kicks2 = np.zeros(N_PARTICLES)
                for k, y2 in enumerate(y2_particles):
                    field = calculate_field(y2, y1_particles)
                    kicks2[k] = DY_KICK * field  # Positive for focusing

                # Apply kicks
                bunch1.apply_kick(i, kicks1)
                bunch2.apply_kick(j, kicks2)


# Initialize bunches
# Electron bunch on the left, moving right (+1)
# Place so that rightmost slice is at Z = -1 (just before overlap)
electron_z_start = -N_SLICES  # Slices at Z = -50, -49, ..., -1
electrons = Bunch(electron_z_start, direction=+1, y_offset=-Y_OFFSET/2)

# Positron bunch on the right, moving left (-1)
# Place so that leftmost slice is at Z = +1
positron_z_start = 1  # Slices at Z = 1, 2, ..., 50
positrons = Bunch(positron_z_start, direction=-1, y_offset=+Y_OFFSET/2)

# Calculate steps for one collision
# Bunches need to pass completely through each other
# Initial gap is 2 (from Z=-1 to Z=+1), plus 2 bunch lengths
# With relative velocity of 2*DZ per step, need (2 + 2*BUNCH_LENGTH) / 2 steps
STEPS_PER_COLLISION = int((2 + 2 * BUNCH_LENGTH) / 2) + 20  # Extra steps to see exit
TOTAL_STEPS = STEPS_PER_COLLISION * 100  # Repeat many times

def reset_bunches():
    """Reset bunches to initial conditions"""
    global electrons, positrons
    electrons = Bunch(electron_z_start, direction=+1, y_offset=-Y_OFFSET/2)
    positrons = Bunch(positron_z_start, direction=-1, y_offset=+Y_OFFSET/2)

print("Beam-Beam Simulation")
print(f"Particles per bunch: {N_SLICES * N_PARTICLES}")
print(f"Bunch length: {BUNCH_LENGTH}")
print(f"Steps per collision: {STEPS_PER_COLLISION}")
print(f"DY_KICK coefficient: {DY_KICK}")
print(f"Y offset: {Y_OFFSET}")
print("(Simulation will repeat continuously)")

# Set up visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Main collision view
ax_main = axes[0, 0]
ax_main.set_xlim(-60, 60)
ax_main.set_ylim(-10, 10)
ax_main.set_xlabel('Z (longitudinal)')
ax_main.set_ylabel('Y (vertical)')
ax_main.set_title('Beam-Beam Collision')
ax_main.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax_main.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

# Scatter plots for particles
electron_scatter = ax_main.scatter([], [], c='blue', s=2, alpha=0.7, label='Electrons')
positron_scatter = ax_main.scatter([], [], c='red', s=2, alpha=0.7, label='Positrons')
ax_main.legend(loc='upper right')

# Y distribution plot (histogram)
ax_hist = axes[0, 1]
ax_hist.set_xlim(-10, 10)
ax_hist.set_ylim(0, 200)
ax_hist.set_xlabel('Y position')
ax_hist.set_ylabel('Count')
ax_hist.set_title('Vertical Distribution')

# Phase space plot (Y vs VY) for central slice
ax_phase = axes[1, 0]
ax_phase.set_xlim(-10, 10)
ax_phase.set_ylim(-0.5, 0.5)
ax_phase.set_xlabel('Y position')
ax_phase.set_ylabel('VY (vertical velocity)')
ax_phase.set_title('Phase Space (Central Slice)')

# Y centroid evolution
ax_centroid = axes[1, 1]
ax_centroid.set_xlim(0, TOTAL_STEPS)
ax_centroid.set_ylim(-2, 2)
ax_centroid.set_xlabel('Time step')
ax_centroid.set_ylabel('Y centroid')
ax_centroid.set_title('Bunch Centroid Evolution')

# Storage for centroid history
time_history = []
electron_centroid_history = []
positron_centroid_history = []

# Time display
time_text = ax_main.text(0.02, 0.95, '', transform=ax_main.transAxes, fontsize=10)
# Parameter display
param_text = ax_main.text(0.98, 0.95, f'DY_KICK = {DY_KICK}\nY_OFFSET = {Y_OFFSET}',
                          transform=ax_main.transAxes, fontsize=10, ha='right', va='top')

def init():
    electron_scatter.set_offsets(np.empty((0, 2)))
    positron_scatter.set_offsets(np.empty((0, 2)))
    return electron_scatter, positron_scatter

def update(frame):
    global electrons, positrons

    # Check if we need to reset for a new collision
    frame_in_cycle = frame % STEPS_PER_COLLISION

    if frame_in_cycle == 0:
        reset_bunches()
        # Clear history for new cycle
        time_history.clear()
        electron_centroid_history.clear()
        positron_centroid_history.clear()

    # Simulation step
    if frame_in_cycle > 0:
        # Move bunches
        electrons.move()
        positrons.move()

        # Calculate and apply beam-beam interaction
        simulate_interaction(electrons, positrons)

        # Update Y positions based on velocities
        electrons.update_y()
        positrons.update_y()

    # Update main scatter plot
    e_positions = np.column_stack((electrons.z.flatten(), electrons.y.flatten()))
    p_positions = np.column_stack((positrons.z.flatten(), positrons.y.flatten()))
    electron_scatter.set_offsets(e_positions)
    positron_scatter.set_offsets(p_positions)

    # Update histogram
    ax_hist.clear()
    ax_hist.hist(electrons.y.flatten(), bins=40, range=(-10, 10), alpha=0.5, color='blue', label='e-')
    ax_hist.hist(positrons.y.flatten(), bins=40, range=(-10, 10), alpha=0.5, color='red', label='e+')
    ax_hist.set_xlim(-10, 10)
    ax_hist.set_ylim(0, 200)
    ax_hist.set_xlabel('Y position')
    ax_hist.set_ylabel('Count')
    ax_hist.set_title('Vertical Distribution')
    ax_hist.legend()

    # Update phase space plot (central slice)
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

    # Track centroids
    time_history.append(frame_in_cycle)
    electron_centroid_history.append(np.mean(electrons.y))
    positron_centroid_history.append(np.mean(positrons.y))

    # Update centroid plot
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

    # Update time display
    collision_num = frame // STEPS_PER_COLLISION + 1
    time_text.set_text(f'Collision #{collision_num}, Step: {frame_in_cycle}/{STEPS_PER_COLLISION}')

    return electron_scatter, positron_scatter

# Create animation
ani = FuncAnimation(fig, update, frames=TOTAL_STEPS, init_func=init,
                    interval=50, blit=False, repeat=False)

plt.tight_layout()
plt.show()
