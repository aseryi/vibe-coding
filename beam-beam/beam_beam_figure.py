# -*- coding: utf-8 -*-
"""
Beam-Beam Collision - Publication Figure
Generates a single snapshot at an intermediate collision moment
"""
import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters (latest values)
N_SLICES = 50
N_PARTICLES = 21
DZ = 1.0
DY = 0.5
DY_KICK = 0.0025
Y_OFFSET = 2.0

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

electrons = Bunch(electron_z_start, direction=+1, y_offset=-Y_OFFSET/2)
positrons = Bunch(positron_z_start, direction=-1, y_offset=+Y_OFFSET/2)

# Run simulation to an intermediate moment (middle of collision)
STEPS_TO_MIDDLE = int((2 + BUNCH_LENGTH) / 2)  # When bunches are overlapping

for step in range(STEPS_TO_MIDDLE):
    electrons.move()
    positrons.move()
    simulate_interaction(electrons, positrons)
    electrons.update_y()
    positrons.update_y()

# Create publication-quality figure
fig, ax = plt.subplots(figsize=(10, 6), facecolor='white')
ax.set_facecolor('white')

# Plot particles
ax.scatter(electrons.z.flatten(), electrons.y.flatten(),
           c='blue', s=8, alpha=0.8, label='Electrons', zorder=2)
ax.scatter(positrons.z.flatten(), positrons.y.flatten(),
           c='red', s=8, alpha=0.8, label='Positrons', zorder=2)

# Configure axes
ax.set_xlim(-60, 60)
ax.set_ylim(-10, 10)
ax.set_xlabel('Z (longitudinal)', fontsize=12, color='black')
ax.set_ylabel('Y (vertical)', fontsize=12, color='black')
ax.set_title('Beam-Beam Collision', fontsize=14, color='black')

# Black axes
ax.tick_params(colors='black', labelsize=10)
for spine in ax.spines.values():
    spine.set_color('black')
    spine.set_linewidth(1)

# Add subtle reference lines
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8, zorder=1)
ax.axvline(x=0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8, zorder=1)

# Legend
ax.legend(loc='upper right', fontsize=10, frameon=True, facecolor='white', edgecolor='black')

plt.tight_layout()

# Save as PNG
output_path = "/Users/seryi/Library/CloudStorage/GoogleDrive-andrei.seryi@gmail.com/My Drive/Claude/vibe-coding/beam-beam/beam_beam_collision.png"
plt.savefig(output_path, dpi=300, facecolor='white', edgecolor='none', bbox_inches='tight')
print(f"Figure saved to: {output_path}")

plt.close()
