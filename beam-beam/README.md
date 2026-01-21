# Beam-Beam Collision Simulation

2D simulation of electron-positron bunch collision with beam-beam electromagnetic forces.

## Initial Prompt

> Please create the code for beam-beam simulations. The bunch of electrons need to collide head-on with the bunch of positrons. We will assume 2D case, i.e. we only consider longitudinal coordinate Z and vertical coordinate Y. Let's assume that the bunches are represented by 50 vertical slices, placed equidistantly along Z axis, and each slice has 21 particles, initially placed equidistantly along Y coordinate. Therefore, there are 50*21 particle in each bunch. Let's work in dimensionless units for the coordinates, and assume that the longitudinal distance between slices is equal to dZ=1, and let's assume that vertical separation between neighboring particles in a slice is also equal dY=0.5. Before collision, please place bunches left and right from zero coordinates, so they would just not overlap. The bunch placed on the left from zero coordinate will move, with each time step, by dZ=1 to the right, and vice versa for the bunch placed on the right. The initial vertical velocity of each particle will be assigned to be zero. When a slice of the bunch will overlap (i.e. will have the same Z coordinate) with the oncoming bunch, the slice of the oncoming bunch will create a focusing force, acting on its particles. (Note that the electrical field of the bunch is not acting on the particle of its own bunch, but only on the opposite bunch). We will use the assumption that in horizontal direction (which we ignore for particle motion) the bunches are very wide, and in this case the vertical electric field is very easy to calculate in the following manner. Let's assume that we need to calculate the electrical field in the location Y1 of a particular slice. The electrical field in the location Y1 will be equal to the sum of the number of particles of the opposite slice which have Y>Y1, minus the number of particles of the opposite slice which have Y<Y1. The electrical field, calculated at the position of each particle of the slice, will change the vertical velocity of this particle, proportionally to the electrical field at the location of the particle and to the constant coefficient DY. Let's try to make this simulation (with visualization) for beams colliding and exiting collisions. We would need to adjust the parameter DY so that the number of oscillations of a particle in the field of the opposite bunch, after the complete collision, is equal about one. My rough estimation is that we can start with DY equal around 0.01. After we will be satisfied with simulations of centered collisions, we will explore the case, when there is an initial vertical offset between bunches before the collision. Please let me know if instructions are sufficiently clear or some clarifications are needed. Thank you! Let's play with colliding beams!

## Follow-up Adjustments

1. **Sign correction**: Changed the sign of DY_KICK to make the beam-beam force focusing (toward center) rather than defocusing
2. **Repeating collisions**: Added continuous loop to repeat collisions with same initial conditions for easier observation
3. **Parameter tuning**: Reduced DY_KICK from 0.01 → 0.005 → 0.0025 to achieve approximately one oscillation during collision
4. **Vertical offset**: Added Y_OFFSET = 2.0 to simulate offset collisions (electrons at Y=-1, positrons at Y=+1)
5. **Parameter display**: Added on-screen display of DY_KICK and Y_OFFSET values

## Features

- 2D beam-beam simulation with 1050 particles per bunch (50 slices × 21 particles)
- Real-time visualization with 4 panels:
  - Main collision view (Z vs Y)
  - Vertical distribution histogram
  - Phase space plot (Y vs VY) for central slice
  - Bunch centroid evolution over time
- Configurable parameters: DY_KICK (focusing strength), Y_OFFSET (initial separation)
- Repeating collisions for extended observation

## Files

- `beam_beam.py` - Interactive simulation with visualization
- `beam_beam_record.py` - Script to record the animation as GIF
- `beam_beam.gif` - Pre-recorded animation (3 collisions with offset)

## Run

```bash
# Interactive simulation
python beam_beam.py

# Record to GIF
python beam_beam_record.py
```

## Requirements

- numpy
- matplotlib

## Physics Background

This simulation models the electromagnetic interaction between colliding particle bunches in a collider. The beam-beam force is one of the fundamental limits on luminosity in particle colliders. Key phenomena visible in the simulation:

- **Focusing effect**: Particles are attracted toward the center of the opposing bunch
- **Beam-beam tune shift**: The focusing force causes particles to oscillate (betatron oscillations)
- **Offset collisions**: When bunches are vertically separated, they experience a net deflection, leading to coherent centroid oscillations

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| N_SLICES | 50 | Number of longitudinal slices per bunch |
| N_PARTICLES | 21 | Particles per slice |
| DZ | 1.0 | Longitudinal spacing between slices |
| DY | 0.5 | Initial vertical spacing between particles |
| DY_KICK | 0.0025 | Velocity kick coefficient |
| Y_OFFSET | 2.0 | Total vertical separation between bunches |
