# N-Body simulation using the adaptive fast multipole method

Implemented with periodic boundary conditions on [-1,1]

## Structure

### cloud.py
This file reads and writes files to load clouds of bodies for the simulation, as well as containing methods to generate distributions of bodies.

### octree.py
This file defines the octree data structure used by the adaptive method

### fast_multipole.py
This file has functions that perform timesteps on the data passed to it using the fast multipole method

