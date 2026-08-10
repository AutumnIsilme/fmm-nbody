# N-Body simulation using the adaptive fast multipole method

Implemented with periodic boundary conditions on [-1,1]

## Structure

### cloud.py
This file reads and writes files to load clouds of bodies for the simulation, as well as containing methods to generate distributions of bodies.

### octree.py
This file defines the octree data structure used by the adaptive method

### fast_multipole.py
This file has functions that perform timesteps on the data passed to it using the fast multipole method

## References

[1] Carrier, J., Greengard, L., & Rokhlin, V. (1988). A fast adaptive multipole algorithm for particle simulations. SIAM journal on scientific and statistical computing, 9(4), 669-686.
[2] Keyframe Codes (2026). The Fastest Gravity Algorithm You've Never Heard Of: Fast Multipole Method. https://www.youtube.com/watch?v=FhMftauQZqU, retrieved 2026-08-10

