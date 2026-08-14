# 2D Fast Multipole Method (FMM) Gravitational N-Body Simulation

## Quick Start

To run the simulation, you must compile the project using `make`, open `Viewer.html` in your browser, and then launch `./fmm-run` with any desired flags.

### 1. Build the Project
Compile the executable using `make`:

```bash
make
```

### 2. Open the Live View
Open **`Viewer.html`** in your browser to view the simulation output:

```bash
# On macOS
open Viewer.html

# On Linux
xdg-open Viewer.html
```

### 3. Run the Simulation
Execute `./fmm-run` with any custom flags you want to pass. Ensure to run it from the same directory as `Viewer.html`

The simulation will now be running in your browser

## Command-Line Interface
The CLI supports both long and short flags to customize simulation runs on the fly:

| Flag | Short | Default | Description |
| :--- | :--- | :--- | :--- |
| `--bodies <int>` | `-n` | `500` | Total number of particles ($N$) |
| `--bodies-per-box <int>` | `-b` | `25` | Maximum particles in a leaf box before splitting |
| `--epsilon <double>` | `-e` | `1e-7` | Multipole precision tolerance |
| `--steps <int>` | `-s` | `5000` | Total number of integration sub-steps |
| `--dt <double>` | — | `1e-3` | Time step delta ($dt$) |
| `--rebuild-every <int>` | `-r` | `1` | Frequency (in steps) to rebuild the Quadtree |
| `--no-boxes` | — | — | Disable background quadtree box rendering |
| `--show-boxes` | — | `enabled` | Enable background quadtree box rendering |
| `--help` | `-h` | — | Display help menu and exit |

---

## Simulation Constants (`include/SimConstants.h`)

Core physical constants and numerical parameters are defined in **`include/SimConstants.h`**. Adjust these constants to tune the numerical behavior of the simulation engine:

* **`kSofteningSquared` (`0.000225` / `0.015^2`):** Softening length scale squared ($h^2$). Passed to the **Cubic Spline Softening** kernel to smoothly cap close-range pairwise forces and prevent singular acceleration spikes.
* **`kTimestepSafetyFactor` (`0.05`):** Adaptive timestepping factor. Smaller values enforce smaller sub-steps for higher integration accuracy during high-velocity interactions.
* **`kMaxSubstepsPerFrame` (`1000`):** Safety ceiling on adaptive sub-steps per frame to prevent execution hangs during tight orbital passes.
* **`kMinSubstepDt` (`1e-7`):** The absolute floor for time-step resolution ($dt$) during sub-stepping iterations.
* **`kEscapeRadiusMultiplier` (`1.5`):** Domain boundary threshold factor. Particles moving beyond this multiple of domain widths are marked as un-bound/escaped.
