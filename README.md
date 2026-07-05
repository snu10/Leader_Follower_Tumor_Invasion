# Model Description

A continuum (partial-differential-equation) model of **collective tumor-cell invasion**
driven by leader–follower dynamics. The tumor is described as a set of coupled cell-density
fields that grow, die, interconvert, and migrate over a 2-D tissue domain in response to a
**chemoattractant (CA) field** and a remodelable **extracellular-matrix (ECM) stiffness field**.

The equations are integrated with an explicit finite-volume / forward-Euler scheme on a
staggered grid, accelerated on the GPU (MATLAB `gpuArray`). The entry point is
[`main_code.m`](main_code.m).


## Domain and discretization

- **Domain:** square arena of half-width `HW = 60` (mm), grid spacing `sc = 0.1`, giving a
  `1201 × 1201` node mesh (`meshgrid`).
- **Time:** step `dt = 0.005`, total time `T = 50`, explicit forward Euler.
- **Staggered grid / finite volume.** Cell densities live at cell centers. Fluxes are built on
  the cell **faces**:
  - [`BoundaryGrad`](BoundaryGrad.m) — face-centered gradients ($\partial_x,\partial_y$ via `diff/sc`).
  - [`BoundaryDen`](BoundaryDen.m) — face-centered densities (average of the two adjacent cells).
  - [`Flux_Div`](Flux_Div.m) — divergence of a face flux back onto cell centers.
  - [`Collapse_X`](Collapse_X.m) / [`Collapse_Y`](Collapse_Y.m) — center-to-face averaging.
- **Boundary conditions:** zero-flux (no-flux / Neumann). The face arrays `BD_X`, `BD_Y` are
  initialized to zero and only their **interior** faces are ever filled, so no mass crosses the
  domain edge. The scheme is mass-conservative up to the reaction terms.


##  Outputs and analysis

### Visualization
With `VisMod = 1`, an 8-tile live figure shows $P$, $C$, ECM, chemoattractant, active leaders
(fixed and adaptive scale), the leader-suppression cue, and dead cells. Set `FileSav = 1` to
write an MP4.

### Quantitative metrics — [`StarShapeAnalyze.m`](StarShapeAnalyze.m)
Enabled by `DurRec = 1` (recorded every `Rec_dT`). For each snapshot it computes:

- **`TumMass_tot` / `TumMass_liv`** — total and living tumor mass ($\sum \rho \cdot sc^2$).
- **`TumArea`** — area of the largest connected tumor component (thresholded mask).
- **`Inv_score`** — invasion index = spatial variance of the living-mass distribution
  ($\text{var}_x+\text{var}_y$); larger = more spread out.
- **Star-shape analysis** — the tumor boundary is converted to polar coordinates about the
  center; radial **peaks (gyri = invasive arms)** and **troughs (sulci = clefts)** are detected.
- **`PeakScore` (starness)** — via [`Starness`](Starness.m), the mean of
  (arm height) / (arm base width) over all arms; a dimensionless measure of how fingered /
  star-shaped the invasion front is.

> Note: `DurRec` and `TimeLapse` default to `0` in the shipped file, and the analysis path only
> runs when `DurRec = 1`. Set it to `1` to record the metrics above during a run.


## File map

**Core loop**
- [`main_code.m`](main_code.m) — parameters, field setup, time integration, visualization.

**Field construction & operators**
- [`WeakGrad.m`](WeakGrad.m), [`StpGrad.m`](StpGrad.m) — chemoattractant profiles.
- [`CA_Weight.m`](CA_Weight.m) — CA-dependent chemotaxis weighting.
- [`CA_Trans.m`](CA_Trans.m) — time-interpolated CA field transitions.
- [`BoundaryGrad.m`](BoundaryGrad.m), [`BoundaryDen.m`](BoundaryDen.m) — staggered-grid gradients / face densities.
- [`Collapse_X.m`](Collapse_X.m), [`Collapse_Y.m`](Collapse_Y.m) — center-to-face averaging.
- [`Flux_Div.m`](Flux_Div.m) — flux divergence.
- [`radial_filter_fft.m`](radial_filter_fft.m) — FFT ring convolution for leader lateral inhibition.

**Analysis**
- [`StarShapeAnalyze.m`](StarShapeAnalyze.m) — mass, area, invasion index, star-shape metrics.
- [`Starness.m`](Starness.m) — arm-height/base-width "starness" score.


## Running

Requires MATLAB with the **Image Processing** and **Parallel Computing** (GPU) toolboxes. Open
[`main_code.m`](main_code.m), set the scenario switches (`CA_cond`, `CTX_cond`, `ECM_cond`) near
the top, and run. The `%%` section headers let you also step through it cell-by-cell in the
MATLAB editor.
