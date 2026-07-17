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


## Data availability

The metrics above, recorded for **282 simulated conditions**, are deposited as
[`RawData.mat`](RawData.mat) (1.2 MB) and can be replotted with
[`DataInvestigate.m`](DataInvestigate.m) without re-running any simulation.

### Reproducing the figures

From this folder, in MATLAB:

```matlab
DataInvestigate
```

Two switches at the top of [`DataInvestigate.m`](DataInvestigate.m) select the panel:

- **`ViewMode`** — which comparison to draw:
  1. Result vs. ECM density (fixed CA and CTX, final time)
  2. Heatmap over ECM density and time (fixed CA and CTX)
  3. Result vs. time (fixed CA, CTX, and ECM)
  4. Result vs. CA configuration (fixed ECM)
- **`CmpRes`** — which quantity to draw: `1` total tumor mass, `2` live tumor mass,
  `3` invasion score, `4` tumor area, `5` peak score.

Which *conditions* are drawn is set by `tg_cds` (ViewModes 1–3) and `tg_ecm` (ViewMode 4);
both are commented in place and meant to be edited. The script also builds a `WholeResults`
table (`282 × 11`) summarizing every run — see its header for the column definitions.

Unlike the model itself, the analysis script needs **no toolboxes and no GPU**.

### Data format

`RawData.mat` holds a single `282 × 4` cell array `DataCell`, one row per simulated condition:

| Column | Variable | Description |
| --- | --- | --- |
| 1 | `ParamTag` | 4-character condition code, `'ABCD'` (see below) |
| 2 | `DurRecMat` | `6 × 101` time record at $t = 0, 0.5, \dots, 50$: rows are time, `TumMass_tot`, `TumMass_liv`, `Inv_score`, `TumArea`, `PeakScore` |
| 3 | `Gyrus_data` | `2 × nGyri`: rows are angle $\theta$ and radius |
| 4 | `Sulcus_data` | `2 × nSulci`: rows are angle $\theta$ and radius |

#### Condition code `'ABCD'`

The code mirrors the scenario switches `CA_cond`, `CTX_cond`, and `ECM_cond` in
[`main_code.m`](main_code.m).

Character **A** — chemoattractant (CA) configuration:

1. high average, mild gradient (HA & LG)
2. medium average, medium gradient (MA & MG)
3. low average, steep gradient (LA & HG)
4. HA & LG → LA & HG over time (anti-angiogenesis)
5. LA & HG → HA & LG over time (rapid angiogenesis)

Character **B** — chemotaxis (CTX) strategy, i.e. the weighting applied by
[`CA_Weight.m`](CA_Weight.m):

1. CA-positive (higher CA enhances chemotaxis)
2. normal (gradient-based)
3. CA-negative (higher CA reduces chemotaxis)

Characters **CD** — ECM density index; ECM density = `0.1 × index`.

**Coverage.** Configurations `A = 1..3` were run across the full ECM sweep, index `01`–`30`
(density 0.1–3.0), for every CTX strategy — 270 runs. The time-varying angiogenesis cases
`A = 4, 5` were run only at ECM index `10` and `20` — 12 runs. Total 282.

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

**Deposited data & figures**
- [`RawData.mat`](RawData.mat) — recorded metrics for all 282 conditions (1.2 MB).
- [`DataInvestigate.m`](DataInvestigate.m) — replots the result figures from `RawData.mat`.


## Running

**The model.** Requires MATLAB with the **Image Processing** and **Parallel Computing** (GPU)
toolboxes. Open [`main_code.m`](main_code.m), set the scenario switches (`CA_cond`, `CTX_cond`,
`ECM_cond`) near the top, and run. The `%%` section headers let you also step through it
cell-by-cell in the MATLAB editor.

**The figures.** To inspect the published results without re-running the model, run
[`DataInvestigate.m`](DataInvestigate.m) against the deposited
[`RawData.mat`](RawData.mat) — no toolboxes or GPU required. See
[Data availability](#data-availability).
