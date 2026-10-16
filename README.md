# Venus HTMC

Code for the MSc thesis *(Thesis_EricHolenstein_final.pdf)*: a 2D hydro-thermo-mechanical-chemical (HTMC) model of hydrated Venusian crust. Rock and water properties come from lookup tables precomputed with MAGEMin.

## Pipeline

```
extractor/extractor.jl  →  table_<name>.dat  →  extractor/dat_to_mat.m  →  databases/database_<name>.mat  →  HTMC.m
      (Julia, MAGEMin)                              (MATLAB)                                                  (MATLAB)
```

You only need steps 1–2 to build a new composition or grid. The six databases are already in `databases/`.

## Files

| Path | What it is |
|---|---|
| `HTMC.m` | Main model: staggered-grid finite differences with marker-in-cell, coupled solid–fluid flow, energy equation, hydration |
| `databases/database_<name>.mat` | P–T–H₂O lookup tables, one per composition (VEGA, VEGA-2, VBasalt, Venera-13, Venera-14, BSE) |
| `extractor/extractor.jl` | Runs MAGEMin (`ume` database) over a P–T–H₂O grid. It computes solid, fluid (Pitzer & Sterner 1994 EoS) and system density, enthalpy, α, β, c_p and bound water |
| `extractor/dat_to_mat.m` | Converts `table_<name>.dat` into the 3D arrays (`np × nt × nwt`) that `HTMC.m` loads |
| `extractor/plot_hydration.jl` | 3D scatter (plus a spinning GIF) of where free water exists in a database |
| `hydration_plots/` | Output of `plot_hydration.jl` for each composition |
| `file.txt` | Restart counter used by `HTMC.m` (see below) |

## 1. Generate a table (Julia)

Requirements: Julia with `MAGEMin_C, DelimitedFiles, GLMakie, Statistics, ForwardDiff, QuadGK`.

Edit the **Inputs** block in `extractor.jl`:

```julia
H₂O_wt  = LinRange(0.0, 5.0, 10)       # added H₂O (wt%)
oX      = C["BSE"]                     # composition key from the dict C
Prange, Trange = (0.0080, 1.0), (740, 1300)   # GPa, K
nP, nT  = 80, 50
generate_plots = true                  # GLMakie; set false on headless machines
```

```bash
julia -t auto extractor.jl
```

The script writes `./+op/table.dat`. Rename it to `table_<name>.dat`, using the same name as the composition key.

## 2. Convert to .mat (MATLAB)

In `dat_to_mat.m`, set `data_name = '<name>'`, then run it next to the `.dat` file. Move the resulting `database_<name>.mat` into `databases/`.

`dat_to_mat.m` reads columns **by position**, so the column order must stay the same as the header the extractor writes.

## 3. Run the model (MATLAB)

Run from the repo root:

```matlab
>> HTMC
```

### Tweakable parameters

Everything you set is in the setup block, roughly lines 1–450 of `HTMC.m`. The tables go in file order.

**Run & composition**

| Variable | Default | Meaning |
|---|---|---|
| `nname` | `'VENERA-14_0.09_0.03_300_'` (not updated to the current values) | Prefix for output files. Convention: `<db>_<XH>_<dXH>_<dT>_` |
| `Venus_composition` | `5` | 1 VEGA, 2 VEGA-2, 3 VBasalt, 4 Venera-13, 5 Venera-14, 6 BSE. Picks the database and the oxide composition `X_Venus_wt` |

**Geometry & resolution**

| Variable | Default | Meaning |
|---|---|---|
| `xsize`, `ysize` | 40 km, 25 km | Model box (y points down) |
| `Nx`, `Ny` | 81, 51 | Grid nodes (500 m spacing) |
| `Nxmc`, `Nymc` | 4, 4 | Markers per cell |
| `gy` | 8.87 | Gravity on Venus, m/s² |

**Rheology & material properties**

| Variable | Default | Meaning |
|---|---|---|
| `flow_law` | `1` | Dislocation creep: 1 plagioclase, 2 diabase, 3 mafic granulite (`A_dis`, `n_dis`, `Ea_dis`, `Va_dis`) |
| `ksolidm`, `hrsolidm`, `kphim0` | per material | Conductivity, radiogenic heat, reference permeability. Columns are material types: 1 crust, 2 anomaly, 3 unused, 4 sticky air |
| `etafluid`, `kfluid` | 1e-3 Pa·s, 0.6 W/m/K | Water viscosity and conductivity |
| `phim0`, `phimin` | 1e-3, 1e-4 | Reference porosity and minimum porosity |
| `etamin`, `etamax` | 1e16, 1e24 Pa·s | Viscosity cut-offs |
| `kphimin`, `kphimax` | 1e-20, 1e-10 m² | Permeability cut-offs |

**Boundary conditions**

| Variable | Default | Meaning |
|---|---|---|
| `tktop` | 740 K | Surface temperature |
| `ptop`, `psurf` | 67.53 bar, 95.6 bar | Pressure at the top of the sticky air (5 km above the surface) and at the surface |
| `bcleft/right/top/bottom` | `-1` | Solid velocity BC: −1 free slip, 1 no slip |
| `bcf*` | `-1` | Fluid velocity BC: −1 free slip, 1 no slip |

**Initial conditions (marker loop)**

| Setting | Default | Meaning |
|---|---|---|
| Sticky air | `ym < 5000` | Top 5 km is air (type 4) |
| Crust temperature | `tktop + ((ym-5000)/20000)*100` | Linear profile, 740 → 840 K |
| `XHtotalm` | `0.02` | Initial bulk water (mol fraction) in the crust |
| `etasolidm` | 1e21 crust, 1e17 air | Initial viscosity |
| `an_width_x/y` | 30 km, 20 km | Size of the Gaussian anomaly |
| `an_x_center`, `an_y_center` | `xsize/2`, `0.9*ysize` | Anomaly centre: horizontally centred, near the bottom |
| `an_dT_gauss` | 300 K | Peak temperature added, scaled by the Gaussian `an_G` (smooth, no sharp edge) |
| `an_dXH_gauss` | 0.04 | Peak water added, scaled by `an_G` |
| `an_G > 0.2` | | Threshold for which markers count as anomaly (type 2) |

**Timestepping & output**

| Variable | Default | Meaning |
|---|---|---|
| `dt`, `dtkoef` | 1e10 s, 1.2 | Initial timestep and its growth factor |
| `dxymax`, `DTmax`, `dphimax` | 0.1, 20 K, 0.01 | Timestep limits: marker displacement (grid steps), temperature change, porosity change |
| `dnfmax` | 100 mol/m³ | Timestep limit: max change in free water per step (set in the time loop, ~line 1656) |
| `titermax`, `pferrmax` | 200, 1e4 Pa | Max iterations and convergence tolerance |
| `savefigstep`, `savematstep` | 50, 50 | How often to save `.jpg` and `.mat` |
| `desktopyn` | 0 | Show figures (1) or not (0) |

Notes:
- **Hard-coded geometry:** the 5 km air layer, the `20000` and the `100` in the temperature profile are written as numbers, not derived from other variables. If you change `ysize` or want a different bottom temperature, edit them too. `tkbottom` is recomputed from the initial profile, so setting it has no effect.
- **Database range:** `XHtotalm + an_dXH_gauss` must stay inside the database's H₂O axis (max ≈ 0.13–0.16 depending on composition). Temperatures must stay within 740–1300 K and pressures below 1 GPa. Values outside these ranges are linearly extrapolated from the edge cell without any warning.
- **Defined but unused:** diffusion creep (`A_dif`, `m_dif`, …), `mu`, `b`, `h` and `rhm` are set but never used. Density comes from the database.

**Restarting:** `file.txt` holds the last saved timestep.
- `0` starts a new run.
- `N` reloads `<nname>N.mat` and continues from there.

Reset it to `0` before starting a new run.

Every 50 steps (and at step 5), the model saves `<nname><step>.mat` and `<nname><step>.jpg` to the working directory. It also creates a folder called `<nname>`, but nothing is written into it.

## 4. Plot free water (optional)

```bash
julia extractor/plot_hydration.jl databases/database_VBasalt.mat
```

This needs `MAT` and `Plots`. It writes `PT_freewater_3D_scatter.png` and a `.gif` to the current directory.
