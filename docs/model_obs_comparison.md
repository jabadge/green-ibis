# Model vs observation comparison tools

These tools compare **observed** surface velocities with **modelled** velocities from ISSM
over Greenland glacier basins. They were developed in the ModelObsCompare (moc) repository
and merged into green-ibis. Every capability exists in both **MATLAB** and **Python**, with
the same conventions and data. Both live in the same `src/` folders.

- **Regions of interest (ROIs):** draw a box or polygon on a velocity backdrop, or derive one
  from an observation file's extent. ROIs save as `.mat`/`.json` and as ISSM `.exp` contours.
- **Spatial difference maps** (model − obs) over an ROI, with both fields regridded onto a
  common grid.
- **Time-series comparison** at a point or along a flowline, with bias, RMSE and correlation.
- **Observed variability maps:** how much a pixel's speed changes over a melt season or a
  whole record, in absolute and speed-normalised terms.
- **Flowlines from archived shapefiles:** given an ROI, take the centre flowline of the
  glacier through it and clip it at the most recently traced terminus.

Bias quantification and bias-corrected differences are included but still in early development.

## MATLAB or Python?

| | MATLAB | Python |
|---|---|---|
| setup | `setup_paths` (repo root) | `src/` on `PYTHONPATH` (see the root README) |
| example call | `load_obs_netcdf(...)` | `from data_io.obs_netcdf import load_obs_netcdf` |
| model side | reads ISSM `.mat` output directly (`load_model_output`) | reads a netCDF written by `export_model.m` (`load_model`) |

Use MATLAB when you are already in an ISSM session and want the model object itself. Use
Python for the observation side, notebooks, and anything you want to run without a MATLAB
license. A typical mixed workflow runs `export_model.m` once per model in MATLAB, then does
all plotting and analysis in Python. The exported netCDF loads in about a second, where the
original ISSM `.mat` can be several GB.

**Related green-ibis code:** `src/data_io/extract_obs_from_model.m` pulls out the observations
that were *used to calibrate* a transient run. The `load_obs_*` readers here read observation
products independently of any calibration.

## Requirements

- **Python:** `numpy scipy h5py xarray matplotlib pandas` (the core dependencies in
  `pyproject.toml`). Optional extras:
  - `shapefile` (`geopandas shapely pyogrio`), needed only for `src/flowline/flowline_shp.py`
  - `lazy` (`dask rioxarray zarr h5netcdf fsspec`), needed only for lazy/cloud reads
  - `interactive` (`ipympl jupyterlab`), for notebooks and interactive ROI drawing

  Everything else imports without the extras.
- **MATLAB:** R2018b+ recommended for `drawrectangle`/`drawpolygon`. Older releases fall back to
  `ginput`.
- **ISSM binary distribution** (optional). It provides `loadmodel`, `InterpFromMeshToGrid`,
  `InterpFromMesh2d`, `SectionValues` and `expread`/`expwrite`. Without it, only the MATLAB
  model-side functions are unavailable. The observation, ROI, variability and flowline tools
  all still work.
- **Mapping Toolbox:** needed only by the MATLAB `src/flowline/` functions (`shaperead`,
  `polyxpoly`).

## Configuring data paths

No path is hard-coded to a machine. MATLAB and Python read and write the same config file,
`~/.config/green-ibis/paths.json`, so answering once configures both. The lookup order is in
the root README. You can also pass a location explicitly, e.g.
`load_flowlines(roi, 'dir', '/data/flowlines')` in MATLAB or
`load_flowlines(roi, directory="/data/flowlines")` in Python.

```matlab
configure_paths('-set')                              % prompt for anything missing
configure_paths('obs_netcdf_dir', '/data/velocity')  % this session
configure_paths('-save')                             % ... and remember it
configure_paths()                                    % show every setting and its source
d   = data_path('model_dir');                        % resolve one (prompts once if needed)
cfg = project_config();  cfg.epsg                    % conventions: epsg, grid_res, obs_bands
```

```bash
PYTHONPATH=src python -m config.config          # show every setting and its source
PYTHONPATH=src python -m config.config --set    # prompt for anything missing
```

```python
from config.config import CONFIG, configure, save_config
configure(obs_netcdf_dir="/data/velocity")
save_config()                 # persist to ~/.config/green-ibis/paths.json
CONFIG.status()               # {'obs_netcdf_dir': ('/data/velocity', 'argument'), ...}
```

A data path might not be configured yet, so in MATLAB get it with `data_path(...)` rather than
assuming a field exists on `project_config()`.

| setting | environment variable | what to point it at |
|---|---|---|
| `issm_dir` | `ISSM_DIR` | root of the ISSM binary distribution (holds `bin/`, `lib/`) |
| `model_dir` | `GREENIBIS_MODEL_DIR` | ISSM model output `.mat` / exported `.nc` files |
| `obs_flowline_dir` | `GREENIBIS_OBS_FLOWLINE_DIR` | per-glacier observed flowline `.mat` files |
| `obs_netcdf_dir` | `GREENIBIS_OBS_NETCDF_DIR` | gridded observed-velocity netCDF files |
| `flowline_shp_dir` | `GREENIBIS_FLOWLINE_SHP_DIR` | Felikson per-glacier flowline shapefiles |
| `termini_shp` | `GREENIBIS_TERMINI_SHP` | Black & Joughin traced-terminus shapefile |
| `icefronts_nc` | `GREENIBIS_ICEFRONTS_NC` | NSIDC-0793 ice fronts/catchments netCDF (used by the shear-margin workflow) |
| `roi_dir` | `GREENIBIS_ROI_DIR` | where ROIs are saved (default: the repo's `data/`, not tracked by git) |

The datasets are not in this repository, because they are large and archived externally.

## Data assumptions

**Model output** (`model_dir/Model_*.mat`) holds saved ISSM `model` objects. The loader
handles both transient runs (`results.TransientSolution` with a `time` array and
`Vx/Vy/Vel/Surface` per step) and snapshot inversions (`results.StressbalanceSolution`).
Velocities are treated as **m/yr**. `load_model_output` extracts only the mesh and the fields
it needs. Use `'trange'` to subset transient runs.

**Observed gridded velocity** (`obs_netcdf_dir/*.nc`, e.g. `Sentinel_Subset_Upernavik.nc`):
- Layout: a 4-D `velocity(time, component, y, x)` array with bands `{vv, vx, vy, ex, ey, dT}`.
  The band order is read from the file's `component` variable when present.
- Coordinates are in EPSG:3413. `time` is stored as CF "days since …".
- The files are netCDF4 (HDF5) and are read with `h5py`, so `netCDF4`/`rioxarray` aren't needed.
- The readers match axes by size, so dimension order doesn't matter. They convert time to
  decimal years.
- **Nodata:** these files declare a useless `_FillValue` *attribute* (NaN) and keep the real
  per-band sentinels in a `_FillValue` *variable*: `-1` on the non-negative bands (`vv`, `ex`,
  `ey`, `dT`) and `-2e9` on `vx`/`vy`. Both readers mask these to NaN. Without the mask,
  15–30% of pixels read as a −1 m/yr speed and corrupt every map, difference and variance
  downstream.

**Observed flowline files** (`obs_flowline_dir/<glacier>.mat`, MATLAB v7), one glacier per file:
- Contents: ordered `(x, y)` points, along-line distance `d`, static geometry (`zb`, `zthick`),
  and velocity time series.
- The velocity series comes from the first of these layouts found: `filteredV` (+ `vti`) →
  `SentinelV` → `SentinelVmonthly` → `u2`.
- Times are **MATLAB datenum** and are converted to decimal years.
- Point/time orientation is detected automatically. For example, `u2.v` is stored time × points
  and gets transposed.

**Flowline shapefiles** (`flowline_shp_dir/glacier*.shp`) are from Felikson et al. (2020). Each
file is one glacier and holds several `Line` features identified by a `flowline` attribute.
Vertices run **from the terminus inland**, so point 1 is the downstream end. `*_iterNN.shp`
files are intermediate versions and are skipped by default.

**Terminus shapefile** (`termini_shp`) holds Black & Joughin traced terminus positions: one
`Line` per glacier per image, with a `SourceDate` (`yyyy-mm-dd`) and a `Quality_Fl`. Only
quality `0` is kept by default. Its `Glacier_ID` is **unrelated** to the flowline collection's
glacier IDs, so the two are paired spatially, never by ID.

Everything is read as-is in EPSG:3413. Nothing is reprojected.

## Model export format (MATLAB → Python)

`export_model.m` writes a neutral netCDF4 file (`Conventions = moc_py-model-1.0`) containing
`x(nv)`, `y(nv)`, `elements(ne,3)` (1-based, tagged `index_base=1`), `time(nt)` in decimal
years, `vel/vx/vy(nv,nt)`, and optionally `surface`/`bed`, with `epsg` in the attributes. On
the Python side, `load_model` (in `data_io/model_reader.py`):
- converts elements to 0-based;
- matches axes by size, so the MATLAB/HDF5 transpose quirk is handled;
- interpolates with `matplotlib.tri.LinearTriInterpolator`, the same P1 linear interpolation as
  ISSM's `InterpFromMesh2d` / `InterpFromMeshToGrid`.

See `examples/example_export_model.m`.

## Quick start

```matlab
cfg = setup_paths();

% spatial difference over an ROI
M   = load_model_output('Model_NW_HindcastRun_TransientInversion.mat');
O   = load_obs_netcdf('Sentinel_Subset_Upernavik.nc', 'bands',{'vv'});
roi = roi_from_obs(O, 'margin',2000);
D   = spatial_diff(M, O, roi, 'time',2019.0, 'obswindow',0.1);
plot_velocity_diff(D, 'roi',roi);

% flowline time series
Fobs = load_obs_flowline('jakobshavn.mat');
Fmod = model_flowline(M, Fobs);
C    = compare_timeseries(M, Fobs, 'dist',5000);
plot_velocity_timeseries(C, 'datetime',true);
```

```python
from data_io.obs_netcdf import load_obs_netcdf
from data_io.obs_flowline import load_obs_flowline
from data_io.model_reader import load_model, model_flowline
from roi.roi import roi_from_obs
from compare.compare import spatial_diff, compare_timeseries
from plotting.plot import plot_diff, plot_compare

M   = load_model("exported/Model_NW_fric1_Snapshot.nc")   # from export_model.m
O   = load_obs_netcdf("Sentinel_Subset_Upernavik.nc", bands=["vv"])
roi = roi_from_obs(O, margin=2000)
D   = spatial_diff(M, O, roi, time=2019.6, obswindow=0.1)
plot_diff(D, roi=roi)

F = load_obs_flowline("jakobshavn.mat")
C = compare_timeseries(M, F, dist=5000)
plot_compare(C)
```

The Python readers return `xarray.Dataset`s with a decimal-year `time` coordinate plus a
`datetime64` `time_dt` coordinate. Full walkthroughs are in `examples/` and `notebooks/`.

## Normalised variability maps

`velocity_variance` reduces the time axis of a gridded observation to a per-pixel map of how
much the speed *changed* over a window, divided by that pixel's own reference speed. The
absolute spread is dominated by the fast trunk simply because it is fast. The normalised field
is dimensionless, so neighbouring pixels are comparable as a **fraction of relative change**.

```matlab
V = velocity_variance(O, roi, 'trange',[2019.4 2019.8]);    % roi may be []
plot_variance(V, 'field','norm',   'roi',roi);   % fractional range [-]
plot_variance(V, 'field','spread', 'roi',roi);   % same window in m/yr
V.stats.normMean, V.stats.normMedian              % ROI-mean relative change
V2 = velocity_variance(O, roi, 'times',[2019.45 2019.75]);   % two individual slices
plot_variance(V2, 'field','normChange', 'roi',roi);          % signed, centred on 0
```

```python
from compare.variance import velocity_variance
from plotting.plot import plot_variance

V = velocity_variance(O, roi, trange=(2019.4, 2019.8))
plot_variance(V, "norm", roi=roi)
V.attrs["norm_mean"], V.attrs["norm_median"]
V2 = velocity_variance(O, roi, times=(2019.45, 2019.75))
plot_variance(V2, "norm_change")
```

The output fields are:
- `spread` and `reference` [m/yr], and `norm = spread / reference`;
- `count`, the number of frames per pixel;
- `change` and the normalised change (`normChange` in MATLAB, `norm_change` in Python), only
  when exactly two frames are used.

Options:
- `statistic`: `range` (default, peak-to-peak) / `std` / `iqr` / `mad`
- `normalize`: `mean` (default) / `median` / `first` / `max` / `none`
- `mincount` / `min_count`: frames a pixel needs (default 2)
- `minspeed` / `min_speed`: guards against huge ratios in near-stagnant ice (off by default)

## Flowlines from the shapefiles

Given an ROI, `flowline_from_roi`:
1. loops over the Felikson collection and keeps the flowlines that reach into the ROI;
2. takes the centre flowline;
3. splits it at the most recent Black & Joughin terminus trace that crosses it, and keeps the
   upstream side;
4. clips the result to the ROI.

```matlab
roi  = load_roi('sverdrups_trunk');          % or a struct, bbox, or np x 2 matrix
FL   = flowline_from_roi(roi, 'spacing',200);
FL.terminus_date                             % the trace it was clipped at
FL.d(1)                                      % 0 at the terminus; d increases inland
Fobs = sample_obs_along_flowline(O, FL);     % adds .vel (npoints x ntime), .time, .source
Fmod = model_flowline(M, FL);
plot_flowline(Fmod, Fobs, 'mode','profile', 'time',2019.6);
```

```python
from flowline.flowline_shp import flowline_from_roi, sample_obs_along_flowline
FL   = flowline_from_roi(roi, spacing=200)
Fobs = sample_obs_along_flowline(O, FL)      # (point, time), like load_obs_flowline
```

The result has the same geometry fields as the observed-flowline readers, so it works with
`model_flowline`, `compare_timeseries` and `plot_flowline`. You can also run the steps
individually (MATLAB shown; the Python names are the same):

```matlab
F = load_flowlines(roi);                 % glacier, flowline, x, y
T = load_termini(roi);                   % date, quality, x, y
c = center_flowline(F, 'glacier','a045');     % fallback ranking: ordered_flowlines
t = latest_terminus(T, 'line',[c.x c.y]);
[xu, yu] = clip_upstream(c.x, c.y, t.x, t.y);
[xc, yc] = clip_to_roi(xu, yu, roi);
FL = flowline_to_struct(xc, yc, 'name','a045_06', 'spacing',200);
```

## Lazy / cloud reads (Python)

## Reading velocities straight from NSIDC

`load_obs_netcdf` reads a whole local file into memory. For remote COG stacks, zarr datacubes
or HTTP netCDF, use `open_obs_dataarray`. It returns a lazy (dask-backed) `DataArray`
normalised to the toolkit's conventions: dims `(y, x, time)`, decimal-year `time` plus
`time_dt`, and an `epsg` attribute. Nothing is read until you slice or compute. It also accepts
an array that is already lazy, such as a GrIMP `myVelSeries.vv`:

```python
from data_io.obs_lazy import open_obs_dataarray, to_obs_dataset
da = open_obs_dataarray(myVelSeries.vv, band="vv")    # -> (y, x, time), lazy
O  = to_obs_dataset(da)                                # downstream fns compute only what they touch
da = open_obs_dataarray("https://.../cube.zarr", band="v", chunks="auto")
da = open_obs_dataarray("s3://bucket/vv_*.tif", band="vv")
```

This needs the `lazy` extra.

## Tests and parity checks

**Smoke tests** (`tests/smoke_test.m`, `python tests/smoke_test.py`) need configured data paths.
They load one of each data type (auto-picking the smallest files), print shapes and ranges, and
check that:
- netCDF bands come back as `(y, x, time)` with `vv ≈ hypot(vx, vy)`, which catches a wrong
  band order;
- flowline velocity is `(point, time)`, which catches the transposed `u2.v`;
- (MATLAB) `InterpFromMeshToGrid` output matches an independent `InterpFromMesh2d` sample,
  which catches a transposed or flipped regrid.

**Synthetic tests** need no data. Each one fabricates its own inputs and checks the results
against known answers:
- `tests/compare/test_variance_synthetic.{m,py}`: a fast-ice/slow-ice grid where the absolute
  spread is 10× larger on the fast side and the normalised spread 2× larger on the slow side.
- `tests/flowline/test_flowline_shp_synthetic.{m,py}`: the flowline pipeline on fabricated
  shapefiles. The Python version skips if geopandas is absent.
- `tests/data_io/test_model_synthetic.py`: a mesh with a known linear field. Interpolation must
  be exact (~1e-13), and a model-vs-itself difference must be ~0.

**Sverdrups reference run.** Sverdrups Glacier is where a gridded subset
(`SentinelMonthly_Sverdrups.nc`, 72 monthly frames, 2015–2020) and a redrawn flowline
(`sverdrups.mat`, 31 points, 0–7.5 km, 2015–2024) cover the same ice.
`examples/example_sverdrups_workflow.{m,py}` and `notebooks/python_sverdrups_workflow.ipynb`
run the whole toolkit on it.
- The MATLAB model section is off by default. Set `RUN_MODEL = true` to include it.
- The Python version runs the model section by default. Use `--no-model` to skip it, or
  `--save-figs` to write the figures.

The per-year table below (ROI, melt season Apr–Nov, `minspeed` 20 m/yr) was reproduced to
every digit by both languages before the merge. It is a useful regression check:

| year | frames | mean relative change | median | mean range |
|------|--------|----------------------|--------|------------|
| 2015 | 7 | 27.1% | 14.5% | 156 m/yr |
| 2016 | 7 | 39.6% | 15.7% | 170 m/yr |
| 2017 | 7 | 34.3% | 17.1% | 211 m/yr |
| 2018 | 7 | 34.8% | 15.2% | 204 m/yr |
| 2019 | 7 | 44.2% | 19.4% | 190 m/yr |
| 2020 | 7 | 38.8% | 18.5% | 183 m/yr |

## Conventions

- Coordinates are in metres, EPSG:3413. Model, observations and flowlines must share this
  projection.
- Time is in decimal years. MATLAB datenum and CF times are converted on load.
- ROIs round-trip between `.json` (Python) and `.exp` (ISSM/MATLAB), so an ROI made in one
  language works in the other and can define a green-ibis analysis domain.
